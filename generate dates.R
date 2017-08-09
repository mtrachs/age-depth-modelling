#-------------------------------------------------------------------------------------------------------------------
#Code Trachsel and Telford, All age-depths models are wrong, but are getting better. The Holocene  
#
#Code in this file generates radiocarbon dates and runs age-depth models
#-------------------------------------------------------------------------------------------------------------------
#First we need to define where several files and modelling routines can be found or stored
cal.curve.path <- '~/age-depth/clam/IntCal13.14C'
elk.path <- '~/age-depth/elk.lake.txt'
holzmaar.path <- '~/age-depth/holzmaar.txt'
iceberg.path <- '~/age-depth/iceberg.txt'
oxcal.path <- '~/OxCal/'
oxcal.ex <- '~/OxCal/OxCal/bin/OxCalLinux'
clam.loc <- '~/age-depth/clam/'
bacon.loc <- '~/age-depth/LinBacon_2.2/'
#save.loc <- 'W://test'
fig.stor <- '~/age-depth/'
install.packages('Bchron')
library(Bchron)




rad.car <- read.table(cal.curve.path)
colnames(rad.car) <- c('cal.year','rad.year','uncertainty')
elk <- read.table(elk.path,header = TRUE)
elk[,'Depth'] <- round(elk[,'Depth']/10,1)
holzmaar <- read.table(holzmaar.path,dec = ',',header = TRUE)
iceberg <-  read.table(iceberg.path,header = TRUE,dec =',')
iceberg[,'Depth'] <- cumsum(iceberg[,'depth'])/10
iceberg [,'Age'] <-  2000- iceberg[,'Year']
iceberg <- iceberg[,c('Age','Depth')]

lakes <- list(elk = elk,holzmaar = holzmaar, iceberg = iceberg)


#get lake data at 0.5 cm increment
for(i in 1:3){
for(x in 1:nrow(lakes[[i]])){  # Daten auf 0.5 runden
if(lakes[[i]][x,'Depth'] %% 1 <= 0.2) lakes[[i]][x,'Depth'] <- lakes[[i]][x,'Depth'] - lakes[[i]][x,'Depth'] %% 1
if((lakes[[i]][x,'Depth'] %% 1 > 0.2) & (lakes[[i]][x,'Depth'] %% 1 <= 0.7)) lakes[[i]][x,'Depth'] <- lakes[[i]][x,'Depth'] - lakes[[i]][x,'Depth'] %% 1 + 0.5
if(lakes[[i]][x,'Depth'] %% 1 > 0.7) lakes[[i]][x,'Depth'] <- lakes[[i]][x,'Depth'] - lakes[[i]][x,'Depth'] %% 1 + 1
}
}






#function to simulate radiocarbon dates
# core: core for which  radiocarbon dates are simulated
# n: number of radiocarbon dates
# uncertainty: uncertainty of the radiocarbon dates
gen.rad.car <- function(core,n,uncertainty) 
{
xmin <- min(core[,'Depth'])
xmax <- max(core[,'Depth'])
depth <- seq(xmin,xmax, length.out = (n+2))[1:(n+1)]
depth1 <- unique(core[core[,'Depth'] %in% depth,'Depth'])
nf <- setdiff(depth,depth1)
if(length(nf)> 0) {
depth.add <- sapply(1:length(nf), function(x) {
core[max(which(core[,'Depth'] - nf[x] <0)),'Depth']
})
depth2 <- sort(c(depth1,depth.add))
}
else {depth2 <- depth1}
ages <- core[core[,'Depth'] %in% depth2,'Age']
ages <- c(ages[diff(ages) > 1], ages[length(ages)])
ages <- unique(ages - ages%%5)
rad.car.age <- round(rad.car[rad.car[,'cal.year'] %in% ages,'rad.year'] + rnorm(length(ages),0,uncertainty))
depth <- depth2
cbind(depth, rad.car.age)
}
#-------------------------------------------------
# 

nam <- names(lakes)
# define the number of radiocarbon dates to be generated
number <- c(5,seq(10,40,10))

#-------------------------------------------------------------------------------
oxcal.n <- c(0.5,0.5,2,1,2)
names(oxcal.n) <- nam

k <- c(0.1,0.25,0.5,1,2,3)
seg.len <-  c(5,5,10,15,20,30) #first entry was changed to five as 1 took very long to fit. 


#-----------------------------------------------------------------------------------
#generate radiocarbon dates
rad.carb <- sapply(nam,function(x){
  sapply(number,function(y){
    gen.rad.car(lakes[[x]],y,30) 
  })
})

colnames(rad.carb) <- nam
rownames(rad.carb) <- as.character(c(6,seq(11,41,10))) 
#-------------------------------------------------------------------------------
#code takes between a couple days and weeks to run :-)
# check that Bacon converged: age-depth model along with trace plots (visual inspection) is stored in the bacon folder 
#if Bacon didn't converge needs new settings, rerun the models manually 
#this code was modified, so that it runs with the current versions of the age-depth modelling routines
models_err <-  lapply(1:length(k) , function(z) { # 
  sapply(nam, function(x) {
    sapply(1:nrow(rad.carb), function(y) {#
    
    
      
      test <- as.data.frame(rad.carb[y,x][[1]])
      n <- nrow(test)
      
      #--------------------------------------------------------------------------------
      #OxCal works for Linux systems, PC needs different settings 
      command1 <- rev(paste('R_Date("',gsub("^\\s+|\\s+$","", 1:nrow(test)),'",',test[,'rad.car.age'],',',rep(30,nrow(test)),'){z=',test[,'depth'],';};',sep=''))
      command<-paste(command1, collapse="")
      zz <- oxcal.n[x]
      command <- paste('P_Sequence("bb",',k[z],',',zz,'){Boundary();',command,'Boundary();};',sep='', collapse="")
      #works for Linux systems, PC needs different settings 
      writeLines(command, con = paste('~/OxCal/oxcal_',x,'_',y,'_',k[z],'.input',sep=''))
      system(paste(oxcal.ex,paste('~/OxCal/oxcal_',x,'_',y,'_',k[z],'.input',sep='')),wait = TRUE)#if using wait = TRUE code gets stuck if OxCal doesn't converge 
      chron.oxcal <- readLines(paste(oxcal.path,paste('oxcal_',x,'_',y,'_',k[z],'.js',sep=''),sep=''))# every user should check this path#system(paste(firefox.loc,'-url', command),wait = FALSE)
      
      #-------------------------------------------------------------------------------
      #Clam preparation
      prep.clam <- data.frame(rad.car.age = test$rad.car.age, cal.age = rep(NA,n),uncertainty =  rep(50,n),reservoir = rep(0,n), depth = test$depth, thickness = rep(1,n))
      write.csv(prep.clam,paste(clam.loc,'Cores/',x,'/',x,'.csv',sep =''))
      ag.lakes <- aggregate(x = lakes[[x]][,'Age'], by  = lakes[[x]]['Depth'], FUN = mean)  #mittleres alter pro Tiefe wird generiert
      colnames(ag.lakes) <- c('Depth','Age')
      ag.lakes[,'Age'] <- round(ag.lakes[,'Age'])
      cor.seq <- seq((min(lakes[[x]][,'Depth'])),(max(lakes[[x]][,'Depth'])),0.5) %in% ag.lakes[,'Depth'] #an welchen Tiefen hat holzmaar Daten?      ist f?r den Spezialfall Holzmaar so umgebaut
      pred.depth <- seq((min(lakes[[x]][,'Depth'])),(max(lakes[[x]][,'Depth'])),0.5)[cor.seq]  # an diesen Tiefen k?nnen Vergleiche angestellt werden
      
      #_-------------------------------------------------------------------------------
      # clam
      setwd(clam.loc)
      source('clam.R')
      clam(x, type=4, smooth = 0.4,  storedat = TRUE, depthseq =seq(min(lakes[[x]][,'Depth']),max(lakes[[x]][,'Depth']),0.5),mixed.effect = FALSE)
      chron.clam <- round(chron)
      #---------------------------------------------------------------------------------
      #Bacon
      setwd(bacon.loc)
      source('Bacon.R')
      prep.bacon <- data.frame(rad.car.age = test$rad.car.age,uncertainty =  rep(30,n), depth = test$depth)
      
      write.csv(prep.bacon,paste(bacon.loc,'Cores/',x,'/',x,'.csv',sep =''))
      
      laenge <- (max(lakes[[x]][,'Depth']) - min(lakes[[x]][,'Depth']))
      zeit <- (max(lakes[[x]][,'Age']) - min(lakes[[x]][,'Age']))
      acc.mean <- round(zeit / laenge)                 # mem.strength=15, mem.mean=0.7
      
      Bacon(x,thick = seg.len[z], acc.shape=1.4, acc.mean=acc.mean, mem.strength=4, mem.mean=0.7,d.min = min(lakes[[x]][,'Depth']),d.max = max(lakes[[x]][,'Depth']),d.by = 20, ask = FALSE, suggest = FALSE)#d.by = ifelse(x ==elk,1,2)dmin und dmax setzten
      lines(lakes[[x]][,c('Depth','Age')],col='blue', lwd = 2)
      
      chron.bacon <- sapply(seq(min(lakes[[x]][,'Depth']),max(lakes[[x]][,'Depth']),ifelse(x =='elk',1,2)), function(x) { 
        Bacon.Age.d(x)
      })
      #for Holzmaar all models were rerun using different settings
      #-------------------------------------------------------------------------------
      #Bchron
      myrun <- Bchronology(test$rad.car.age,rep(30,nrow(test)),test$depth,rep(1,nrow(test)),predictPositions = seq(min(test$depth), max(test$depth), 0.5),iterations = 110000, burn = 10000, thin = 100)
      plot(myrun)
      chron.bchron <- myrun$thetaPredict
      
      #-------------------------------------------------------------------------------
      chrons <- list(clam = chron.clam, bacon = chron.bacon,bchron = chron.bchron,oxcal = chron.oxcal)
    })
  })
#in case you have a lovely IT-department that reboots the system every now and again, you might want to store an intermediate output :-)  
#  save.image(paste(save.loc,x,k[z],'.RData', sep =''))
})
#-------------------------------------------------------------------------------
models_err1 <- cbind(models_err[[1]],models_err[[2]],models_err[[3]],models_err[[4]],models_err[[5]],models_err[[6]])
colnames(models_err1) <- paste(rep(c('elk','holzmaar','iceberg'),6),rep(1:6,each = 3),sep='')
rownames(models_err1) <- paste(rep(c('clam','Bacon','Bchron','OxCal'),5),rep(number+1,each = 4),sep='')
#alternative
#rownames(models_err1) <- rep(c('clam','Bacon','Bchron','OxCal'),5)
