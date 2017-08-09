#---------------------------------------------------------------------------------------------------------------------------
#Calculate continuous ranked probability score (crps)
#-------------------------------------------------------------------------------------------------------------------
rownames(models_err1) <- rep(c('clam','Bacon','Bchron','OxCal'),5)
# 
crps <- sapply(colnames(models_err1), function(x) {
  sapply(1:20, function(y) {
    xx <- substring(x,1,nchar(x)-1)                   # column names are holzmaar1 asf now we only want the lake name
    if((rownames(models_err1)[y] == 'Bacon') | (rownames(models_err1)[y] == 'clam') | (rownames(models_err1)[y] == 'Bchron')) { # analysis carried out for models that are not OxCal
      fa <- ceiling(y/4)  # determine the number of radiocarbon dates depending on the row
      pred.depth <- seq(min(lakes[[xx]][,'Depth']),max(rad.carb[fa,x][[1]][,'depth']),2) # damit man nicht die Extrapolationen vergleicht, generate depths at which reconstructions models and truth are compared, we only use every second cm because OxCal did only predict every second cm for Elk Lake 
      real.depth <- pred.depth %in% unique(lakes[[xx]][,'Depth'])   # check which of the potential depths are in the real depths in the core
      ag.lakes <- aggregate(lakes[[xx]][,'Age'], by = as.data.frame(lakes[[xx]][,'Depth']), FUN = mean) # find the mean age of all varves assigned to one depth
      colnames(ag.lakes) <- c('Depth','Age')
      ag.lakes <- ag.lakes[ag.lakes[,'Depth']%in% pred.depth[real.depth],]  # make sure we only use data at depth that are available in reality
      
      
      #This differential selection of depths is needed as we want to use the same depths in clam Bacon Bchron and OxCal
      #OxCal at maximum returns the results at 1000 depths per core      
      if(xx == 'elk') {
        if(rownames(models_err1)[y] == 'clam')  chron <- round(models_err1[y,x][[1]][seq(1,nrow(models_err1[y,x][[1]]),4),])
        if(rownames(models_err1)[y] == 'Bacon') chron <- t(round(models_err1[y,x][[1]]))[seq(1,ncol(models_err1[y,x][[1]]),2),round(seq(1,nrow(models_err1[y,x][[1]]),length.out = 1000))]      
        if(rownames(models_err1)[y] == 'Bchron') chron <- t(matrix(unlist(models_err1[y,x]),nrow=1000))[seq(1,(length(unlist(models_err1[y,x][[1]]))/1000),4),]
      }  
      
      if(xx %in% c('holzmaar','iceberg')) {
        if(rownames(models_err1)[y]== 'clam') chron <- round(models_err1[y,x][[1]][seq(1,nrow(models_err1[y,x][[1]]),4),])
        if(rownames(models_err1)[y]== 'Bacon') chron <- t(round(models_err1[y,x][[1]]))[seq(1,ncol(models_err1[y,x][[1]])),round(seq(1,nrow(models_err1[y,x][[1]]),length.out = 1000))]
        if(rownames(models_err1)[y]== 'Bchron') chron <-  t(matrix(unlist(models_err1[y,x]),nrow=1000))[seq(1,(length(unlist(models_err1[y,x][[1]]))/1000),4),]
      }
      
      
      ens.tot <- chron[which(real.depth == TRUE),] 
      real.age <- ag.lakes[,'Age']
      
      # second part of the term
      
      crps <- sapply(1:nrow(ens.tot), function(x) {
       ens <- ens.tot[x,]
        true.age <- real.age[x]
        K <- length(ens)
        ((1/K *sum(outer(ens,true.age,FUN = function(x,y) abs(x-y))))# first term in equation 14 of Tipton et al. (2015)
        - (1/(2*K^2)*sum(outer(ens,ens,FUN = function(x,y) abs(x-y))))) # second term in equation 14 of Tipton et al. (2015)
      })
      #mcrps <- avgCRPS(ens.tot,real.age)
      mcrps <- mean(crps)
    }
    
    
    
    
    
    #-------------------------------------------------------------------------------
    #oxcal Murks
    
    if(rownames(models_err1)[y] == 'OxCal') {
      js <- models_err1[y,x][[1]]
      z<-grep("data.z", js)#, extended=F)
      js1<-js[z]
      psq.depths<-as.vector(sapply(js1,function(b){
        bb<-unlist(strsplit(b,""))
        as.numeric(substring(b,which(bb=="=")+1, nchar(b)-1))
      }))
      # psq.depths<-unique(psq.depths[-(1:65)])
      
      
      
      sta<-grep("start", js)#, extended=F)
      js1<-js[sta]
      psq.sta<-as.vector(sapply(js1,function(b){
        bb<-unlist(strsplit(b,""))
        as.numeric(substring(b,which(bb=="=")+1, nchar(b)-1))
      }))
      
      ll <- which(is.na(psq.sta)== TRUE)
      psq.sta1 <- psq.sta[-(1:ll)]
      psq.depths<-unique(psq.depths[-(1:ll)])
      
      
      pro<-     grep("posterior.prob=", js)#, extended=F)
      js1<-js[pro]
      psq.pro<-as.vector(sapply(js1,function(b){
        bb<-unlist(strsplit(b,""))
        cc <- substring(b,which(bb=="=")+2,nchar(b)-2)
        as.numeric(unlist(strsplit(cc,",")))
      }))
      
      
      oxcal.elk.31 <- sapply((ll+1):length(psq.sta1), function(x){
        post.prob <- psq.pro[[x]] / sum(psq.pro[[x]])
        year <- seq(from=round(abs(psq.sta1[x]-1950)), by = -5, length.out =length(psq.pro[[x]]))
        list(post.prob=post.prob,year = year)
      })
      
     
      ag.lakes <- aggregate(lakes[[xx]][,'Age'], by = as.data.frame(lakes[[xx]][,'Depth']), FUN = mean)
      colnames(ag.lakes) <- c('Depth','Age')
      int.lakes <- approx(ag.lakes[,'Depth'],ag.lakes[,'Age'],xout = psq.depths[psq.depths%in%ag.lakes$Depth])
      age.of.int <-   int.lakes$y
      
      
      #---------------------------------------
      #crps schätzen
      # ensemble generieren
      if((xx == 'elk')|(xx == 'holzmaar')){
      ens.tot <- sapply((1:ncol(oxcal.elk.31))[psq.depths%in%ag.lakes$Depth],function(aa) {
        sort(sample(oxcal.elk.31[,aa]$year,1000,prob =oxcal.elk.31[,aa]$post.prob,replace = TRUE))
      })
      }
      if(xx == 'iceberg'){
          ens.tot <- sapply((1:ncol(oxcal.elk.31))[psq.depths%in%ag.lakes$Depth],function(aa) {
          sort(sample(oxcal.elk.31[,aa]$year,1000,prob =oxcal.elk.31[,aa]$post.prob,replace = TRUE))
        })
      ens.tot <- ens.tot[,seq(1,ncol(ens.tot),2)]   
      age.of.int <- age.of.int[seq(1,length(age.of.int),2)]
      }
      
      
      #----------------------------------------
      crps <- sapply(1:ncol(ens.tot), function(x) {
        ens <- ens.tot[,x]
       true.age <- age.of.int[x]
        K <- length(ens)
        ((1/K *sum(outer(ens,true.age,FUN = function(x,y) abs(x-y))))# first term in equation 14 of Tipton et al. (2015)
        - (1/(2*K^2)*sum(outer(ens,ens,FUN = function(x,y) abs(x-y))))) # second term in equation 14 of Tipton et al. (2015)
      })
      #mcrps <- avgCRPS(t(ens.tot),age.of.int)
      mcrps <- mean(crps)
    }
    mcrps
  })
})



crps <- round(crps,2)

nam.crps <- rep(c('clam','Bacon','Bchron','OxCal'),5)
crps.tot<- crps
rownames(crps.tot) <- paste(rep(c('clam','Bacon','Bchron','OxCal'),5),rep(1:5,each = 4), sep ='')
colnames(crps.tot) <- rep(c('elk1','holzmaar1','iceberg1'),6)



rownames(crps.tot) <- nam.crps 

rn1 <- sapply(colnames(crps.tot),function(x) substr(x,1,nchar(x)-1))

crps2 <- crps.tot
rownames(crps2) <- paste(rownames(crps.tot),rep(1:5,each=4),sep='')
colnames(crps2) <- colnames(crps.tot) 
#----------------------------------------------
#generate legend
ncal <- c(6,11,21,31,41)
leg <- sapply(1:5, function(z) {
  paste('n = ',ncal[z],',',sep='')
})

legtest <- as.matrix(unlist(strsplit(leg,',')))


#---------------------------
#plot CRPS
jpeg(paste(fig.stor,'Fi5_crps.jpg',sep=''), width = 12, height = 12,units='in',res = 600)
layout(rbind(1:3,4:6,7:9),widths = c(3,3,1),heights=rep(3,3))
par(oma = c(3,3,2,3),mar=c(0.5,0.5,0.5,0.5),cex = 1.25) # mfrow=c(3,2)


plot(c(0.1,0.25,0.5,1,2,3),crps2['OxCal1',colnames(crps2) %in% 'elk1'],type='b', ylim =c(0,80), axes = F,ylab ='',xlab='')#, main ='OxCal Elk Lake')
sapply(2:5, function(z) {
  points(c(0.1,0.25,0.5,1,2,3),crps2[paste('OxCal',z,sep=''),colnames(crps2) %in% 'elk1'],type='b',pch = z, col =z)
})
axis(2)
box()
legend('bottomright',legend =legtest, col = 1:5, lwd = 2,cex = 0.9)
mtext(cex = 1.1,side = 2, "CRPS", line = 2.2, font = 2)
mtext(cex = 1.1,side = 3, 'OxCal', line = 0.5, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'a)',font = 2)



#---------
# Bacon
plot(c(5,10,15,20),crps2['Bacon1',colnames(crps2) %in% 'elk1'][2:5],type='b', ylim =c(0,80), axes = F,ylab ='',xlab='')#, main ='Bacon Elk Lake')
sapply(2:5, function(z) {
  points(c(5,10,15,20),crps2[paste('Bacon',z,sep=''),colnames(crps2) %in% 'elk1'][2:5],type='b',pch = z, col =z)
})
mtext (side =3, line = -1.5,adj = 1, 'b)',font = 2)
box()
#
mtext(cex = 1.1,side = 3, 'Bacon', line = 0.5, font = 2)

#-----------
# Clam 
plot(rep(1.2,5),apply(crps2[nam.crps %in%"clam",colnames(crps2) %in% 'elk1'],1,median), col = 1:5,pch = 16,cex =2,xlim=c(1,2),ylim=c(0,80),axes = FALSE)
#
box()
mtext(cex = 1.1,side = 3, 'Clam and Bchron', line = 0.5, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'c)',font = 2)
#------
# Bchron
points(seq(1.5,1.9,length.out = 5),apply(crps2[nam.crps%in%"Bchron",colnames(crps2) %in% 'elk1'],1,median), col = 1:5,pch = 15,cex =1.5,xlim=c(1,2),ylim=c(0,1))
axis(4)
#
box()
mtext(cex = 1.1,side = 4, "CRPS", line = 2.2, font = 2)

#----------------------------------------------------------------------------------------------------------------
#---------------
#holzmaar
plot(c(0.1,0.25,0.5,1,2,3),crps2['OxCal1',colnames(crps2) %in% 'holzmaar1'],type='b', ylim =c(0,160), axes = F,ylab ='',xlab='')#, main ='OxCal holzmaar Lake')
sapply(2:5, function(z) {
  points(c(0.1,0.25,0.5,1,2,3),crps2[paste('OxCal',z,sep=''),colnames(crps2) %in% 'holzmaar1'],type='b',pch = z, col =z)
})
#legend('bottomleft',legend =legtest, col = 1:5, lwd = 1)
axis(2)
box()
mtext(cex = 1.1,side = 2, "CRPS", line = 2.2, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'd)',font = 2)



#---------
# Bacon
plot(c(5,10,15,20),crps2['Bacon1',colnames(crps2) %in% 'holzmaar1'][2:5],type='b', ylim =c(0,160), axes = F,ylab ='',xlab='')#, main ='Bacon holzmaar Lake')
sapply(2:5, function(z) {
  points(c(5,10,15,20),crps2[paste('Bacon',z,sep=''),colnames(crps2) %in% 'holzmaar1'][2:5],type='b',pch = z, col =z)
})
#legend('bottomleft',legend =legtest, col = 1:5, lwd = 1)
mtext (side =3, line = -1.5,adj = 1, 'e)',font = 2)
box()

#-----------
# Clam 
plot(rep(1.2,5),apply(crps2[nam.crps%in%"clam",colnames(crps2) %in% 'holzmaar1'],1,median), col = 1:5,pch = 16,cex =2,xlim=c(1,2),ylim=c(0,160),axes = FALSE)

box()
#mtext(cex = 1.1,side = 3, 'Clam', line = 0.5, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'f)',font = 2)
#------
# Bchron
points(seq(1.5,1.9,length.out = 5),apply(crps2[nam.crps%in%"Bchron",colnames(crps2) %in% 'holzmaar1'],1,median), col = 1:5,pch = 15,cex =1.5,xlim=c(1,2),ylim=c(0,1))
axis(4)
box()
mtext(cex = 1.1,side = 4, "CRPS", line = 2.2, font = 2)

#-----------------------------------------------
#------------------#------------------------------------------------------------------------------------------
#iceberg
plot(c(0.1,0.25,0.5,1,2,3),crps2['OxCal1',colnames(crps2) %in% 'iceberg1'],type='b', ylim =c(0,30), axes = F,ylab ='',xlab='')#, main ='OxCal iceberg Lake')
sapply(2:5, function(z) {
  points(c(0.1,0.25,0.5,1,2,3),crps2[paste('OxCal',z,sep=''),colnames(crps2) %in% 'iceberg1'],type='b',pch = z, col =z)
})
axis(2)
axis(1, at = c(0.1,0.5))
axis(1, at = c(0.25),labels = FALSE)
axis(1, at = c(1,2,3))
box()
mtext(cex = 1.1,side =1, 'k',line = 2.2, font =2)
mtext(cex = 1.1,side = 2, "CRPS", line = 2.2, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'g)',font = 2)



#---------
# Bacon
plot(c(5,10,15,20),crps2['Bacon1',colnames(crps2) %in% 'iceberg1'][2:5],type='b', ylim =c(0,30), axes = F,ylab ='',xlab='')#, main ='Bacon iceberg Lake')
sapply(2:5, function(z) {
  points(c(5,10,15,20),crps2[paste('Bacon',z,sep=''),colnames(crps2) %in% 'iceberg1'][2:5],type='b',pch = z, col =z)
})
mtext (side =3, line = -1.5,adj = 1, 'h)',font = 2)
#axis(2)
axis(1, at = c(5,10,15,20,30))
box()
mtext(cex = 1.1,side =1, 'segment length [cm]',line = 2.2, font =2)

#-----------
#------
# Bchron
plot(seq(1.5,1.9,length.out = 5),apply(crps2[nam.crps%in%"Bchron",colnames(crps2) %in% 'iceberg1'],1,median), col = 1:5,pch = 15,cex =1.5,xlim=c(1,2),ylim=c(0,30),axes=FALSE)
axis(4)
box()
mtext(cex = 1.1,side = 4, "CRPS", line = 2.2, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'i)',font = 2)
dev.off()

