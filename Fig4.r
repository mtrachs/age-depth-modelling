#--------------------------------------------------------------------------------------------------------------------
#Calculate coverage rates and produce figure 4
#--------------------------------------------------------------------------------------------------------------------

colnames(models_err1) <- paste('',rep(c('elk','holzmaar','iceberg'),6),'',rep(1:6,each = 3),sep='')
colnames(rad.carb) <- colnames(models_err1)



coverage <- sapply(colnames(models_err1), function(x) {   # loop going through the columns of the data containing age-depth models
  sapply(1:20, function(y) {                            # loop going through the rows of the data containing age-depth models
      xx <- substring(x,1,nchar(x)-1)                   # column names are elk1 asf now we only want the lake name
      if((rownames(models_err1)[y] == 'Bacon') | (rownames(models_err1)[y] == 'clam') | (rownames(models_err1)[y] == 'Bchron')) { # analysis carried out for models that are not OxCal
      fa <- ceiling(y/4)  # determine the number of radiocarbon dates depending on the row
      pred.depth <- seq(min(lakes[[xx]][,'Depth']),max(rad.carb[fa,x][[1]][,'depth']),2) # damit man nicht die Extrapolationen vergleicht, generate depths at which reconstructions models and truth are compared, we only use every second cm because OxCal did only predict every second cm for Elk Lake 
      real.depth <- pred.depth %in% unique(lakes[[xx]][,'Depth'])   # check which of the potential depths are in the real depths in the core
      ag.lakes <- aggregate(lakes[[xx]][,'Age'], by = as.data.frame(lakes[[xx]][,'Depth']), FUN = mean) # find the mean age of all varves assigned to one depth
      colnames(ag.lakes) <- c('Depth','Age')
      ag.lakes <- ag.lakes[ag.lakes[,'Depth']%in% pred.depth[real.depth],]  # make sure we only use data at depth that are available in reality


      #This differential selection of depths is needed as we want to use the (approximately) same depths in clam Bacon Bchron and OxCal, OxCal starting at the bottom of the core
      #and the three other routines starting at the top, we do not always get the very same depths (but within 1 cm of each other)
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

cdfs <- lapply(which(real.depth == TRUE), function(z) {   # function for all real.depths,prdocue the cdfs (distribution functions)
Fn <- ecdf((chron)[z,])                                   # cdf of data
positions <- knots(Fn)                                    # positions of steps of cdf 
values <- Fn(knots(Fn))                                   # values at the steps 
cbind(positions,values)                                   # steps and values are returned
})

test <- sapply(1:sum(real.depth),function(z) {
minimum <- cdfs[[z]][1,'positions']                        # find the minimum value in the ensemble
maximum <- max(cdfs[[z]][,'positions'])                    # find the maximum value in the ensemble
if(ag.lakes[z,'Age']  %in%  cdfs[[z]][,'positions'])  {    # if the real age is directly found in the cdf 
pos <- cdfs[[z]][,'positions'] == ag.lakes[z,'Age']        # which is the position of this year
perc <- mean(cdfs[[z]][pos,'values'])                      # find the cdf value at this position
}
if((minimum < ag.lakes[z,'Age']) & (ag.lakes[z,'Age'] < maximum)) {  # if the real age is between the minimum and maxium of the cdf
diffpos <- abs(cdfs[[z]][,'positions'] - ag.lakes[z,'Age'])          # generate the absolute difference between positions in the ecdf and the true age
perc <- mean(cdfs[[z]][diffpos == min(diffpos),'values'])            # find the minimum difference, if two years have the same distance take the average of the cdf value of the two ages
}
else {
perc <- ifelse(ag.lakes[z,'Age'] < minimum, 0,1)                     # if the true age is outside the model, check if it is below or above and assign 0 or 1 accordingly
}
})
#zero <- test == 0
#one <- test == 1

#test1 <- test[(zero ==FALSE) & (one == FALSE)]

hist.bin <- hist(1-unlist(test),freq = FALSE, main = x, xlab = 'Percentile',breaks=seq(0,1,0.01), ylim =c(0,6))
lines(density(1-unlist(test), from = 0, to = 1, adjust= 0.5),col='blue')
box()
#res <- c(sum(zero),hist.bin$density,sum(one))
test <- 1-test
histo <- hist.bin$density
test
#list(perc = errors, histo = histo, dens = dens)
}





#-------------------------------------------------------------------------------
#Code for OxCal

if(rownames(models_err1)[y] == 'OxCal') {
js <- models_err1[y,x][[1]]
    z<-grep("data.z", js)#, extended=F)
  js1<-js[z]
  #find depths
  psq.depths<-as.vector(sapply(js1,function(b){
    bb<-unlist(strsplit(b,""))
    as.numeric(substring(b,which(bb=="=")+1, nchar(b)-1))
  }))
 # psq.depths<-unique(psq.depths[-(1:65)])


#find starting years of the pdf
  sta<-grep("start", js)#, extended=F)
  js1<-js[sta]
  psq.sta<-as.vector(sapply(js1,function(b){
    bb<-unlist(strsplit(b,""))
    as.numeric(substring(b,which(bb=="=")+1, nchar(b)-1))
  }))

  #find values that are going to be used. 
ll <- which(is.na(psq.sta)== TRUE)
#starting year of the PDF produced by P-sequence (PSQ) 
psq.sta1 <- psq.sta[-(1:ll)]
psq.depths<-unique(psq.depths[-(1:ll)])

#find the entire posterior probability
    pro<-     grep("posterior.prob=", js)#, extended=F)
  js1<-js[pro]
  psq.pro<-as.vector(sapply(js1,function(b){
    bb<-unlist(strsplit(b,""))
    cc <- substring(b,which(bb=="=")+2,nchar(b)-2)
    as.numeric(unlist(strsplit(cc,",")))
  }))


oxcal.elk.31 <- sapply((ll+1):length(psq.sta1), function(x){
        post.prob <- psq.pro[[x]] / sum(psq.pro[[x]])
        year <- seq(from=round(abs(psq.sta1[x]-1950)), by = -5, length.out =length(psq.pro[[x]]))#values of PDF ar given at 5 year increment
list(post.prob=post.prob,year = year)
})


ag.lakes <- aggregate(lakes[[xx]][,'Age'], by = as.data.frame(lakes[[xx]][,'Depth']), FUN = mean)
colnames(ag.lakes) <- c('Depth','Age')
int.lakes <- approx(ag.lakes[,'Depth'],ag.lakes[,'Age'],xout = psq.depths[psq.depths%in%ag.lakes$Depth])
age.of.int <-   int.lakes$y
oxcal.elk.31 <- oxcal.elk.31[,psq.depths%in%ag.lakes$Depth]

#find probablities and according years
probs <- sapply((1:(length(oxcal.elk.31)/2)), function(z) {
  #find year closest to the true year
  zahl <- which(abs(oxcal.elk.31[[2*z]] -  age.of.int[z]) == min(abs(oxcal.elk.31[[2*z]] -  age.of.int[z])))
  #find probability mass below the true value
  prob <- min(cumsum(oxcal.elk.31[[2*z-1]])[zahl])#min used zahl is not one index but more than one index
})


hist.bin <- hist(unlist(probs),freq = FALSE, main = x, xlab = 'Percentile',breaks=seq(0,1,0.01), ylim =c(0,6))
lines(density(unlist(probs), from = 0, to = 1, adjust = 0.5),col='blue')
box()
test <- unlist(probs)
if(xx =='iceberg')
  {
  test <- test[seq(1,length(test),2)]
  }
}
test
})
})


#-------------------------------------------------------------------------------
# proportion of varves between 16% and 83%
g.of.fit3 <- sapply(colnames(models_err1), function(x) {
        sapply(1:20, function(y) {
              dens <- ((coverage[y,x][[1]] > 0.166) & (coverage[y,x][[1]] < 0.833))
              sum(dens)/length(dens)
              })
})


rownames(g.of.fit3) <- paste(rep(rownames(models_err1)[1:4],5),rep(1:5,each = 4), sep ='')
colnames(g.of.fit3) <- rep(colnames(models_err1)[1:3],6)
#-------------------------------------------------------------------------------
aa <- substring(rownames(g.of.fit3),1,nchar(rownames(g.of.fit3))-1)
#find clam models
gof.clam.el <- matrix(ncol=5,unlist(g.of.fit3[aa%in%'clam',colnames(g.of.fit3)%in%'elk1'][,2:6]))
gof.clam.ho <- matrix(ncol=5,unlist(g.of.fit3[aa%in%'clam',colnames(g.of.fit3)%in%'holzmaar1'][,2:6]))
gof.clam.ic <- matrix(ncol=5,unlist(g.of.fit3[aa%in%'clam',colnames(g.of.fit3)%in%'iceberg1'][,2:6]))

#find bchron models
gof.bchron.el <- matrix(ncol=5,unlist(g.of.fit3[aa%in%'Bchron',colnames(g.of.fit3)%in%'elk1'][,2:6]))
gof.bchron.ho <- matrix(ncol=5,unlist(g.of.fit3[aa%in%'Bchron',colnames(g.of.fit3)%in%'holzmaar1'][,2:6]))
gof.bchron.ic <- matrix(ncol=5,unlist(g.of.fit3[aa%in%'Bchron',colnames(g.of.fit3)%in%'iceberg1'][,2:6]))

#median of clam and bchron models 
cl.el <- apply(gof.clam.el,1,median)
cl.ho <- apply(gof.clam.ho,1,median)
cl.ic <- apply(gof.clam.ic,1,median)
bc.el <- apply(gof.bchron.el,1,median)
bc.ho <- apply(gof.bchron.ho,1,median)
bc.ic <- apply(gof.bchron.ic,1,median)



#-------------------------------------------------------------------------------
#make figure 4
jpeg(paste(fig.stor,'Fig4_coverage.jpg',sep=''), width = 12, height = 12,units='in',res = 600)
layout(rbind(1:3,4:6,7:9),widths = c(3,3,1),heights=rep(3,3))
par(oma = c(3,3,2,3),mar=c(0.5,0.5,0.5,0.5),cex = 1.25) # mfrow=c(3,2)


#-------------------------------------------------------------------------------
# Elk Lake
#-----
# OxCal
plot(c(0.1,0.25,0.5,1,2,3),g.of.fit3['OxCal1',colnames(g.of.fit3) %in% 'elk1'],type='b', ylim =c(0,1), axes = F,ylab ='',xlab='')#, main ='OxCal Elk Lake')
sapply(2:5, function(z) {
points(c(0.1,0.25,0.5,1,2,3),g.of.fit3[paste('OxCal',z,sep=''),colnames(g.of.fit3) %in% 'elk1'],type='b',pch = z, col =z)
})
#legend('bottomleft',legend =legtest, col = 1:5, lwd = 1)
axis(2)
#axis(1, at = c(0.1,0.25,0.5,1,2,3))
box()
#mtext(cex = 1.1,side =1, 'k',line = 2.2, font =2)
mtext(cex = 1.1,side = 2, 'Proportion', line = 2.2, font = 2)
mtext(cex = 1.1,side = 3, 'OxCal', line = 0.5, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'a)',font = 2)
abline( h = 0.66,lty=2)


#---------
# Bacon
plot(c(5,10,15,20),g.of.fit3['Bacon1',colnames(g.of.fit3) %in% 'elk1'][2:5],type='b', ylim =c(0,1), axes = F,ylab ='',xlab='')#, main ='Bacon Elk Lake')
sapply(2:5, function(z) {
points(c(5,10,15,20),g.of.fit3[paste('Bacon',z,sep=''),colnames(g.of.fit3) %in% 'elk1'][2:5],type='b',pch = z, col =z)
})
#legend('bottomleft',legend =legtest, col = 1:5, lwd = 1)
mtext (side =3, line = -1.5,adj = 1, 'b)',font = 2)
#axis(2)
#axis(1, at = c(5,10,15,20,30))
box()
#mtext(cex = 1.1,side =1, 'segment length',line = 2.2, font =2)
#mtext(cex = 1.1,side = 2, 'Proportion', line = 2.2, font = 2)
abline( h = 0.66,lty=2)
mtext(cex = 1.1,side = 3, 'Bacon', line = 0.5, font = 2)

#-----------
# Clam 
plot(rep(1.2,5),cl.el, col = 1:5,pch = 16,cex =2,xlim=c(1,2),ylim=c(0,1),axes = FALSE)
abline( h = 0.66,lty=2)
box()
mtext(cex = 1.1,side = 3, 'Clam and Bchron', line = 0.5, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'c)',font = 2)
#------
# Bchron
points(seq(1.5,1.9,length.out = 5),bc.el, col = 1:5,pch = 15,cex =1.5,xlim=c(1,2),ylim=c(0,1))
axis(4)
abline( h = 0.66,lty=2)
box()
mtext(cex = 1.1,side = 4, 'Proportion', line = 2.2, font = 2)
#mtext(cex = 1.1,side = 3, 'Bchron', line = 0.5, font = 2)
#mtext (side =3, line = -1.5,adj = 1, 'd)',font = 2)


#-------------------------------------------------------------------------------
# Holzmaar
#-------------
# OxCal
plot(c(0.1,0.25,0.5,1,2,3),g.of.fit3['OxCal1',colnames(g.of.fit3) %in% 'holzmaar1'],type='b', ylim =c(0,1), axes = F,ylab ='',xlab='')#, main ='OxCal Holzmaar')
sapply(2:5, function(z) {
points(c(0.1,0.25,0.5,1,2,3),g.of.fit3[paste('OxCal',z,sep=''),colnames(g.of.fit3) %in% 'holzmaar1'],type='b',pch = z, col =z)
})
#legend('bottomleft',legend =legtest, col = 1:5, lwd = 1)
axis(2)
#axis(1, at = c(0.1,0.25,0.5,1,2,3))
box()
#mtext(cex = 1.1,side =1, 'k',line = 2.2, font =2)
mtext(cex = 1.1,side = 2, 'Proportion', line = 2.2, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'd)',font = 2)
abline( h = 0.66,lty=2)


plot(c(5,10,15,20),g.of.fit3['Bacon1',colnames(g.of.fit3) %in% 'holzmaar1'][1:4],type='b', ylim =c(0,1), axes = F,ylab ='',xlab='')#, main ='Bacon Holzmaar')
sapply(2:5, function(z) {
points(c(5,10,15,20),g.of.fit3[paste('Bacon',z,sep=''),colnames(g.of.fit3) %in% 'holzmaar1'][1:4],type='b',pch = z, col =z)
})
#legend('bottomleft',legend =legtest, col = 1:5, lwd = 1)
mtext (side =3, line = -1.5,adj = 1, 'e)',font = 2)
#axis(2)
#axis(1, at = c(5,10,15,20,30))
box()
#mtext(cex = 1.1,side =1, 'segment length',line = 2.2, font =2)
#mtext(cex = 1.1,side = 2, 'Proportion', line = 2.2, font = 2)
abline( h = 0.66,lty=2)
#-----------------
#Clam
plot(rep(1.2, 5),cl.ho, col = 1:5,pch = 16,cex =2,ylim=c(0,1),xlim=c(1,2),axes = FALSE)
box()
abline(h = 0.66,lty=2)
mtext (side =3, line = -1.5,adj = 1, 'f)',font = 2)
#-----------------
# Bchron
points(seq(1.5,1.9,length.out = 5),bc.ho, col = 1:5,pch = 15,cex =1.5,ylim=c(0,1),xlim=c(1,2))
box()
axis(4)
abline(h = 0.66,lty=2)
#mtext (side =3, line = -1.5,adj = 1, 'h)',font = 2)
mtext(cex = 1.1,side = 4, 'Proportion', line = 2.2, font = 2)

#-------------------------------------------------------------------------------
# Iceberg Lake
#-------------
#OxCal
plot(c(0.1,0.25,0.5,1,2,3),g.of.fit3['OxCal1',colnames(g.of.fit3) %in% 'iceberg1'],type='b', ylim =c(0,1), axes = F,ylab ='',xlab='')#, main ='OxCal Iceberg')
sapply(2:5, function(z) {
points(c(0.1,0.25,0.5,1,2,3),g.of.fit3[paste('OxCal',z,sep=''),colnames(g.of.fit3) %in% 'iceberg1'],type='b',pch = z, col =z)
})
#legend('bottomleft',legend =legtest, col = 1:5, lwd = 1)
axis(2)
axis(1, at = c(0.1,0.5))
axis(1, at = c(0.25),labels = FALSE)
axis(1, at = c(1,2,3))
box()
mtext(cex = 1.1,side =1, 'k',line = 2.2, font =2)
mtext(cex = 1.1,side = 2, 'Proportion', line = 2.2, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'g)',font = 2)
abline( h = 0.66,lty=2)
#------------
#Bacon
plot(c(5,10,15,20),g.of.fit3['Bacon1',colnames(g.of.fit3) %in% 'iceberg1'][1:4],type='b', ylim =c(0,1), axes = F,ylab ='',xlab='')#, main ='Bacon Iceberg')
sapply(2:5, function(z) {
points(c(5,10,15,20),g.of.fit3[paste('Bacon',z,sep=''),colnames(g.of.fit3) %in% 'iceberg1'][1:4],type='b',pch = z, col =z)
})
legend('bottomleft',legend =legtest, col = 1:5, lwd = 2)
#axis(2)
axis(1, at = c(5,10,15,20,30))
box()
mtext(cex = 1.1,side =1, 'segment length  [cm]',line = 2.2, font =2)
mtext (side =3, line = -1.5,adj = 1, 'h)',font = 2)
#mtext(cex = 1.1,side = 2, 'Proportion', line = 2.2, font = 2)
abline( h = 0.66,lty=2)
#-------------------------------------------------------------------------------
# Bchron
plot(seq(1.2,1.8,length.out = 5),bc.ic, col = 1:5,pch = 15,cex =1.5,axes = FALSE, ylim=c(0,1),xlim=c(1,2))
box()
axis(4)
abline( h = 0.66,lty=2)
mtext(cex = 1.1,side = 4, 'Proportion', line = 2.2, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'i)',font = 2)
dev.off()