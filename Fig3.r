#---------------------------------------------------------------------------------------------------------------------------
#Calculate RMSEP and produce Figure 3
#------------------------------------------------------------------------------------------------------------------------

rownames(models_err1) <- rep(c('clam','Bacon','Bchron','OxCal'),5)

#calculate the root mean square error of prediction for all models
fehler <- sapply(colnames(models_err1), function(x) {
  sapply(1:20, function(y) {
      xx <- substring(x,1,nchar(x)-1)
      if((rownames(models_err1)[y] == 'Bacon') | (rownames(models_err1)[y] == 'clam')|(rownames(models_err1)[y] == 'Bchron')) {
        #define variable that is needed to access radiocarbon ages
        fa <- ceiling(y/4)
        #generate depths at which predictions of the age-depth model were made
        pred.depth <- seq(min(lakes[[xx]][,'Depth']),max(rad.carb[fa,x][[1]][,'depth']),2) # damit man nicht die Extrapolationen vergleicht
        # find depths that do exist
        real.depth <- pred.depth %in% unique(lakes[[xx]][,'Depth'])
        # make a data series where one depth only occurs once (ages with same depths are averaged)
        ag.lakes <- aggregate(lakes[[xx]][,'Age'], by = as.data.frame(lakes[[xx]][,'Depth']), FUN = mean)
      colnames(ag.lakes) <- c('Depth','Age')
      # subset the refined series so that only existing depths are used
      ag.lakes <- ag.lakes[ag.lakes[,'Depth']%in% pred.depth[real.depth],]

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

#calculate the mean
mittel <- rowMeans(chron[1:length(pred.depth),])[real.depth]
rmsep <-   sqrt(mean((ag.lakes[,'Age']-mittel)^2))
mad <-   mean(abs(ag.lakes[,'Age']-mittel))
max.f <-  max(abs(ag.lakes[,'Age']-mittel))
}

#-------------------------------------------------------------------------------
#Oxcal
if(rownames(models_err1)[y] == 'OxCal') {
  js <- models_err1[y,x][[1]]
  # find depths 
  z<-grep("data.z", js)#
  js1<-js[z]
  psq.depths<-as.vector(sapply(js1,function(b){
    bb<-unlist(strsplit(b,""))
    as.numeric(substring(b,which(bb=="=")+1, nchar(b)-1))
  }))
 
  #starting years
  sta<-grep("start", js)#, extended=F)
  js1<-js[sta]
  psq.sta<-as.vector(sapply(js1,function(b){
  bb<-unlist(strsplit(b,""))
  as.numeric(substring(b,which(bb=="=")+1, nchar(b)-1))
  }))
  #find depths at which no real data is available
  ll <- which(is.na(psq.sta)== TRUE)
 
 
 #find mean of the age depth model
 z<-grep("posterior.mean", js)
  js1<-js[z]
  mittel<-as.vector(sapply(js1,function(b){
    bb<-unlist(strsplit(b,""))
    as.numeric(substring(b,which(bb=="=")+1, nchar(b)-1))
  }))

  
  mittel <- mittel[-(1:ll)]
  mittel <- abs(mittel-1950)
  psq.depths <- psq.depths[-(1:ll)]
  
  
ag.lakes <- aggregate(lakes[[xx]][,'Age'], by = as.data.frame(lakes[[xx]][,'Depth']), FUN = mean)
colnames(ag.lakes) <- c('Depth','Age')
#approx is used to make sure the depth scale is inverted properly (OxCal has lowermost values first)
int.lakes <- approx(ag.lakes[,'Depth'],ag.lakes[,'Age'],xout = psq.depths[psq.depths%in%ag.lakes$Depth])
age.of.int <-   int.lakes$y

rmsep <- sqrt(mean((age.of.int-mittel[psq.depths%in%ag.lakes$Depth])^2))
mad <-   mean(abs(age.of.int-mittel[psq.depths%in%ag.lakes$Depth]))
max.f <-  max(abs(age.of.int-mittel[psq.depths%in%ag.lakes$Depth]))


if(xx=='iceberg') {
  rmsep <- sqrt(mean(((age.of.int-mittel[psq.depths%in%ag.lakes$Depth])[seq(1,length(age.of.int),2)])^2))
  mad <-   mean(abs((age.of.int-mittel[psq.depths%in%ag.lakes$Depth])[seq(1,length(age.of.int),2)]))
  max.f <-  max(abs((age.of.int-mittel[psq.depths%in%ag.lakes$Depth])[seq(1,length(age.of.int),2)]))
}
}
list(rmsep = round(rmsep), mad = round(mad),max.f = round(max.f))
})
})

rownames(fehler) <-      paste(rep(c('rmsep','mad','max'),20),rep(rep(c('clam','Bacon','Bchron','OxCal'), each = 3),5),rep(c(60,11,21,31,41),each = 12), sep ='.')
colnames(fehler) <- rep(colnames(models_err1)[1:3],6)

#-------------------------------------------------------------------------------
aa <- substring(rownames(fehler),1,nchar(rownames(fehler))-3)
#find all values for clam models
fehler.clam.el <- matrix(ncol=5,unlist(fehler[aa%in%'rmsep.clam',colnames(fehler)%in%'elk1'][,2:6]))
fehler.clam.ho <- matrix(ncol=5,unlist(fehler[aa%in%'rmsep.clam',colnames(fehler)%in%'holzmaar1'][,2:6]))
fehler.clam.ic <- matrix(ncol=5,unlist(fehler[aa%in%'rmsep.clam',colnames(fehler)%in%'iceberg1'][,2:6]))

#find all values for bchron models
fehler.bchron.el <- matrix(ncol=5,unlist(fehler[aa%in%'rmsep.Bchron',colnames(fehler)%in%'elk1'][,2:6]))
fehler.bchron.ho <- matrix(ncol=5,unlist(fehler[aa%in%'rmsep.Bchron',colnames(fehler)%in%'holzmaar1'][,2:6]))
fehler.bchron.ic <- matrix(ncol=5,unlist(fehler[aa%in%'rmsep.Bchron',colnames(fehler)%in%'iceberg1'][,2:6]))

#take median of clam and bchron models
cl.el <- apply(fehler.clam.el,1,median)
cl.ho <- apply(fehler.clam.ho,1,median)
cl.ic <- apply(fehler.clam.ic,1,median)
bc.el <- apply(fehler.bchron.el,1,median)
bc.ho <- apply(fehler.bchron.ho,1,median)
bc.ic <- apply(fehler.bchron.ic,1,median)

ncal <- c(6,11,21,31,41)
leg <- sapply(1:5, function(z) {
              paste('n = ',ncal[z],',',sep='')
       })

legtest <- as.matrix(unlist(strsplit(leg,',')))
#change rownames from 6 to 60 
rownames(fehler) <- paste(rep(c('rmsep','mad','max'),20),rep(rep(c('clam','Bacon','Bchron','OxCal'), each = 3),5),rep(c(6,11,21,31,41),each = 12), sep ='.')


jpeg(paste(fig.stor,'Fig3_rmsep.jpg',sep=''), width = 12, height = 12,units='in',res=600)
layout(rbind(1:3,4:6,7:9),widths = c(3,3,1),heights=rep(3,3))
par(oma = c(3,3,2,3),mar=c(0.5,0.5,0.5,0.5),cex = 1.15) # mfrow=c(3,2)


plot(c(0.1,0.25,0.5,1,2,3),fehler['rmsep.OxCal.6',colnames(fehler) %in% 'elk1'],type='b', ylim =c(0,170), axes = F,ylab ='',xlab='')#, main ='OxCal Elk Lake')
sapply(2:5, function(z) {
points(c(0.1,0.25,0.5,1,2,3),fehler[paste('rmsep.OxCal.',ncal[z],sep=''),colnames(fehler) %in% 'elk1'],type='b',pch = z, col =z)
})
legend('topleft',legend =legtest, col = 1:5, lwd = 2)
axis(2)
box()
mtext(cex =1.1,side = 2, 'RMSEP [Year]', line = 2.2, font = 2)
mtext(cex =1.1,side = 3, 'OxCal', line = 0.5, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'a)',font = 2)

plot(c(5,10,15,20),fehler['rmsep.Bacon.6',colnames(fehler) %in% 'elk1'][2:5],type='b', ylim =c(0,170), axes = F,ylab ='',xlab='')#, main ='Bacon Elk Lake')
sapply(2:5, function(z) {
points(c(5,10,15,20),fehler[paste('rmsep.Bacon.',ncal[z],sep=''),colnames(fehler) %in% 'elk1'][2:5],type='b',pch = z, col =z)
})
points(seq(22,28,length.out = 5),bc.el, col = 1:5,pch = 15,cex =2)
box()
mtext(cex =1.1,side = 3, 'Bacon', line = 0.5, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'b)',font = 2)

#--------
#clam
plot(rep(1.3, 5),cl.el, col = 1:5,pch = 16,cex =2,ylim=c(0,170),xlim=c(1,2),axes =FALSE)
box()
mtext(cex =1.1,side = 3, 'Clam and Bchron', line = 0.5, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'c)',font = 2)
#-------------
#Bchron
points(rep(1.7, 5),bc.el, col = 1:5,pch = 15,cex =1.5,ylim=c(0,170),xlim=c(1,2))
mtext(cex =1.1,side = 4, 'RMSEP [Year]', line = 2.2, font = 2)
box()
axis(4)


#-------------------------------------------------------------------------------
# Holzmaar
plot(c(0.1,0.25,0.5,1,2,3),fehler['rmsep.OxCal.6',colnames(fehler) %in% 'holzmaar1'],type='b', ylim =c(0,300), axes = F,ylab ='',xlab='')#, main ='OxCal Holzmaar')
sapply(2:5, function(z) {
points(c(0.1,0.25,0.5,1,2,3),fehler[paste('rmsep.OxCal.',ncal[z],sep=''),colnames(fehler) %in% 'holzmaar1'],type='b',pch = z, col =z)
})
axis(2)
box()
mtext(cex =1.1,side = 2, 'RMSEP [Year]', line = 2.2, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'd)',font = 2)


plot(c(5,10,15,20),fehler['rmsep.Bacon.6',colnames(fehler) %in% 'holzmaar1'][1:4],type='b', ylim =c(0,300), axes = F,ylab ='',xlab='')#, main ='Bacon Holzmaar')
sapply(2:5, function(z) {
points(c(5,10,15,20),fehler[paste('rmsep.Bacon.',ncal[z],sep=''),colnames(fehler) %in% 'holzmaar1'][1:4],type='b',pch = z, col =z)
})
box()
mtext (side =3, line = -1.5,adj = 1, 'e)',font = 2)
#------
#--------
#clam
plot(rep(1.3, 5),cl.ho, col = 1:5,pch = 16,cex =2,ylim=c(0,300),xlim=c(1,2),axes =FALSE)
box()
mtext (side =3, line = -1.5,adj = 1, 'f)',font = 2)
#-------------
#Bchron
points(rep(1.7,5),bc.ho, col = 1:5,pch = 15,cex =1.5,ylim=c(0,300),xlim=c(1,2))
mtext(cex =1.1,side = 4, 'RMSEP [Year]', line = 2.2, font = 2)
box()
axis(4)



#-------------------------------------------------------------------------------
#iceberg
plot(c(0.1,0.25,0.5,1,2,3),fehler['rmsep.OxCal.6',colnames(fehler) %in% 'iceberg1'],type='b', ylim =c(0,55), axes = F,ylab ='',xlab='')#, main ='OxCal Iceberg')
sapply(2:5, function(z) {
points(c(0.1,0.25,0.5,1,2,3),fehler[paste('rmsep.OxCal.',ncal[z],sep=''),colnames(fehler) %in% 'iceberg1'],type='b',pch = z, col =z)
})
axis(2)
axis(1, at = c(0.1,0.5))
axis(1, at = c(0.25),labels = FALSE)
axis(1, at = c(1,2,3))
box()
mtext(cex =1.1,side =1, 'k',line = 2.2, font =2)
mtext(cex =1.1,side = 2, 'RMSEP [Year]', line = 2.2, font = 2)
mtext (side =3, line = -1.5,adj = 1, 'g)',font = 2)


plot(c(5,10,15,20),fehler['rmsep.Bacon.6',colnames(fehler) %in% 'iceberg1'][1:4],type='b', ylim =c(0,55), axes = F,ylab ='',xlab='')#, main ='Bacon Iceberg')
sapply(2:5, function(z) {
points(c(5,10,15,20),fehler[paste('rmsep.Bacon.',ncal[z],sep=''),colnames(fehler) %in% 'iceberg1'][1:4],type='b',pch = z, col =z)
})
axis(1, at = c(5,10,15,20,30,40))
box()
mtext(cex =1.1,side =1, 'segment length [cm]',line = 2.2, font =2)
mtext (side =3, line = -1.5,adj = 1, 'h)',font = 2)
#-------------
#Bchron
plot(rep(1.5, 5),bc.ic, col = 1:5,pch = 15,cex =1.5,ylim=c(0,55),xlim=c(1,2),axes = FALSE)
mtext(cex =1.1,side = 4, 'RMSEP [Year]', line = 2.2, font = 2)
box()
axis(4)
mtext (side =3, line = -1.5,adj = 1, 'i)',font = 2)
dev.off()