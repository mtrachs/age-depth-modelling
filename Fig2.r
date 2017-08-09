#--------------------------------------------------------------------------------------------------------------------
#Plot some age-depth models
#--------------------------------------------------------------------------------------------------------------------
# To start, either load the data provided or run the file generate dates
#-----------------------------------------------------------------------------------------------------------------

load('~/age-depth/codeupload/data.RData')
fig.stor <- '~/age-depth/'

#matrix needed to reshuffle the order of the output of age-depth models in object models err1
# we want the order clam, Bchron, Bacon,bacon, OxCal,OxCal,OxCal, this is achieved using this ratheer strang ematrix
dd <- cbind(rownames(models_err1)[c(1,3,2,4)],rownames(models_err1)[c(9,11,10,12)])
rownames(dd) <- c('clam','Bchron','Bacon','OxCal')
#access object from 
ee <- c(5,2)

plot.name <- c('holzmaar_6', 'holzmaar_21')
#par(mfrow=c(1,3))
jpeg(paste(fig.stor, 'Fig2.jpg',sep=''), width = 13, height = 16,units='in', res = 300)
par(mfrow=c(2,1),cex=2,mar=c(1,2,1,1),oma=c(1,1,0,0))
sapply(1:2,function(bb) {

#construct the data we actually want to use, the order of this data is different from the order in which we want to plot 
#dd and ee are matrices constructed above,   
aa <- rbind(as.matrix(models_err1[dd[1:3,bb],'holzmaar1' ]),as.matrix(models_err1[dd['Bacon',bb],"holzmaar4"]),as.matrix(models_err1[dd['OxCal',bb],'holzmaar2']),as.matrix(models_err1[dd['OxCal',bb],'holzmaar4']), as.matrix(oxcal.01[ee[bb],2]))
models_err_xx <- cbind(aa,aa)
rownames(models_err_xx) <- c('clam','Bchron',paste(rep('Bacon',2),c(10,20),sep='_'),paste(rep('OxCal',3),c(0.25,2,'prior'),sep='_'))

#extract the radiocarbon dates and find the depth of the lowest 14C date
radiocarbon <- rad.carb[,'holzmaar1'] 
if (bb == 1)   depth.rc <- (as.data.frame(radiocarbon[[1]])$depth)
if (bb == 2)   depth.rc <- (as.data.frame(radiocarbon[[2]])$depth)
#if (bb == 3)   depth.rc <- (as.data.frame(radiocarbon[[4]])$depth)
max.depth <- max (depth.rc)
colnames(models_err_xx) <- rep('holzmaar',2)


sapply('holzmaar', function(x) {
#x11(13,8)
        sapply(1:7, function(y1) {

y <- rownames(models_err_xx)[y1] 

#some data handling
if(y %in% 'clam') {chron <- round(models_err_xx[y,x][[1]])}
if(strsplit(y,'_')[[1]][1] %in% 'Bacon') {chron <- t(round(models_err_xx[y,x][[1]]))}
if(y %in% 'Bchron') {chron <- t(matrix(unlist(models_err_xx[y,x]),nrow=1000))}
max.plot.depth <- which(lakes[[x]][,'Depth'] == max.depth)

pred.depth <- seq(min(lakes[[x]][,'Depth']),max.depth,2)

if((y %in% 'clam') | (y %in% 'Bchron')) chron <- chron [seq(1,4*length(pred.depth),4),]
if(strsplit(y,'_')[[1]][1] %in% 'Bacon') chron <- chron [1:length(pred.depth),]

if((y %in% 'clam') | (y %in% 'Bchron')| (strsplit(y,'_')[[1]][1] %in% 'Bacon')) {
bounds <- apply(chron,1, function(z) quantile(z,probs = c(0.025,0.975,0.166,0.833),na.rm = TRUE))
mittel <- rowMeans(chron)

par(mgp=c(3,1,0))
if(y1 >1) {
par(new = TRUE)
}
plot(lakes[[x]][,'Age'],lakes[[x]][,'Depth'],type='n', ylim =rev(range(pred.depth)), xlim =c(11000, min(lakes[[x]][,'Age'])-14000), xlab ='', ylab ='',col='red',axes = FALSE)
polygon(c(bounds[1,]-2000*(y1),rev(bounds[2,]-2000*(y1))),c(pred.depth,rev(pred.depth)),col='grey',border = NA)
polygon(c(bounds[3,]-2000*(y1),rev(bounds[4,]-2000*(y1))),c(pred.depth,rev(pred.depth)),col='black',border = NA)
lines((lakes[[x]][1:max.plot.depth,'Age']-2000*(y1)),lakes[[x]][1:max.plot.depth,'Depth'],col ='blue',lwd = 2)
lines(mittel-2000*(y1),pred.depth,col ='red',lwd = 2)
mtext(side = 1, font = 2, 'Age [cal BP]',line =0.1,cex=2)
mtext(side = 2, font = 2, 'Depth [cm]',line =2.2,cex=2)
if (y1 == 1){ abline(h = depth.rc, lty = 2)
axis(2, at=seq(100,900,100))
#mtext('a)',side=3, adj=0)
par(mgp=c(3,-13,-14))
axis(1, at=seq(10000,0,-2000))
box()
}
if(bb==1)mtext('a)',side=3, adj=0,cex=2)
if(bb==2)mtext('b)',side=3, adj=0,cex=2)
}




#code for OxCal results
if(strsplit(y,'_')[[1]][1] %in% 'OxCal') {
  js <- models_err_xx[y,x][[1]]
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
psq.sta1 <- psq.sta[(ll+1):length(psq.sta)]
psq.depths<-unique(psq.depths[-(1:(ll))])


    pro<-     grep("posterior.prob=", js)#, extended=F)
  js1<-js[pro]
  psq.pro<-as.vector(sapply(js1,function(b){
    bb<-unlist(strsplit(b,""))
    cc <- substring(b,which(bb=="=")+2,nchar(b)-2)
    as.numeric(unlist(strsplit(cc,",")))
  }))

 z<-grep("posterior.mean", js)#, extended=F)
  js1<-js[z]
  mittel<-as.vector(sapply(js1,function(b){
    bb<-unlist(strsplit(b,""))
    as.numeric(substring(b,which(bb=="=")+1, nchar(b)-1))
  }))

  mittel <- mittel[-(1:ll)]
  mittel <- abs(mittel-1950)


oxcal.elk.31 <- sapply(1:length(psq.sta1), function(x){
        post.prob <- psq.pro[[x]] / sum(psq.pro[[x]])
        year <- seq(from=round(abs(psq.sta1[x]-1950)), by = -5, length.out =length(psq.pro[[x]]))
list(post.prob=post.prob,year = year)
})

bounds <- sapply(seq(2,length(oxcal.elk.31),2),function(x) {
    cs <- cumsum(oxcal.elk.31[[x-1]])
    lower <- oxcal.elk.31[[x]][min(which(cs > 0.025))]
    upper <- oxcal.elk.31[[x]][max(which(cs < 0.975))]
    lower1 <- oxcal.elk.31[[x]][min(which(cs > 0.166))]
    upper1 <- oxcal.elk.31[[x]][max(which(cs < 0.833))]
    c(lower,upper,lower1,upper1)
})
par(new = TRUE)
plot(lakes[[x]][,'Age'],lakes[[x]][,'Depth'],type='n', ylim =rev(range(pred.depth)), xlim =c(11000, min(lakes[[x]][,'Age'])-14000),xlab ='', ylab ='',col='red', axes = FALSE)
polygon(c(bounds[1,-c(1:(ll))]-2000*(y1),rev(bounds[2,-c(1:(ll))]-2000*(y1))),c(psq.depths,rev(psq.depths)), col = 'grey' ,border = NA)
polygon(c(bounds[3,-c(1:(ll))]-2000*(y1),rev(bounds[4,-c(1:(ll))]-2000*(y1))),c(psq.depths,rev(psq.depths)), col = 'black' ,border = NA)
lines(lakes[[x]][1:max.plot.depth,'Age']-2000*(y1),lakes[[x]][1:max.plot.depth,'Depth'],col ='blue',lwd=2)
lines(mittel-2000*(y1),psq.depths,col ='red',lwd=2)
}
})
})
})
dev.off()

