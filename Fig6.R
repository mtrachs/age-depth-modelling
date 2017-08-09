#-----------------------------------------------------------------------------------------------------------------
#Distributions of different age-depth models for selected depths of Elk Lake
#-----------------------------------------------------------------------------------------------------------------

jpeg(paste(fig.stor,'Fig6.jpg',sep=''),units='in',res=600, width = 13, height = 6.5)
par(mfrow=c(1,2))
sapply(c(1,9),function(x) {
#OxCal

oxcal.res <- lapply(c(4,13,111),function(t){
if(t!=111) {js <- models_err1[x+3,t][[1]]}
#open oxcal with prior therefore the random number 111 is needed
if ((t == 111)&(x==1)) {js <- as.matrix(oxcal.01[5,'elk'][[1]])}
if ((t == 111)&(x==9)) {js <- as.matrix(oxcal.01[2,'elk'][[1]])}


  z<-grep("data.z", js)#, extended=F)
  js1<-js[z]
  psq.depths<-as.vector(sapply(js1,function(b){
    bb<-unlist(strsplit(b,""))
    as.numeric(substring(b,which(bb=="=")+1, nchar(b)-1))
  }))

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
list(post.prob=post.prob,year = year,depths = psq.depths)
})
})


x.oxcal <- ifelse(x == 1,552,1332/2) #depth 500 cm in OxCal
oxcal.25 <- oxcal.res[[1]]
oxcal2 <- oxcal.res[[2]]
oxcal.prior <- oxcal.res[[3]]

oxcal.years.25 <- unlist(sapply(1:length(oxcal.25[,x.oxcal]$year),function(a) rep(oxcal.25[,x.oxcal]$year[a], ceiling(oxcal.25[,x.oxcal]$post.prob[a]*10^4)+1)))
oxcal.years2 <- unlist(sapply(1:length(oxcal2[,x.oxcal]$year),function(a) rep(oxcal2[,x.oxcal]$year[a], ceiling(oxcal2[,x.oxcal]$post.prob[a]*10^4)+1)))
oxcal.years.prior <- unlist(sapply(1:length(oxcal.prior[,x.oxcal]$year),function(a) rep(oxcal.prior[,x.oxcal]$year[a], ceiling(oxcal.prior[,x.oxcal]$post.prob[a]*10^4)+1)))


ft <- range(c(oxcal.years.25,oxcal.years2,oxcal.years.prior,models_err1[x,1][[1]][1000,],models_err1[x+1,1][[1]][,500],models_err1[x+2,1][[1]][,1000]))

density.clam <- density(models_err1[x,1][[1]][1000,],from = ft[1],to=ft[2])
density.OxCal.25 <- density(oxcal.years.25,from = ft[1],to=ft[2])
density.OxCal2 <- density(oxcal.years2,from = ft[1],to=ft[2])
density.OxCal.prior <- density(oxcal.years.prior,from = ft[1],to=ft[2])
density.Bacon5 <- density(models_err1[x+1,1][[1]][,500],from = ft[1],to=ft[2])
density.Bacon15 <- density(models_err1[x+1,10][[1]][,500],from = ft[1],to=ft[2])
density.Bchron <- density(models_err1[x+2,4][[1]][,1000],from = ft[1],to=ft[2])

y.max <- max(c(density.OxCal.25$y,density.OxCal2$y,density.OxCal.prior$y,density.clam$y,density.Bchron$y,density.Bacon5$y,density.Bacon15$y))


tt <- ifelse(x ==1,'5 radiocarbon dates','21 radiocarbon dates')

plot(density.clam,type='l',xlim=c(1650,3200),ylim=c(0,y.max),main = tt,xlab='Year cal BP',lwd=2)
lines(density.OxCal.25,col=2,lty=2,lwd=2)
lines(density.OxCal2,col='orange',lty=2,lwd=2)
lines(density.OxCal.prior,col='brown',lty=2,lwd=2)
lines(density.Bacon5,type='l',col='cyan',lty=3,lwd=2)
lines(density.Bacon15,type='l',col=4,lty=3,lwd=2)
lines(density.Bchron,type='l',col=3,lwd=2,lty=4)
if(x==9) {
legend('topleft',lty=c(1,2,2,2,3,3,4),lwd = 2,col=c(1:2,'orange','brown','cyan',4,3),c('clam','OxCal k = 0.25','OxCal k = 2','OxCal prior','Bacon seg.len = 5 cm','Bacon seg.len = 15 cm','Bchron'))
}                                                            
})
dev.off()


