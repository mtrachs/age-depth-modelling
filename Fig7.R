#----------------------------------------------------------------------------------------------------------------------
#accumulation rates
#----------------------------------------------------------------------------------------------------------------------

colnames(models_err1) <- paste('',rep(c('elk','holzmaar','iceberg'),6),'',rep(1:6,each = 3),sep='')
colnames(rad.carb) <- colnames(models_err1)
rownames(models_err1) <- rep(c('clam','Bacon','Bchron','OxCal'),5)

chrons<- sapply(colnames(models_err1)[1:3], function(x) {
  sapply(1:20, function(y) {
      xx <- substring(x,1,nchar(x)-1)
      if((rownames(models_err1)[y] == 'Bacon') | (rownames(models_err1)[y] == 'clam') | (rownames(models_err1)[y] == 'Bchron')) {
      #define variable that is needed to access radiocarbon ages
      fa <- ceiling(y/4)
      #generate depths at which predictions of the age-depth model were made
      pred.depth <- seq(min(lakes[[xx]][,'Depth']),max(rad.carb[fa,x][[1]][,'depth']),4) 
      # find depths that do exist
      real.depth <- pred.depth %in% unique(lakes[[xx]][,'Depth'])
      # make a data series where one depth only occurs once (ages with same depths are averaged)
      ag.lakes <- aggregate(lakes[[xx]][,'Age'], by = as.data.frame(lakes[[xx]][,'Depth']), FUN = mean)
      colnames(ag.lakes) <- c('Depth','Age')
      # subset the refined series so that only existing depths are used
      ag.lakes <- ag.lakes[ag.lakes[,'Depth']%in% pred.depth[real.depth],]

      #extract age-depth models at 4 cm increment
      if(xx == 'elk') {
        if(rownames(models_err1)[y] == 'clam')  chron <- round(models_err1[y,x][[1]][seq(1,nrow(models_err1[y,x][[1]]),8),])
        if(rownames(models_err1)[y] == 'Bacon') chron <- t(round(models_err1[y,x][[1]]))[seq(1,ncol(models_err1[y,x][[1]]),4),round(seq(1,nrow(models_err1[y,x][[1]]),length.out = 1000))]     
        if(rownames(models_err1)[y] == 'Bchron') chron <- t(matrix(unlist(models_err1[y,x]),nrow=1000))[seq(1,(length(unlist(models_err1[y,x][[1]]))/1000),8),] 
      }        
      
      if(xx %in% c('holzmaar','iceberg')) {
        if(rownames(models_err1)[y] == 'clam')  chron <- round(models_err1[y,x][[1]][seq(1,nrow(models_err1[y,x][[1]]),8),])
        if(rownames(models_err1)[y] == 'Bacon') chron <- t(round(models_err1[y,x][[1]]))[seq(1,ncol(models_err1[y,x][[1]]),2),round(seq(1,nrow(models_err1[y,x][[1]]),length.out = 1000))]
        if(rownames(models_err1)[y] == 'Bchron') chron <- t(matrix(unlist(models_err1[y,x]),nrow=1000))[seq(1,(length(unlist(models_err1[y,x][[1]]))/1000),8),]
      }
chron 
  }
  })
})


colnames(chrons) <- colnames(models_err1)[1:3]
rownames(chrons) <- rownames(models_err1)

chrons <- chrons[(rownames(chrons) %in% 'OxCal')==FALSE,]

rownames(chrons) <- paste(rownames(chrons),rep(c(6,11,21,31,41),each = 3), sep ='')


acc.rates <- lapply(colnames(chrons), function(y) {
  lapply(rownames(chrons), function(x) {
          acc.rate <- apply(chrons[x,y][[1]],2,diff) #calculate accumulation rates in year/four cm
          bounds <- apply(acc.rate,1,function(x) quantile(x, probs = c(0.025,0.5,0.975))) 
  })
})



#-------------------------------------------------------------------------------
#find the reall accumulation rates for each lake (and going to the depth covered with radiocarbon dates, only needed once for each combination)
real.acc <- sapply(colnames(models_err1)[1:3],function(x) {
  sapply(seq(1,20,4), function(y) {
    xx <-  substring(x,1,nchar(x)-1)
    fa <- ceiling(y/4)
    pred.depth <- seq(min(lakes[[xx]][,'Depth']),max(rad.carb[fa,x][[1]][,'depth']),4) # damit man nicht die Extrapolationen vergleicht
    real.depth <- pred.depth %in% unique(lakes[[xx]][,'Depth'])
    ag.lakes <- aggregate(lakes[[xx]][,'Age'], by = as.data.frame(lakes[[xx]][,'Depth']), FUN = mean)
    colnames(ag.lakes) <- c('Depth','Age')
    ag.lakes <- ag.lakes[ag.lakes[,'Depth']%in% pred.depth[real.depth],]
    diff(ag.lakes$Age)
  })
})

#-------------------------------------------------------------------------------
jpeg(paste(fig.stor,'Fig6_accumulation_rate.jpg',sep=''), width = 18, height = 9,units='in',res = 600)
par(mfrow=c(1,9))
par(oma=c(4,4,0,1), mar = c(0,0,2,0),cex = 1.25)
titel <- paste(rep(c('clam','Bacon','Bchron'),5),rep(c(6,11,21,31,41),each = 3))

sapply(2,function(y) { #only plot holzmaar results
    sapply(c(1:3,7:9,13:15), function(x) { # plot accumulation rates for 6,21 and 41 radiocarbon dates
     laenge <- length(real.acc[ceiling(x/3),ifelse(y%%3 == 0,3,y%%3)][[1]])
     xlims <- c(0,300)
     ylims <- c(900,121)
     plot(acc.rates[[y]][[x]][2,1:laenge],121 + 4*(1:laenge),type='l',xlim=xlims,ylim=ylims, axes = FALSE, main = titel[x])
     axis(1,at=c(100,200,300))
      axis(1,at=0,labels = FALSE)
     if(x ==1) { axis(2, at = seq(100,900,100))
     mtext(cex = 1.25,side = 2, 'Depth [cm]',line = 2.2, font = 2)
     }
    
   
     box()
     mtext(cex = 1.25,side = 1,'#Year / 4cm',line = 2.2 ,font = 2)
     polygon(c(acc.rates[[y]][[x]][1,1:laenge],rev(acc.rates[[y]][[x]][3,1:laenge])),121 +4*c(1:laenge,rev(1:laenge)),col='lightblue')
     lines(acc.rates[[y]][[x]][2,1:laenge], 121 +4*1:laenge,col='blue', lwd = 2)
     lines(real.acc[ceiling(x/3),ifelse(y%%3 == 0,3,y%%3)][[1]],121 +4*(1:laenge),lwd = 2)
     abline(h = rad.carb[ceiling(x/3),y][[1]][,'depth'],lty=2)
     if(x %in% 1:3) mtext(side = 3, adj = 0,line=-1.15,paste(letters[x],')',sep=''),font=2,cex = 1.3)
     if(x %in% 7:9) mtext(side = 3, adj = 0,line=-1.15,paste(letters[x-3],')',sep=''),font=2,cex = 1.3)
     if(x %in% 13:15) mtext(side = 3, adj = 0,line=-1.15,paste(letters[x-6],')',sep=''),font=2,cex = 1.3)
  })
})
dev.off()




