#Load data from simulations
agm1<-with(datoutlst[[1]], aggregate(cbind(obs, est), list(sp=sp, sr=plantedsr), function(x) exp(mean(log(x[is.finite(log(x))]), na.rm=T))*(sum(is.finite(x) & x!=0))/length(x)))
agsd1<-with(datoutlst[[1]], aggregate(cbind(obs, est), list(sp=sp, sr=plantedsr), function(x) sd(log(x), na.rm=T)))
agm2<-with(datoutlst[[2]], aggregate(cbind(obs, est), list(sp=sp, sr=plantedsr), function(x) exp(mean(log(x[is.finite(log(x))]), na.rm=T))*(sum(is.finite(x) & x!=0))/length(x)))
agsd2<-with(datoutlst[[2]], aggregate(cbind(obs, est), list(sp=sp, sr=plantedsr), function(x) sd(log(x), na.rm=T)))

#plot mean +/- sd for each species
for(i in 1:length(splst)) {
  datout_ed<-datout
  
  #get mean and sd for each 
  agg_dat<-with(datout_ed, aggregate(cbind(obs, est), list(sp=sp, sr=plantedsr), hurdlemodel))
  agg_dat_se<-with(datout_ed, aggregate(cbind(obs, est), list(sp=sp, sr=plantedsr), function(x) sd(log(x), na.rm=T)))
  
  subsmu<-which(agg_dat$sp==splst[i])
  subssd<-which(agg_dat_se$sp==splst[i])
  
  sduse<-agg_dat_se[subssd,c("obs", "est")]
  sduse[is.na(sduse)]<-0
  
  #get min and max for plotting
  mx<-max(exp(log(agg_dat[subsmu,]$obs)+sduse$obs), na.rm=T)
  mn<-min(exp(log(agg_dat[subsmu,]$obs)-sduse$obs), na.rm=T)
  mn<-min(c(mn, min(exp(log(agm1[subsmu,]$est)-agsd1[subssd,]$est), na.rm=T)))
  mn<-min(c(mn, min(exp(log(agm2[subsmu,]$est)-agsd2[subssd,]$est), na.rm=T)))
  
  plot(c(0,4), c(mn, mx), axes=F, type="n", xlab="Planted Richness", ylab=expression(paste("Aboveground Biomass, g m"^2)), main=splst[i], log="y")
  if(i==2) {
    axis(2); axis(1, 0:4, c("1", "0", "0", "0", "16")); box()
  } else {
    axis(2); axis(1, 0:4, c(1, 2, 4, 8, 16)); box()
  }
  
  with(agg_dat[subsmu,], lines(log(sr,2), obs, lwd=2))
  polygon(c(log(agg_dat[subsmu,]$sr,2), rev(log(agg_dat[subsmu,]$sr,2))),
          c(exp(log(agg_dat[subsmu,]$obs)+sduse$obs),
            rev(exp(log(agg_dat[subsmu,]$obs)-sduse$obs))),
          col=adjustcolor("black", alpha.f = 0.5), border=NA)
  
  #plot with and without snapping
  for(usetr in 1:2) {
    datout<-datoutlst[[usetr]]
    agg_dat<-with(datout, aggregate(cbind(obs, est), list(sp=sp, sr=plantedsr), function(x) exp(mean(log(x[is.finite(log(x))]), na.rm=T))*(sum(is.finite(x) & x!=0))/length(x)))
    agg_dat_se<-with(datout, aggregate(cbind(obs, est), list(sp=sp, sr=plantedsr), function(x) sd(log(x), na.rm=T)))
    
    subsmum1<-subsmu
    
    upperpoly<-exp(log(agg_dat[subsmum1,]$est)+sduse[,]$est)
    lowerpoly<-exp(log(agg_dat[subsmum1,]$est)-sduse[,]$est)
    
    upperpoly<-upperpoly[is.finite(agg_dat[subsmum1,]$est)]
    lowerpoly<-lowerpoly[is.finite(agg_dat[subsmum1,]$est)]
    
    spprint<-log((agg_dat[subsmum1,]$sr)[is.finite(agg_dat[subsmum1,]$est)],2)
    
    with(agg_dat[subsmum1,], lines(log(sr,2), est, col=c("red", "blue")[usetr]))
    polygon(c(spprint, rev(spprint)),
            c(upperpoly, rev(lowerpoly)),
            col=adjustcolor(c("red", "blue")[usetr], alpha.f = 0.5), border=NA)
  }
}

par(mfrow=c(1,1))
mtext("Planted Richness", 1, line=3.5)
mtext(expression(paste("Aboveground Biomass, g m"^2)), 2, line=3.5)

