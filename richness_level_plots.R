##########################################
#Richness-level plots
##########################################
par(mfrow=c(2,2), oma=c(1,2,2,0))
for(isp in c(2, 4, 8, 16)) {
  subs<-which(datout$plantedsr==isp) #plot only for non-monocultures
  subs2<-which(is.finite(log(datout[subs,]$est)) & is.finite(log(datout[subs,]$obs)))
  
  plot(log10(obs)~log10(est), datout[subs,][subs2,],
       main=paste("Planted Richness =", isp), type="n",
       axes=F, xlab="", ylab="",
       cex.axis=1.5)
  put.fig.letter(label = c("A.", "B.", "C.", "D.")[which(c(2,4,8,16)==isp)], location = "topleft", cex=2, x=0.04, y=0.95)
  
  sq<-c(0.1, 1, 10, 100)
  axis(1, log10(sq), sq, cex.axis=1.5)
  axis(2, log10(sq), sq, cex.axis=1.5)
  box()
  abline(a=0, b=1, col="darkgrey", lwd=2, lty=3)
  
  #fit ranged major axis regression
  obs<-log10(datout[subs,]$obs)
  est<-log10(datout[subs,]$est)
  obs[!is.finite(obs)]<-quantile(obs[is.finite(obs)], 0.025)
  est[!is.finite(est)]<-quantile(est[is.finite(est)], 0.025)
  
  lmd<-lmodel2(obs~est, range.y="interval", range.x="interval", nperm=nrep)
  abline(a=lmd$regression.results[4,2], b=lmd$regression.results[4,3], lwd=2)
  abline(a=mean(obs)-mean(est)*lmd$confidence.intervals[4,5],
         b=lmd$confidence.intervals[4,5], lwd=2, lty=2)
  abline(a=mean(obs)-mean(est)*lmd$confidence.intervals[4,4],
         b=lmd$confidence.intervals[4,4], lwd=2, lty=2)
  
  points(log10(datout[subs,]$est[subs2]), log10(datout[subs,]$obs[subs2]),
         pch=16,
         col=grey.colors(4, start=0, end=1)[as.numeric(as.factor(datout[subs,]$plantedsr))],
         lwd=2, cex=0.8)
  points(log10(datout[subs,]$est[subs2]), log10(datout[subs,]$obs[subs2]),
         pch=1,
         col=1,
         lwd=1, cex=0.8)
  
  #plot CD
  rsqest<-lmd$rsquare
  
  r<-"R"
  p<-"p"
  rd<-paste(" =", round(rsqest, 2))
  pv<-paste(" <", ceiling(lmd$regression.results[4,5]*1000)/1000)
  
  z<-log10(datout[subs,]$obs[subs2])
  zfit<-log10(datout[subs,]$est[subs2])
  text((max(zfit)-min(zfit))*0.1+min(zfit), (max(z)-min(z))*0.95+min(z),
       bquote(.(r[1])^2 ~ .(rd[1])), cex=1.5, pos=4)
  text((max(zfit)-min(zfit))*0.1+min(zfit), (max(z)-min(z))*0.85+min(z),
       bquote(.(p[1]) ~ .(pv[1])), cex=1.5, pos=4)
}
par(mfrow=c(1,1))
mtext("Species-Level Biomass", 3, line=4, cex=1.5)
mtext(expression(paste("Predicted Biomass, g m"^"-2", sep="")), 1, line=4, cex=1.5)
mtext(expression(paste("Observed Biomass, g m"^"-2", sep="")), 2, line=4, cex=1.5)







par(mfrow=c(2,2), oma=c(1,2,2,0))
for(isp in c(2, 4, 8, 16)) {
  agdat_plot<-aggregate(cbind(est=datout$est, obs=datout$obs),
                        list(plantedsr=datout$plantedsr, plt=datout$plt),
                        function(x) sum(x, na.rm=T))
  agdat_plot<-agdat_plot[agdat_plot$plantedsr==isp,] #plot only for non-monocultures
  
  plot(log10(obs)~log10(est), agdat_plot,
       main=paste("Planted Richness =", isp), type="n",
       xlab="", ylab="",
       axes=F,
       cex.axis=1.5)
  put.fig.letter(label = c("A.", "B.", "C.", "D.")[which(c(2,4,8,16)==isp)], location = "topleft", cex=2, x=0.04, y=0.95)
  
  sq<-c(25, 50, 100, 200, 400)
  axis(1, log10(sq), sq, cex.axis=1.5)
  axis(2, log10(sq), sq, cex.axis=1.5)
  box()
  
  abline(a=0, b=1, col="darkgrey", lwd=2, lty=3)
  
  #fit ranged major axis regression
  obs<-log10(agdat_plot$obs)
  est<-log10(agdat_plot$est)
  obs[!is.finite(obs)]<-quantile(obs[is.finite(obs)], 0.025)
  est[!is.finite(est)]<-quantile(est[is.finite(est)], 0.025)
  
  lmd<-lmodel2(obs~est, range.y="interval", range.x="interval", nperm=nrep)
  abline(a=lmd$regression.results[4,2], b=lmd$regression.results[4,3], lwd=2)
  abline(a=mean(obs)-mean(est)*lmd$confidence.intervals[4,5],
         b=lmd$confidence.intervals[4,5], lwd=2, lty=2)
  abline(a=mean(obs)-mean(est)*lmd$confidence.intervals[4,4],
         b=lmd$confidence.intervals[4,4], lwd=2, lty=2)
  
  with(agdat_plot, points(log10(est), log10(obs),
                          pch=16,
                          col=grey.colors(4, start=0, end=1)[as.numeric(as.factor(plantedsr))],
                          lwd=2, cex=0.8))
  with(agdat_plot, points(log10(est), log10(obs),
                          pch=1,
                          col=1,
                          lwd=1, cex=0.8))
  
  #plot CD
  rsqest<-lmd$rsquare
  
  r<-"R"
  p<-"p"
  rd<-paste(" =", round(rsqest, 2))
  pv<-paste(" <", ceiling(lmd$regression.results[4,5]*1000)/1000)
  
  z<-log10(agdat_plot$obs)
  zfit<-log10(agdat_plot$est)
  text((max(zfit)-min(zfit))*0.1+min(zfit), (max(z)-min(z))*0.95+min(z),
       bquote(.(r[1])^2 ~ .(rd[1])), cex=1.5, pos=4)
  text((max(zfit)-min(zfit))*0.1+min(zfit), (max(z)-min(z))*0.85+min(z),
       bquote(.(p[1]) ~ .(pv[1])), cex=1.5, pos=4)
  
}

par(mfrow=c(1,1))
mtext("Total Plot Biomass", 3, line=4, cex=1.5)
mtext(expression(paste("Predicted Biomass, g m"^"-2", sep="")), 1, line=4, cex=1.5)
mtext(expression(paste("Observed Biomass, g m"^"-2", sep="")), 2, line=4, cex=1.5)