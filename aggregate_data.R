#Set up plot
m<-rbind(c(1,1,1,1,2,2,2,2),
         c(1,1,1,1,2,2,2,2),
         c(1,1,1,1,2,2,2,2),
         c(1,1,1,1,2,2,2,2),
         c(3,3,3,3,4,4,4,4),
         c(3,3,3,3,4,4,4,4),
         c(3,3,3,3,4,4,4,4),
         c(3,3,3,3,4,4,4,4),
         c(5,5,5,5,6,6,6,6),
         c(7,7,7,7,8,8,8,8),
         c(7,7,7,7,8,8,8,8),
         c(7,7,7,7,8,8,8,8),
         c(7,7,7,7,8,8,8,8))
layout(m)

par(mar=c(2,2.5,2,2), oma=c(3,3,1,1))


#Plot data
for(usetr in 1:2) {
  require(lmodel2)
  
  #lists for storing r squred values
  rsq_outlst_log<-NULL
  rsq_outlst_log_tot<-NULL
  
  #load data from simulation
  datout<-datoutlst[[usetr]]
  ssout<-ssoutlst[[usetr]]
  
  #get planted species richness for each plot
  datout$plantedsr<-pltloc$plantedsr[match(datout$plt, pltloc$plot)]
  ssout$plantedsr<-pltloc$plantedsr[match(ssout$plot, pltloc$plot)]
  
  ##########################
  #species-level biomass across all plots
  ##########################
  subs<-which(datout$sr>1) #plot only for non-monocultures
  subs2<-which(is.finite(log(datout[subs,]$est)) & is.finite(log(datout[subs,]$obs))) #only plot non-zero values
  
  plot(log10(obs)~log10(est), datout[subs,][subs2,],
       type="n",
       axes=F, xlab="", ylab="",
       cex.axis=1.3)
  if(usetr==2) {
    put.fig.letter(label = "C.", location = "topleft", cex=2, x=0.04, y=0.98)
  } else {
    put.fig.letter(label = "A.", location = "topleft", cex=2, x=0.04, y=0.98)
  }
  sq<-c(0.1, 1, 10, 100)
  axis(1, log10(sq), sq, cex.axis=1.3)
  axis(2, log10(sq), sq, cex.axis=1.3)
  box()
  abline(a=0, b=1, col="darkgrey", lwd=1.5, lty=3)
  
  #fit ranged major axis regression
  obs<-log10(datout[subs,]$obs)
  est<-log10(datout[subs,]$est)
  
  #for calculating fits in log space, replace zeros with lower 2.5% quantile.
  obs[!is.finite(obs)]<-quantile(obs[is.finite(obs)], 0.025)
  est[!is.finite(est)]<-quantile(est[is.finite(est)], 0.025)
  
  lmd<-lmodel2(obs~est, range.y="interval", range.x="interval", nperm=nrep)
  abline(a=lmd$regression.results[4,2], b=lmd$regression.results[4,3], lwd=1.5)
  
  #add 95% CI for major axis regression slope
  abline(a=mean(obs)-mean(est)*lmd$confidence.intervals[4,5],
         b=lmd$confidence.intervals[4,5], lwd=1.5, lty=2)
  abline(a=mean(obs)-mean(est)*lmd$confidence.intervals[4,4],
         b=lmd$confidence.intervals[4,4], lwd=1.5, lty=2)
  
  points(log10(datout[subs,]$est[subs2]), log10(datout[subs,]$obs[subs2]),
         pch=16,
         col=grey.colors(4, start=0, end=1)[as.numeric(as.factor(datout[subs,]$plantedsr))],
         lwd=1.5, cex=0.8)
  points(log10(datout[subs,]$est[subs2]), log10(datout[subs,]$obs[subs2]),
         pch=1,
         col=1,
         lwd=1, cex=0.8)
  
  #plot R-squared from major axis regression
  rsqest<-lmd$rsquare
  
  r<-"R"
  p<-"p"
  rd<-paste(" =", round(rsqest, 2))
  pv<-paste(" <", ceiling(lmd$regression.results[4,5]*1000)/1000)
  
  z<-log10(datout[subs,]$obs[subs2])
  zfit<-log10(datout[subs,]$est[subs2])
  text((max(zfit)-min(zfit))*0.79+min(zfit), (max(z)-min(z))*0.12+min(z),
       bquote(.(r[1])^2 ~ .(rd[1])), cex=1.3, pos=4)
  text((max(zfit)-min(zfit))*0.79+min(zfit), (max(z)-min(z))*0.02+min(z),
       bquote(.(p[1]) ~ .(pv[1])), cex=1.3, pos=4)
  
  #calculate CD based on observed and predicted values
  if(bootr2) {
    #calculate fit only for non-monocultures
    subs2<-which(ssout$plot%in%pltloc$plot[pltloc$plantedsr>1])
    
    agg_res<-ssout[subs2,]
    sp_pltr<-paste(agg_res$sp, agg_res$plantedsr)
    usp_pltr<-sort(unique(sp_pltr))
    
    #for calculating fits in log space, replace zeros with lower 2.5% quantile.
    agg_res$obs[!is.finite(log(agg_res$obs))]<-quantile(agg_res$obs[is.finite(log(agg_res$obs))], 0.025)
    agg_res$est[!is.finite(log(agg_res$est))]<-quantile(agg_res$est[is.finite(log(agg_res$est))], 0.025)
    
    
    #calculate residuals
    agg_res$res<-(log10(agg_res$obs)-log10(agg_res$est))^2
    #calcualte deviation from mean
    sstot_res<-aggregate(cbind(res=agg_res$obs), list(iter=agg_res$iter), function(x) sum((log10(x)-mean(log10(x), na.rm=T))^2,na.rm=T))
    
    #Null model: B*mono*(1/n)
    #agg_res$null<-e120_abv[match(agg_res$sp, splst)]/agg_res$plantedsr
    #agg_res$null_dev<-(agg_res$obs-agg_res$null)^2
    #sstot_res<-aggregate(cbind(res=agg_res$null_dev), list(iter=agg_res$iter), function(x) sum(x,na.rm=T))
    #agg_res<-aggregate(cbind(res=agg_res$res, sstotrest=agg_res$sstotrest), list(iter=agg_res$iter), function(x) sum(x,na.rm=T))
    
    #CD from mean
    rsest<-1-agg_res$res/sstot_res$res
    rsq_out<-quantile(rsest, c(0.025, pnorm(-1, 0, 1), 0.5, pnorm(1, 0, 1), 0.975), na.rm=T) #relative to mean mass
 
    #save outputs
    rsq_outlst_log<-rbind(rsq_outlst_log,
                          data.frame(rsq=c(unname(rsq_out), mean(rsest)),
                                     level=c(names(rsq_out), "mu"),
                                     sp="plottot"))
    
    rsq_outlst_log_tot<-rbind(rsq_outlst_log_tot,
                              data.frame(rsq=rsest,
                                         sp=99))
    #"plottot" and "99" indicate saved CD for species-level biomass
  }
  
  ##########################
  #total plot biomass
  ##########################
  #aggregate total biomass across iterations for plotting
  agdat_plot<-aggregate(cbind(est=datout$est, obs=datout$obs),
                        list(plantedsr=datout$plantedsr, plt=datout$plt),
                        function(x) sum(x, na.rm=T))
  agdat_plot<-agdat_plot[agdat_plot$plantedsr>1,] #plot only for non-monocultures
  
  #get mean plot-level biomass by planted richness treatment
  agdat_plot_mu<-aggregate(cbind(est=agdat_plot$est, obs=agdat_plot$obs),
                           list(plantedsr=agdat_plot$plantedsr),
                           hurdlemodel)
  agdat_plot_sd<-aggregate(cbind(est=agdat_plot$est, obs=agdat_plot$obs),
                           list(plantedsr=agdat_plot$plantedsr),
                           function(x) sd(log(x), na.rm=T))
  
  #only plot for model with tradeoff
  plot(log10(obs)~log10(est), agdat_plot,
       type="n",
       xlab="", ylab="",
       axes=F,
       cex.axis=1.3)
  if(usetr==2) {
    put.fig.letter(label = "D.", location = "topleft", cex=2, x=0.04, y=0.98)
  } else {
    put.fig.letter(label = "B.", location = "topleft", cex=2, x=0.04, y=0.98)
  }
  
  sq<-c(25, 50, 100, 200, 400)
  axis(1, log10(sq), sq, cex.axis=1.3)
  axis(2, log10(sq), sq, cex.axis=1.3)
  box()
  
  abline(a=0, b=1, col="darkgrey", lwd=1.5, lty=3)
  
  #fit ranged major axis regression
  obs<-log10(agdat_plot$obs)
  est<-log10(agdat_plot$est)
  
  #for calculating fits in log space, replace zeros with lower 2.5% quantile.
  obs[!is.finite(obs)]<-quantile(obs[is.finite(obs)], 0.025)
  est[!is.finite(est)]<-quantile(est[is.finite(est)], 0.025)
  
  lmd<-lmodel2(obs~est, range.y="interval", range.x="interval", nperm=nrep)
  abline(a=lmd$regression.results[4,2], b=lmd$regression.results[4,3], lwd=1.5)
  
  #add 95% CI for major axis regression slope
  abline(a=mean(obs)-mean(est)*lmd$confidence.intervals[4,5],
         b=lmd$confidence.intervals[4,5], lwd=1.5, lty=2)
  abline(a=mean(obs)-mean(est)*lmd$confidence.intervals[4,4],
         b=lmd$confidence.intervals[4,4], lwd=1.5, lty=2)
  
  with(agdat_plot, points(log10(est), log10(obs),
                          pch=16,
                          col=grey.colors(4, start=0, end=1)[as.numeric(as.factor(plantedsr))],
                          lwd=1.5, cex=0.8))
  with(agdat_plot, points(log10(est), log10(obs),
                          pch=1,
                          col=1,
                          lwd=1, cex=0.8))
  
  #plot R-squared from major axis regression
  rsqest<-lmd$rsquare
  
  r<-"R"
  p<-"p"
  rd<-paste(" =", round(rsqest, 2))
  pv<-paste(" <", ceiling(lmd$regression.results[4,5]*1000)/1000)
  
  z<-log10(agdat_plot$obs)
  zfit<-log10(agdat_plot$est)
  text((max(zfit)-min(zfit))*0.79+min(zfit), (max(z)-min(z))*0.12+min(z),
       bquote(.(r[1])^2 ~ .(rd[1])), cex=1.3, pos=4)
  text((max(zfit)-min(zfit))*0.79+min(zfit), (max(z)-min(z))*0.02+min(z),
       bquote(.(p[1]) ~ .(pv[1])), cex=1.3, pos=4)
  
  #calculate CD based on observed and predicted values
  if(bootr2) {
    subs2<-which(ssout$plot%in%pltloc$plot[pltloc$plantedsr>1])
    
    agg_res<-with(ssout[subs2,], aggregate(cbind(obs, est), list(plot=plot, iter=iter, plantedsr=plantedsr), function(x) sum(x, na.rm=T)))
    agg_res$res<-(log10(agg_res$obs)-log10(agg_res$est))^2
    
    for(ii in c(2,4,8,16)) {
      #record fit for each diversity level
      #deviation from predictions
      ssres<-with(agg_res[agg_res$plantedsr==ii,], tapply(res, iter, sum))
      #deviation from the mean
      sstot2<-sum((log10(agdat_plot[agdat_plot$plantedsr==ii,]$obs)-mean(log10(agdat_plot$obs), na.rm=T))^2, na.rm=T)
      
      #null model: sum(B_obs)|D
      #sstot2<-sum((log10(agdat_plot[agdat_plot$plantedsr==ii,]$obs)-mean(log10(agdat_plot[agdat_plot$plantedsr==ii,]$obs), na.rm=T))^2, na.rm=T)
      
      #calculate CD for predictions
      rsqest<-1-ssres/sstot2
      
      rsq_out<-quantile(rsqest, c(0.025, pnorm(-1, 0, 1), 0.5, pnorm(1, 0, 1), 0.975), na.rm=T)
      
      #save outputs
      rsq_outlst_log<-rbind(rsq_outlst_log,
                            data.frame(rsq=c(unname(rsq_out), mean(rsqest)),
                                       level=c(names(rsq_out), "mu"),
                                       sp=paste("bmtot", ii)))
      
      rsq_outlst_log_tot<-rbind(rsq_outlst_log_tot,
                                data.frame(rsq=rsqest,
                                           sp=paste(999, ii)))
      
      #"bmtot" and "999" indicate saved R2 for total plot-level biomass
    }
    
    #record total fit across all diversity levels
    #deviation from predictions
    ssres<-with(agg_res, tapply(res, iter, sum))
    #deviation from mean
    sstot2<-sum((log10(agdat_plot$obs)-mean(log10(agdat_plot$obs), na.rm=T))^2, na.rm=T)
    
    #calculate R2 for predictions and null
    rsqest<-1-ssres/sstot2
    
    #record output
    rsq_out<-quantile(rsqest, c(0.025, pnorm(-1, 0, 1), 0.5, pnorm(1, 0, 1), 0.975), na.rm=T)
    rsq_outlst_log<-rbind(rsq_outlst_log,
                          data.frame(rsq=c(unname(rsq_out), mean(rsqest)),
                                     level=c(names(rsq_out), "mu"),
                                     sp=paste("bmtotall")))
    
    rsq_outlst_log_tot<-rbind(rsq_outlst_log_tot,
                              data.frame(rsq=rsqest,
                                         sp=999))
  }
  
  ##########################
  #species-level biomass across each diversity level
  ##########################
  for(i in c(2,4,8,16)) {
    if(bootr2) {
      #record CD
      subs2<-which(ssout$plot%in%pltloc$plot[pltloc$plantedsr==i])
      
      agg_res<-ssout[subs2,]
      
      #for calculating fits in log space, replace zeros with lower 2.5% quantile.
      agg_res$obs[!is.finite(log(agg_res$obs))]<-quantile(agg_res$obs[is.finite(log(agg_res$obs))], 0.025)
      agg_res$est[!is.finite(log(agg_res$est))]<-quantile(agg_res$est[is.finite(log(agg_res$est))], 0.025)
      
      #deviation from predictions
      agg_res$res<-(log10(agg_res$obs)-log10(agg_res$est))^2
      #deviations from mean
      sstot_res<-aggregate(cbind(res=agg_res$obs), list(iter=agg_res$iter), function(x) sum((log10(x)-mean(log10(x), na.rm=T))^2,na.rm=T))
      #aggregate acrodd iterations
      agg_res<-aggregate(cbind(res=agg_res$res, restot=agg_res$restot), list(iter=agg_res$iter), function(x) sum(x,na.rm=T))
      
      #R2 for predictions
      rsqest<-1-agg_res$res/sstot_res$res
      rsq_out<-quantile(rsqest, c(0.025, pnorm(-1, 0, 1), 0.5, pnorm(1, 0, 1), 0.975), na.rm=T)
      
      rsq_outlst_log<-rbind(rsq_outlst_log,
                            data.frame(rsq=c(unname(rsq_out), mean(rsqest)),
                                       level=c(names(rsq_out), "mu"),
                                       sp=as.character(i)))
      
      rsq_outlst_log_tot<-rbind(rsq_outlst_log_tot,
                                data.frame(rsq=1-agg_res$res/sstot_res$res,
                                           sp=as.character(i)))
    }
  }
  
  if(bootr2) {
    write.csv(rsq_outlst_log, paste("data/data_products/rsqlst_obs_tr_", usetr, "_log.csv", sep=""), row.names=F)
    write.csv(rsq_outlst_log_tot, paste("data/data_products/rsqlst_obs_tot_tr_", usetr, "_log.csv", sep=""), row.names=F)
  }
}