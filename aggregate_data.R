#Set up plot
m<-rbind(c(5,5,5,5,6,6,6,6),
         c(5,5,5,5,6,6,6,6),
         c(5,5,5,5,6,6,6,6),
         c(5,5,5,5,6,6,6,6),
         c(7,7,7,7,8,8,8,8),
         c(1,1,1,1,2,2,2,2),
         c(1,1,1,1,2,2,2,2),
         c(1,1,1,1,2,2,2,2),
         c(1,1,1,1,2,2,2,2),
         c(3,3,3,3,4,4,4,4),
         c(3,3,3,3,4,4,4,4),
         c(3,3,3,3,4,4,4,4),
         c(3,3,3,3,4,4,4,4))

layout(m)
replacelower<-log10(0.02) #0.0025 quantile - 99.9% interval

par(mar=c(2,2.5,2,2), oma=c(3,3,1,0))


#Plot data
for(usetr in 1:2) {
  require(lmodel2)
  
  #lists for storing r squred values
  rsq_outlst_log<-NULL
  rsq_outlst_log_tot<-NULL
  
  #load data from simulation
  datout<-datoutlst[[usetr]]
  
  #get planted species richness for each plot
  datout$plantedsr<-pltloc$plantedsr[match(datout$plt, pltloc$plot)]
  datout_old<-datout
  
  datout$obs[!is.finite(log(datout$obs)) | log10(datout$obs)<replacelower]<-10^replacelower
  datout$est[!is.finite(log(datout$est)) | log10(datout$est)<replacelower]<-10^replacelower
  
  ##########################
  #species-level biomass across all plots
  ##########################
  subs<-which(datout$sr>1) #plot only for non-monocultures
  
  plot(log10(obs)~log10(est), datout[subs,],
       type="n",
       axes=F, xlab="", ylab="",
       cex.axis=1.3,
       xlim=c(replacelower, log10(max(datout[subs,]$est))),
       ylim=c(replacelower, log10(max(datout[subs,]$obs))))
  text(replacelower, log10(max(datout[subs,]$obs)),
       c("without snapping", "with snapping")[usetr], cex=1.8, adj=c(0,1))
  
  if(usetr==2) {
    put.fig.letter(label = "E.", location = "topleft", cex=2, x=0.04, y=0.98)
  } else {
    put.fig.letter(label = "C.", location = "topleft", cex=2, x=0.04, y=0.98)
  }
  
  sq0<-c(0.1, 1, 10, 100)
  sq0<-sq0[sq0>10^replacelower]
  sq<-c(10^replacelower, sq0)
  sq_name<-c(paste("<", 10^replacelower), sq0)
  axis(1, log10(sq), sq_name, cex.axis=1.3)
  axis(2, log10(sq), sq_name, cex.axis=1.3)
  box()
  abline(a=0, b=1, col="darkgrey", lwd=1.5, lty=3)
  abline(h=replacelower, v=replacelower, col="darkgrey", lwd=1.5, lty=3)
  
  #fit ranged major axis regression
  obs<-log10(datout[subs,]$obs)
  est<-log10(datout[subs,]$est)
  
  lmd<-lmodel2(obs~est, range.y="interval", range.x="interval", nperm=nrep)
  abline(a=lmd$regression.results[4,2], b=lmd$regression.results[4,3], lwd=1.5)
  
  #add 95% CI for major axis regression slope
  abline(a=mean(obs)-mean(est)*lmd$confidence.intervals[4,5],
         b=lmd$confidence.intervals[4,5], lwd=1.5, lty=2)
  abline(a=mean(obs)-mean(est)*lmd$confidence.intervals[4,4],
         b=lmd$confidence.intervals[4,4], lwd=1.5, lty=2)
  
  #add points
  obs_plot<-log10(datout[subs,]$obs)
  est_plot<-log10(datout[subs,]$est)
  
  points(est_plot, obs_plot,
         pch=16,
         col=grey.colors(4, start=0, end=1)[as.numeric(as.factor(datout[subs,]$plantedsr))],
         lwd=1.5, cex=0.8)
  points(est_plot, obs_plot,
         pch=1,
         col=1,
         lwd=1, cex=0.8)
  
  #plot R-squared from major axis regression
  rsqest<-lmd$rsquare
  
  r<-"R"
  p<-"p"
  rd<-paste(" =", round(rsqest, 2))
  pv<-paste(" <", ceiling(lmd$regression.results[4,5]*1000)/1000)
  
  z<-log10(c(10^replacelower, datout[subs,]$obs))
  zfit<-log10(c(10^replacelower, datout[subs,]$est))
  text((max(zfit)-min(zfit))*0.79+min(zfit), (max(z)-min(z))*0.12+min(z),
       bquote(.(r[1])^2 ~ .(rd[1])), cex=1.3, pos=4)
  text((max(zfit)-min(zfit))*0.79+min(zfit), (max(z)-min(z))*0.02+min(z),
       bquote(.(p[1]) ~ .(pv[1])), cex=1.3, pos=4)
  
  #calculate MAE based on observed and predicted values
  if(bootr2) {
    #calculate fit only for non-monocultures
    subs2<-which(datout$plt%in%pltloc$plot[pltloc$plantedsr>1])
    
    agg_res<-datout[subs2,]
    sp_pltr<-paste(agg_res$sp, agg_res$plantedsr)
    usp_pltr<-sort(unique(sp_pltr))
        
    #calculate residuals
    agg_res$res<-abs(log10(agg_res$obs)-log10(agg_res$est))
        
    #MAE
    #mean within plots (across all species)
    agg_res<-with(agg_res, aggregate(cbind(res=res), list(plot=plt, plantedsr=plantedsr), function(x) mean(x,na.rm=T)))
    
    mae_out<-t(matrix(nrow=6, unlist(tapply(agg_res$res, agg_res$plantedsr, function(x) {
      mux<-mean(x)
      mulogx<-mean(log(x))
      sdx<-sd(log(x))/sqrt(length(x))
      c(mux, mulogx, exp(log(mux)-sdx), exp(log(mux)+sdx), exp(log(mux)-sdx*qnorm(0.975, 0, 1)), exp(log(mux)+sdx*qnorm(0.975, 0, 1)))
    }))))
    
    mux<-mean(agg_res$res)
    mulogx<-mean(log(agg_res$res))
    sdx<-sd(log(agg_res$res))/sqrt(length(agg_res$res))
    
    mae_out<-rbind(mae_out, c(mux, mulogx, exp(log(mux)-sdx), exp(log(mux)+sdx), exp(log(mux)-sdx*qnorm(0.975, 0, 1)), exp(log(mux)+sdx*qnorm(0.975, 0, 1))))
    colnames(mae_out)<-c("mu", "mulogx", "mumsd", "mupsd", "l025", "u975")
    
    #save outputs
    rsq_outlst_log<-rbind(rsq_outlst_log,
                          data.frame(mae_out,
                                     level=c(2,4,8,16, "mu"),
                                     type="species"))
    
    rsq_outlst_log_tot<-rbind(rsq_outlst_log_tot,
                              data.frame(mae=agg_res$res,
                                         level=agg_res$plantedsr,
                                         type="species"))
  }
  
  ##########################
  #total plot biomass
  ##########################
  datout<-datout_old #remove "fixed" biomasses
  
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
  
  #plot
  plot(log10(obs)~log10(est), agdat_plot,
       type="n",
       xlab="", ylab="",
       axes=F,
       cex.axis=1.3)
  text(min(log10(agdat_plot$est)), max(log10(agdat_plot$obs)),
       c("without snapping", "with snapping")[usetr], cex=1.8, adj=c(0,1))
  
  if(usetr==2) {
    put.fig.letter(label = "F.", location = "topleft", cex=2, x=0.04, y=0.98)
  } else {
    put.fig.letter(label = "D.", location = "topleft", cex=2, x=0.04, y=0.98)
  }
  
  sq<-c(25, 50, 100, 200, 400)
  axis(1, log10(sq), sq, cex.axis=1.3)
  axis(2, log10(sq), sq, cex.axis=1.3)
  box()
  
  abline(a=0, b=1, col="darkgrey", lwd=1.5, lty=3)
  
  #fit ranged major axis regression
  obs<-log10(agdat_plot$obs)
  est<-log10(agdat_plot$est)
    
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
  
  #calculate MAE based on observed and predicted values
  if(bootr2) {
    #calculate fit only for non-monocultures
    subs2<-which(agdat_plot$plt%in%pltloc$plot[pltloc$plantedsr>1])
    
    agg_res<-agdat_plot[subs2,]
    sp_pltr<-paste(agg_res$sp, agg_res$plantedsr)
    usp_pltr<-sort(unique(sp_pltr))
        
    #calculate residuals
    agg_res$res<-abs(log10(agg_res$obs)-log10(agg_res$est))
    agg_res$plot<-agg_res$plt
    
    mae_out<-t(matrix(nrow=6, unlist(tapply(agg_res$res, agg_res$plantedsr, function(x) {
      mux<-mean(x)
      mulogx<-mean(log(x))
      sdx<-sd(log(x))/sqrt(length(x))
      c(mux, mulogx, exp(log(mux)-sdx), exp(log(mux)+sdx), exp(log(mux)-sdx*qnorm(0.975, 0, 1)), exp(log(mux)+sdx*qnorm(0.975, 0, 1)))
    }))))
    
    mux<-mean(agg_res$res)
    mulogx<-mean(log(agg_res$res))
    sdx<-sd(log(agg_res$res))/sqrt(length(agg_res$res))
    
    mae_out<-rbind(mae_out, c(mux, mulogx, exp(log(mux)-sdx), exp(log(mux)+sdx), exp(log(mux)-sdx*qnorm(0.975, 0, 1)), exp(log(mux)+sdx*qnorm(0.975, 0, 1))))
    colnames(mae_out)<-c("mu", "mulogx", "mumsd", "mupsd", "l025", "u975")
    
    #save outputs
    rsq_outlst_log<-rbind(rsq_outlst_log,
                          data.frame(mae_out,
                                     level=c(2,4,8,16, "mu"),
                                     type="totalnpp"))
    
    rsq_outlst_log_tot<-rbind(rsq_outlst_log_tot,
                              data.frame(mae=agg_res$res,
                                         level=agg_res$plantedsr,
                                         type="totalnpp"))
  }
  
  if(bootr2) {
    write.csv(rsq_outlst_log, paste("data/data_products/rsqlst_obs_tr_", usetr, "_log.csv", sep=""), row.names=F)
    write.csv(rsq_outlst_log_tot, paste("data/data_products/rsqlst_obs_tot_tr_", usetr, "_log.csv", sep=""), row.names=F)
  }
}