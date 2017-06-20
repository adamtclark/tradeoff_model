richmat_est_lst<-NULL
richlst<-NULL

#calc richness differences
for(iternumber in 1:4) {
  #load data from simulation
  datout<-datoutlst[[iternumber]]
  
  #zeros
  datout$est[!is.finite(datout$est)]<-0
  datout$obs[!is.finite(datout$obs)]<-0

  #get diversity
  rich<-with(datout[datout$plantedsr>1,], aggregate(cbind(obs, est), list(plt=plt, plantedsr=plantedsr), function(x) sum(x>0)))
  
  plrichlst<-sort(unique(rich$plantedsr))
  
  if(iternumber==1) {
    richmat_obs<-matrix(nrow=4, ncol=5)
  }
  richmat_est<-matrix(nrow=4, ncol=5)

  for(i in 1:4) {
    if(iternumber==1) {
      #observed
      subs<-which(rich$plantedsr==plrichlst[i])
      mux<-exp(mean(log(rich[subs,]$obs)))
      sdx<-(sd(log(rich[subs,]$obs)))
      
      richmat_obs[i,]<-c(mux, exp(log(mux)-sdx), exp(log(mux)+sdx),
                     exp(log(mux)+sdx*qnorm(0.025,0,1)), exp(log(mux)+sdx*qnorm(0.975,0,1)))
    }
    
    #estimated
    subs<-which(rich$plantedsr==plrichlst[i])
    mux<-exp(mean(log(rich[subs,]$est)))
    sdx<-(sd(log(rich[subs,]$est)))
    
    richmat_est[i,]<-c(mux, exp(log(mux)-sdx), exp(log(mux)+sdx),
                       exp(log(mux)+sdx*qnorm(0.025,0,1)), exp(log(mux)+sdx*qnorm(0.975,0,1)))
  }
  
  richlst[[iternumber]]<-rich
  richmat_est_lst[[iternumber]]<-richmat_est
}


#plot
collst<-c("red", "blue")
adj<-c(-0.325, 0.1625, -0.1625, 0.325, 0)
srichlst<-c(2,4,8,16)
pvallevels<-c(0.1, 0.05, 0.01, 0.001)
pvalmark<-c("", "", "*", "**", "***")

ylims1<-min(c(richmat_obs, unlist(richmat_est)))
ylims2<-max(c(richmat_obs, unlist(richmat_est)))

for(plti in 1:2) {
  plot(c(0.5,4.5), c(1,ylims2), axes=F, type="n", xlab="", ylab="",
       cex.axis=1.3,
       main="")
  put.fig.letter(label = c("A.", "B.")[plti], location = "topleft", cex=2, x=0.04, y=0.98)  
  axis(2, cex.axis=1.3)
  axis(1, 1:4, c("2", "4", "8", "16"), cex.axis=1.3)
  abline(h=seq(2, 14, by=2), v=c(0.5, 1.5, 2.5, 3.5, 4.5), col="lightgrey")
  box()
  
  #observed richness
  sppos<-5
  colpos<-2
  
  segments(1:4+adj[sppos], richmat_obs[,4], 1:4+adj[sppos], richmat_obs[,5], lwd=3.5, lend=2, col="black")
  segments(1:4+adj[sppos], richmat_obs[,2], 1:4+adj[sppos], richmat_obs[,3], lwd=5, lend=2, col="black")
  segments(1:4+adj[sppos]-0.05, richmat_obs[,1], 1:4+adj[sppos]+0.05, richmat_obs[,1], lwd=5, lend=2, col="black")
  
  lines(1:4+adj[sppos], richmat_obs[,1], lwd=3.5, col="black")
  
  segments(1:4+adj[sppos], richmat_obs[,4], 1:4+adj[sppos], richmat_obs[,5], lwd=1, lend=2, col="black")
  segments(1:4+adj[sppos], richmat_obs[,2], 1:4+adj[sppos], richmat_obs[,3], lwd=2, lend=2, col="black")
  segments(1:4+adj[sppos]-0.05, richmat_obs[,1], 1:4+adj[sppos]+0.05, richmat_obs[,1], lwd=2, lend=2, col="black")
  
  lines(1:4+adj[sppos], richmat_obs[,1], lwd=1, col="black")
  
  #est
  clst<-rbind(c(3,4), c(1,2))
  for(iternumber in clst[plti,]) {
    richmat_est<-richmat_est_lst[[iternumber]]
    sppos<-c(1,4,3,2)[iternumber]
    colpos<-c(1,2,1,2)[iternumber]
    linepos<-c(1,1,3,3)[iternumber]
    
    segments(1:4+adj[sppos], richmat_est[,2], 1:4+adj[sppos], richmat_est[,3], lwd=3, lend=2, col=collst[colpos], lty=1)
    segments(1:4+adj[sppos], richmat_est[,3], 1:4+adj[sppos], richmat_est[,5], lwd=1.5, lend=2, col=collst[colpos], lty=1)
    segments(1:4+adj[sppos], richmat_est[,2], 1:4+adj[sppos], richmat_est[,4], lwd=1.5, lend=2, col=collst[colpos], lty=1)
    segments(1:4+adj[sppos]-0.05, richmat_est[,1], 1:4+adj[sppos]+0.05, richmat_est[,1], lwd=3, lend=2, col=collst[colpos], lty=1)
    
    #lines(1:4+adj[sppos], richmat_est[,1], lwd=1.5, col=collst[colpos], lty=1)
    
    #p-value
    tmp<-richlst[[iternumber]]
    pv<-numeric(4)
    for(j in 1:4) {
      subs<-which(tmp$plantedsr==srichlst[j])
      pv<-wilcox.test(tmp$obs[subs], tmp$est[subs], paired=TRUE, exact = FALSE)$p.value
      
      text(j+adj[sppos], richmat_est[j,4], pvalmark[sum(pv<pvallevels)+1], col=collst[colpos], adj=c(c(1, 0, 1, 0)[iternumber], 1))
    }
  }
}
