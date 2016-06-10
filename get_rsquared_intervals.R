#load data
pnums<-colSums(table(e120dat$Plot[e120dat$Year==2014], e120dat$NumSp[e120dat$Year==2014])>0)
rsqlst_1<-read.csv("data/data_products/rsqlst_obs_tr_1_log.csv")
rsqlst_2<-read.csv("data/data_products/rsqlst_obs_tr_2_log.csv")

rsqlst_1_tot<-read.csv("data/data_products/rsqlst_obs_tot_tr_1_log.csv")
rsqlst_2_tot<-read.csv("data/data_products/rsqlst_obs_tot_tr_2_log.csv")
  
#levels for significance tests
pvallevels<-c(0.1, 0.05, 0.01, 0.001)
pvalmark<-c("", "", "*", "**", "***")

#Run significance tests between snapping and no snapping
#levels indicated in R2 tables
spllvls<-c(2,4,8,16,"mu")
#list for outputting significance tests

for(itype in 1:2) {
  type<-c("species", "totalnpp")[itype]
  pvallst<-NULL
  
  for(i in 1:length(spllvls)) {
    if(spllvls[i]!="mu") {
      tmp1<-unlist(rsqlst_1_tot[(rsqlst_1_tot$level==spllvls[i]) & rsqlst_1_tot$type==type,]$mae)
      tmp2<-unlist(rsqlst_2_tot[(rsqlst_2_tot$level==spllvls[i]) & rsqlst_2_tot$type==type,]$mae)
    } else {
      tmp1<-unlist(rsqlst_1_tot[rsqlst_1_tot$type==type,]$mae)
      tmp2<-unlist(rsqlst_2_tot[rsqlst_2_tot$type==type,]$mae)
    }
    
    #snapping vs. no snapping
    trnotr<-suppressWarnings(wilcox.test((tmp1), (tmp2), paired = TRUE)$p.value)
    #trnotr<-t.test((tmp1), (tmp2), paired = TRUE)$p.value
    
    pvallst<-rbind(pvallst, data.frame(trnotr=trnotr, sp=spllvls[i]))
  }
  
  #extract p-vales
  pvl<-rep("", 5)
  for(i in 1:5) {
    pvl[i]<-pvalmark[sum(pvallst[i,"trnotr"]<pvallevels)+1]
  }
  
  #make plot
  ylims2<-max(c(rsqlst_1[rsqlst_1$type==type,]$u975, rsqlst_2[rsqlst_1$type==type,]$u975))
  ylims1<-min(c(rsqlst_1[rsqlst_1$type==type,]$l025, rsqlst_2[rsqlst_1$type==type,]$l025))
  plot(c(0.5,5.5), c(ylims1,ylims2), axes=F, type="n", xlab="", ylab="",
       cex.axis=1.3)
  put.fig.letter(label = c("C.", "D.")[itype], location = "topleft", cex=2, x=0.04, y=0.98)
  
  if(itype==1) {
    sq<-log10(c(1,2,3,5,9))
  } else {
    sq<-log10(seq(0.1, 12, by=c(0.1)))
  }
  abline(h=sq[-which(sq==0)], col="lightgrey")
  abline(v=seq(0.5, 5.5, by=1), col="lightgrey")
  
  axis(2, sq, round(10^(sq),1)-1, cex.axis=1.3)
  
  axis(1, 1:5, c("2", "4", "8", "16", "2-16"), cex.axis=1.3)
  text(1:5+c(rep(0.18, 3), 0.22, 0.4), rep(ylims1-(ylims2-ylims1)*0.11, 5), pvl[1:5], xpd=NA)
  box()
  
  adj<-c(-0.2, 0.1, -0.1, 0.2)
  collst<-c("darkgrey", "black")
  
  #No snapping
  rsuse<-rsqlst_1[rsqlst_1$type==type,]
  colpos<-1
  sppos<-1
  segments(1:6+adj[sppos], rsuse$l025, 1:6+adj[sppos], rsuse$u975, col=collst[colpos], lwd=1.5, lend=2)
  segments(1:6+adj[sppos], rsuse$mumsd, 1:6+adj[sppos], rsuse$mupsd, col=collst[colpos], lwd=3, lend=2)
  segments(1:6+adj[sppos]-0.05, rsuse$mu, 1:6+adj[sppos]+0.05, rsuse$mu, col=collst[colpos], lwd=3, lend=2)
  
  lines(1:4+adj[sppos], rsuse$mu[1:4], col=collst[colpos], lwd=1.5, lty=1)
  
  #snapping
  rsuse<-rsqlst_2[rsqlst_1$type==type,]
  colpos<-2
  sppos<-4
  segments(1:6+adj[sppos], rsuse$l025, 1:6+adj[sppos], rsuse$u975, col=collst[colpos], lwd=1.5, lend=2)
  segments(1:6+adj[sppos], rsuse$mumsd, 1:6+adj[sppos], rsuse$mupsd, col=collst[colpos], lwd=3, lend=2)
  segments(1:6+adj[sppos]-0.05, rsuse$mu, 1:6+adj[sppos]+0.05, rsuse$mu, col=collst[colpos], lwd=3, lend=2)
  
  lines(1:4+adj[sppos], rsuse$mu[1:4], col=collst[colpos], lwd=1.5, lty=1)
}
