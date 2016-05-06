#plot blank space
plot(0,0,type="n", axes=F, xlab="", ylab=""); plot(0,0,type="n", axes=F, xlab="", ylab="")


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
    trnotr<-t.test(log(tmp1), log(tmp2), paired = TRUE)$p.value
    
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
  put.fig.letter(label = "E.", location = "topleft", cex=2, x=0.04, y=0.98)
  
  sq<-log(seq(0, 6, by=c(0.5, 0.1)[itype]))
  abline(h=sq[-which(sq==0)], col="grey")
  abline(v=seq(1, 5, by=1), col="grey")
  
  axis(2, sq, round(exp(sq),1)-1, cex.axis=1.3)
  
  axis(1, 1:5, c("2", "4", "8", "16", "2-16"), cex.axis=1.3)
  text(1:5+c(rep(0.18, 3), 0.22, 0.4), rep(ylims1-(ylims2-ylims1)*0.125, 5), pvl[1:5], xpd=NA)
  box()
  
  adj<-c(-0.2, 0.1, -0.1, 0.2)
  collst<-c("dimgrey", "black")
  
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
  sppos<-2
  segments(1:6+adj[sppos], rsuse$l025, 1:6+adj[sppos], rsuse$u975, col=collst[colpos], lwd=1.5, lend=2)
  segments(1:6+adj[sppos], rsuse$mumsd, 1:6+adj[sppos], rsuse$mupsd, col=collst[colpos], lwd=3, lend=2)
  segments(1:6+adj[sppos]-0.05, rsuse$mu, 1:6+adj[sppos]+0.05, rsuse$mu, col=collst[colpos], lwd=3, lend=2)
  
  lines(1:4+adj[sppos], rsuse$mu[1:4], col=collst[colpos], lwd=1.5, lty=1)
}


#make figure labels
par(mfrow=c(1,1), new=TRUE)
mtext("Planted Richness", 1, line=2, cex=1.5, outer=T)
mtext(expression(paste("Observed Biomass, g m"^"-2", sep="")), 2, line=1, cex=1.5, outer=T, adj=0.8)
mtext(expression(paste("MAE, fold change", sep="")), 2, line=1, cex=1.5, outer=T, adj=0)  

put.fig.letter(expression(paste("Predicted Biomass, g m"^"-2", sep="")), x=0.5, y=4.5/13, cex=1.5)
put.fig.letter("Species-Level Biomass", x=0.25, y=1, cex=1.3)
put.fig.letter("Total Plot Biomass", x=0.75, y=1, cex=1.3)

put.fig.letter("With Snapping", x=1, y=7/13, cex=1.1, srt=90+180)
put.fig.letter("Without Snapping", x=1, y=11/13, cex=1.1, srt=90+180)