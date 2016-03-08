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


##########################
#species-level tests
##########################
#Run significance tests between snapping and no snapping
#levels indicated in R2 tables
spllvls<-c(2,4,8,16,99)
#list for outputting significance tests
pvallst<-NULL
nperm<-mean(table(rsqlst_1_tot$sp))
for(i in 1:length(spllvls)) {
  #snapping vs. no snapping
  trnotr<-sum(rsqlst_1_tot[rsqlst_1_tot$sp==spllvls[i],]$rsq>
        rsqlst_2_tot[rsqlst_2_tot$sp==spllvls[i],]$rsq)/nperm
  
  pvallst<-rbind(pvallst, data.frame(trnotr=trnotr, sp=spllvls[i]))
}

#extract p-vales
pvl<-rep("", 5)
for(i in 1:5) {
  pvl[i]<-pvalmark[sum(pvallst[i,"trnotr"]<pvallevels)+1]
}

#make plot
plot(c(0.5,5.5), c(-1.35,0.8), axes=F, type="n", xlab="", ylab="",
     cex.axis=1.3)
put.fig.letter(label = "E.", location = "topleft", cex=2, x=0.04, y=0.98)

sq<-round(seq(-1.2, 0.8, by=0.4),1)
abline(h=sq[-which(sq==0)], col="grey")
abline(v=seq(1, 5, by=1), col="grey")

axis(2, sq, sq, cex.axis=1.3)

axis(1, 1:5, c("2", "4", "8", "16", "2-16"), cex.axis=1.3)
text(1:5+c(rep(0.18, 3), 0.22, 0.4), rep(-0.23-1.35, 5), pvl[1:5], xpd=NA)
box()
abline(h=c(0), col="darkgrey", lwd=1.5, lty=3)

adj<-c(-0.2, 0.1, -0.1, 0.2)
psadj<-which(rsqlst_1$level=="2.5%")[match(c("2", "4", "8", "16", "plottot", "bmtot"), rsqlst_1$sp[which(rsqlst_1$level=="2.5%")])]
collst<-c("dimgrey", "black")

#No snapping
rsuse<-rsqlst_1
colpos<-1
sppos<-1
segments(1:6+adj[sppos], rsuse[psadj,]$rsq, 1:6+adj[sppos], rsuse[psadj+4,]$rsq, col=collst[colpos], lwd=1.5)
segments(1:6+adj[sppos], rsuse[psadj+1,]$rsq, 1:6+adj[sppos], rsuse[psadj+3,]$rsq, col=collst[colpos], lwd=3)
segments(1:6+adj[sppos]-0.05, rsuse[psadj+2,]$rsq, 1:6+adj[sppos]+0.05, rsuse[psadj+2,]$rsq, col=collst[colpos], lwd=3)

lines(1:4+adj[sppos], rsuse[psadj+2,]$rsq[1:4], col=collst[1], lwd=1.5, lty=1)

#snapping
rsuse<-rsqlst_2
colpos<-2
sppos<-2
segments(1:6+adj[sppos], rsuse[psadj,]$rsq, 1:6+adj[sppos], rsuse[psadj+4,]$rsq, col=collst[colpos], lwd=1.5)
segments(1:6+adj[sppos], rsuse[psadj+1,]$rsq, 1:6+adj[sppos], rsuse[psadj+3,]$rsq, col=collst[colpos], lwd=3)
segments(1:6+adj[sppos]-0.05, rsuse[psadj+2,]$rsq, 1:6+adj[sppos]+0.05, rsuse[psadj+2,]$rsq, col=collst[colpos], lwd=3)

lines(1:4+adj[sppos], rsuse[psadj+2,]$rsq[1:4], col=collst[2], lwd=1.5)
p_r2_species_level<-recordPlot()

##########################
#total plot biomass
##########################
#load correct positions in saved matrix
spllvls<-levels(rsqlst_2_tot$sp)[grep("999", levels(rsqlst_2_tot$sp))]
#output list
pvallst<-NULL
nperm<-mean(table(rsqlst_1_tot$sp))
for(i in 1:length(spllvls)) {
  #snapping vs. no snapping
  trnotr<-sum(rsqlst_1_tot[rsqlst_1_tot$sp==spllvls[i],]$rsq>
                rsqlst_2_tot[rsqlst_2_tot$sp==spllvls[i],]$rsq)/nperm
  
  pvallst<-rbind(pvallst, data.frame(trnotr=trnotr, sp=spllvls[i]))
}
pvallst<-pvallst[order(as.numeric(gsub("999 ", "", as.character(pvallst$sp), fixed=T))),]

#extract p-values
pvl<-rep("", 5)
for(i in 1:5) {
  pvl[i]<-pvalmark[sum(pvallst[i,"trnotr"]<pvallevels)+1]
}

#make plot
plot(c(0.5,5.5), c(-0.25,0.8), axes=F, type="n", xlab="", ylab="",
     cex.axis=1.3)

sq<-round(seq(-0.4, 0.8, by=0.2),2)

put.fig.letter(label = "F.", location = "topleft", cex=2, x=0.04, y=0.98)
abline(h=c(0), col="darkgrey", lwd=1.5, lty=3)
abline(h=sq[-which(sq==0)], col="grey")
abline(v=seq(1, 5, by=1), col="grey")

axis(2, sq, sq, cex.axis=1.3)
axis(1, 1:5, c("2", "4", "8", "16", "2-16"), cex.axis=1.3)
text(1:5+c(rep(0.18, 3), 0.22, 0.4), rep(-0.11-0.25, 5), pvl[1:5], xpd=NA)
box()

adj<-c(-0.2, 0.1, -0.1, 0.2)
psadj<-which(rsqlst_1$level=="2.5%")[match(c("bmtot 2", "bmtot 4", "bmtot 8", "bmtot 16", "bmtotall"), rsqlst_1$sp[which(rsqlst_1$level=="2.5%")])]
collst<-c("dimgrey", "black")

#No snapping
rsuse<-rsqlst_1
colpos<-1
sppos<-1
segments(1:6+adj[sppos], rsuse[psadj,]$rsq, 1:6+adj[sppos], rsuse[psadj+4,]$rsq, col=collst[colpos], lwd=1.5)
segments(1:6+adj[sppos], rsuse[psadj+1,]$rsq, 1:6+adj[sppos], rsuse[psadj+3,]$rsq, col=collst[colpos], lwd=3)
segments(1:6+adj[sppos]-0.05, rsuse[psadj+2,]$rsq, 1:6+adj[sppos]+0.05, rsuse[psadj+2,]$rsq, col=collst[colpos], lwd=3)

lines(1:4+adj[sppos], rsuse[psadj+2,]$rsq[1:4], col=collst[1], lwd=1.5, lty=1)

#snapping
rsuse<-rsqlst_2
colpos<-2
sppos<-2
segments(1:6+adj[sppos], rsuse[psadj,]$rsq, 1:6+adj[sppos], rsuse[psadj+4,]$rsq, col=collst[colpos], lwd=1.5)
segments(1:6+adj[sppos], rsuse[psadj+1,]$rsq, 1:6+adj[sppos], rsuse[psadj+3,]$rsq, col=collst[colpos], lwd=3)
segments(1:6+adj[sppos]-0.05, rsuse[psadj+2,]$rsq, 1:6+adj[sppos]+0.05, rsuse[psadj+2,]$rsq, col=collst[colpos], lwd=3)

lines(1:4+adj[sppos], rsuse[psadj+2,]$rsq[1:4], col=collst[2], lwd=1.5)


#make figure labels
par(mfrow=c(1,1), new=TRUE)
mtext("Planted Richness", 1, line=2, cex=1.5, outer=T)
mtext(expression(paste("Observed Biomass, g m"^"-2", sep="")), 2, line=1, cex=1.5, outer=T, adj=0.8)
mtext(expression(paste("Coef. of Determination (CD)", sep="")), 2, line=1, cex=1.5, outer=T, adj=-0.12)  

put.fig.letter(expression(paste("Predicted Biomass, g m"^"-2", sep="")), x=0.5, y=4.5/13, cex=1.5)
put.fig.letter("Species-Level Biomass", x=0.25, y=1, cex=1.3)
put.fig.letter("Total Plot Biomass", x=0.75, y=1, cex=1.3)

put.fig.letter("With Snapping", x=1, y=7/13, cex=1.1, srt=90+180)
put.fig.letter("Without Snapping", x=1, y=11/13, cex=1.1, srt=90+180)