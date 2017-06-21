m<-rbind(c(1,2),
         c(3,4))
layout(m)
ofs2<-c(0.078, -0.01)

par(mar=c(4,4,1,1), oma=c(1,1,1,0))
spcoexdat<-spcoexdatsave
spnumnames<-c("2 Species Mixtures",
              "4 Species Mixtures",
              "8 Species Mixtures",
              "16 Species Mixtures")

dsq<-c(0, 0.02, -0.02)

collst_tr2<-adjustcolor(c("black", "red", "blue"), alpha.f = 0.9)


for(i in 1:length(sprichlst)) {
  sbs<-which(spcoexdat$plsr==sprichlst[i])
  
  plot(1-pzero_notr.3~I(qB_notr+dsq[1]), data=spcoexdat[sbs,][order(spcoexdat[sbs,]$qB_notr),], type="b", col=collst_tr[2],
       xlab="", ylab="", ylim=c(0,1), lwd=1.5, cex=0.5)
  
  segments(spcoexdat[sbs,]$qB_notr+dsq[1],
           1-spcoexdat[sbs,]$pzero_notr.1,
           spcoexdat[sbs,]$qB_notr+dsq[1],
           1-spcoexdat[sbs,]$pzero_notr.5, lwd=1.2, lend=2, col=collst_tr[2])
  
  segments(spcoexdat[sbs,]$qB_notr+dsq[1],
           1-spcoexdat[sbs,]$pzero_notr.2,
           spcoexdat[sbs,]$qB_notr+dsq[1],
           1-spcoexdat[sbs,]$pzero_notr.4, lwd=2.2, lend=2, col=collst_tr2[2])
  
  lines(1-pzero_tr.3~I(qB_notr+dsq[2]), data=spcoexdat[sbs,][order(spcoexdat[sbs,]$qB_tr),], type="b", col=collst_tr[3], lwd=1.5, cex=0.5)
  
  segments(spcoexdat[sbs,]$qB_notr+dsq[2],
           1-spcoexdat[sbs,]$pzero_tr.1,
           spcoexdat[sbs,]$qB_notr+dsq[2],
           1-spcoexdat[sbs,]$pzero_tr.5, lwd=1.2, lend=2, col=collst_tr[3])
  
  segments(spcoexdat[sbs,]$qB_notr+dsq[2],
           1-spcoexdat[sbs,]$pzero_tr.2,
           spcoexdat[sbs,]$qB_notr+dsq[2],
           1-spcoexdat[sbs,]$pzero_tr.4, lwd=2.2, lend=2, col=collst_tr2[3])
  
  
  
  lines(pzero_obs.3~I(qB_notr+dsq[3]), data=spcoexdat[sbs,][order(spcoexdat[sbs,]$qB_tr),], type="b", col=collst_tr[1], lwd=1.5, cex=0.5)
  
  segments(spcoexdat[sbs,]$qB_notr+dsq[3],
           spcoexdat[sbs,]$pzero_obs.1,
           spcoexdat[sbs,]$qB_notr+dsq[3],
           spcoexdat[sbs,]$pzero_obs.5, lwd=1.2, lend=2, col=collst_tr[1])
  
  segments(spcoexdat[sbs,]$qB_notr+dsq[3],
           spcoexdat[sbs,]$pzero_obs.2,
           spcoexdat[sbs,]$qB_notr+dsq[3],
           spcoexdat[sbs,]$pzero_obs.4, lwd=2.2, lend=2, col=collst_tr2[1])
  
  text((max(spcoexdat[sbs,]$qB_notr)-min(spcoexdat[sbs,]$qB_notr))*0.7+min(spcoexdat[sbs,]$qB_notr),
       0.2,
       spnumnames[i], cex=1.3)
  put.fig.letter(paste(LETTERS[i], ".", sep=""), "topleft", offset=ofs2, cex=1.4)
  
  text(spcoexdat[sbs,]$qB_notr, 1-spcoexdat[sbs,]$pzero_notr.5, spcoexdat[sbs,]$sp, xpd=NA, srt=90, cex=0.8, adj = c(1.2))
}

mtext(expression(paste(italic("qB"), "*"[italic(mono)], ", g m"^-1, sep="")), 1, cex=1.2, outer=T, line=-0.5)
mtext("Frequency of Stable Coexistence", 2, cex=1.2, line=-0.5, outer=T)
