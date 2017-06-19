#1 is no tradeoff, 2 is tradeoff
sprichlst<-sort(unique(datoutlst[[1]]$plantedsr[datoutlst[[1]]$plantedsr>1]))

bootfrac_fun<-function(tfdat) {
  if(sum(tfdat)==length(tfdat)) {
    return(c(m2sd=1, m1sd=1, mu=1, p1sd=1, p2sd=1))
  } else {
    blst<-numeric(nrep_traits)
    for(ibf in 1:nrep_traits) {
      smp<-sample(1:length(tfdat), rep=T)
      blst[ibf]<-sum(tfdat[smp])
    }
    blst<-blst/length(tfdat)
    qout<-quantile(blst, c(pnorm(-2, 0, 1), pnorm(-1, 0, 1), 0.5, pnorm(1, 0, 1), pnorm(2, 0, 1)))
    names(qout)<-c("m2sd", "m1sd", "mu", "p1sd", "p2sd")
    
    return(qout)
  }
}

bootpobs_fun<-function(pdat) {
  if(sum(pdat==1)==length(pdat)) {
    return(c(m2sd=1, m1sd=1, mu=1, p1sd=1, p2sd=1))
  } else {
    blst<-numeric(nrep_traits)
    for(ibf in 1:nrep_traits) {
      smp<-sample(1:length(pdat), rep=T)
      blst[ibf]<-mean(pdat[smp])
    }
    qout<-quantile(blst, c(pnorm(-2, 0, 1), pnorm(-1, 0, 1), 0.5, pnorm(1, 0, 1), pnorm(2, 0, 1)))
    names(qout)<-c("m2sd", "m1sd", "mu", "p1sd", "p2sd")
    
    return(qout)
  }
}

coexdat<-NULL
spcoexdat<-NULL

for(i in 1:length(sprichlst)) {
  ps_notr<-which(datoutlst[[1]]$plantedsr==sprichlst[i])
  ps_tr<-which(datoutlst[[2]]$plantedsr==sprichlst[i])
  
  #get observed coex frequency
  coex<-!is.na(datoutlst[[1]]$obs[ps_notr])
  #sum(coex)/length(coex)
  fcoex<-bootfrac_fun(coex)
  #fcoex<-bootpobs_fun(1-datoutlst[[1]]$pzero_obs)
  
  #get estimated coex frequency, for both tradeoff and non-tradeoff
  pz_notr<-datoutlst[[1]]$pzero_est[ps_notr]
  pz_tr<-datoutlst[[2]]$pzero_est[ps_tr]
  f_coex_est_notr<-bootpobs_fun(1-pz_notr)
  f_coex_est_tr<-bootpobs_fun(1-pz_tr)
  
  #pz_notr<-!is.na(datoutlst[[1]]$est[ps_notr])
  #pz_tr<-!is.na(datoutlst[[2]]$est[ps_tr])
  #f_coex_est_notr<-bootfrac_fun(pz_notr)
  #f_coex_est_tr<-bootfrac_fun(pz_tr)
  
  coextmp<-data.frame(fcoex, f_coex_est_notr, f_coex_est_tr, plsr=sprichlst[i])
  
  coexdat<-rbind(coexdat, coextmp)
  
  tmp<-data.frame(sp_notr=datoutlst[[1]]$sp[ps_notr],
             qB_notr=datoutlst[[1]]$abv[ps_notr]*datoutlst[[1]]$ptisn[ps_notr],
             Rs_notr=datoutlst[[1]]$no3[ps_notr],
             pzero_obs_notr=datoutlst[[1]]$obs[ps_notr],
             pzero_est_notr=datoutlst[[1]]$pzero_est[ps_notr],
             plantedsr_notr=datoutlst[[1]]$plantedsr[ps_notr],
             sp_tr=datoutlst[[2]]$sp[ps_notr],
             qB_tr=datoutlst[[2]]$abv[ps_notr]*datoutlst[[1]]$ptisn[ps_notr],
             Rs_tr=datoutlst[[2]]$no3[ps_notr],
             pzero_obs_tr=datoutlst[[2]]$obs[ps_notr],
             pzero_est_tr=datoutlst[[2]]$pzero_est[ps_notr],
             plantedsr_tr=datoutlst[[2]]$plantedsr[ps_notr])
  
  spcoexdat<-rbind(spcoexdat,
    data.frame(sp=sort(unique(tmp$sp_notr)),
      qB_notr=unname(tapply(tmp$qB_notr, tmp$sp_notr, mean)),
      Rs_notr=unname(tapply(tmp$Rs_notr, tmp$sp_notr, mean)),
      pzero_notr=t(matrix(nrow=5, unlist(tapply(tmp$pzero_est_notr, tmp$sp_notr, bootpobs_fun)))),
      qB_tr=unname(tapply(tmp$qB_tr, tmp$sp_notr, mean)),
      Rs_tr=unname(tapply(tmp$Rs_tr, tmp$sp_notr, mean)),
      pzero_tr=t(matrix(nrow=5, unlist(tapply(tmp$pzero_est_tr, tmp$sp_tr, bootpobs_fun)))),
      pzero_obs=t(matrix(nrow=5, unlist(tapply(!is.na(tmp$pzero_obs_tr), tmp$sp_tr, bootfrac_fun)))),
      plsr=sprichlst[i]))
}
coexdatsave<-coexdat
spcoexdatsave<-spcoexdat












#plot
m<-rbind(c(1,1),
         c(1,1),
         c(2,3),
         c(4,5))
layout(m)
ofs1<-c(0.022, -0.06)
ofs2<-c(0.055, -0.02)

coexdat<-coexdatsave
coexdat[,1:3]<-logit(unlist(coexdat[,1:3]))
mxcx<-logit(0.999)
coexdat[,1:3][is.na(coexdat[,1:3]) | coexdat[,1:3]>mxcx]<-mxcx

yrng<-range(coexdat[,1:3])
collst_tr<-adjustcolor(c("black", "red", "blue"), alpha.f = 0.8)


par(mar=c(4,2,2,2), oma=c(2,3,0,0))
plot(c(0.6, length(sprichlst)+0.4), yrng, xlab="", ylab="", type="n", axes=F, cex.lab=1.2)
axis(1, at=1:length(sprichlst), sprichlst)
ax2lst<-c(0.6, 0.8, 0.9, 0.95, 0.99, 0.999)
ax2lst_ps<-logit(ax2lst)
ax2lst[length(ax2lst)]<-(">0.999")
axis(2, at=ax2lst_ps, ax2lst, las=2)
abline(h=ax2lst_ps, col=adjustcolor("black", 0.2))
box()

mtext("Planted Richness", 1, cex=1.2, line=2.4)
mtext("Coexistence Frequency", 2, cex=1.2, line=1.2, outer=T)

psq<-c(-0.26, 0, 0.26)
dsq<-0.12

for(i in 1:length(sprichlst)) {
  pstmp<-coexdat$plsr==sprichlst[i]
  
  for(j in 1:3) {
    polygon(c(i+psq[j]-dsq, i+psq[j]-dsq,
              i+psq[j]+dsq, i+psq[j]+dsq),
            c(0, coexdat[pstmp,][3,j], coexdat[pstmp,][3,j], 0),
            col=collst_tr[j])
    
  }
  
  segments(i+psq,
           unlist(coexdat[pstmp,][1,1:3]),
           i+psq,
           unlist(coexdat[pstmp,][5,1:3]),
           lwd=1.5, lend=2)
  
  segments(i+psq,
           unlist(coexdat[pstmp,][2,1:3]),
           i+psq,
           unlist(coexdat[pstmp,][4,1:3]),
           lwd=3, lend=2)
  
  #segments(i+psq-dsq,
  #         unlist(coexdat[pstmp,][3,1:3]),
  #         i+psq+dsq,
  #         unlist(coexdat[pstmp,][3,1:3]),
  #         lwd=3, lend=2)
}
put.fig.letter("A.", "topleft", offset=ofs1, cex=1.4)

par(mar=c(2,2,1,1))
spcoexdat<-spcoexdatsave
spnumnames<-c("2 Species Mixtures",
              "4 Species Mixtures",
              "8 Species Mixtures",
              "16 Species Mixtures")

dsq<-c(0, 0.02, -0.02)

for(i in 1:length(sprichlst)) {
  sbs<-which(spcoexdat$plsr==sprichlst[i])
  
  plot(1-pzero_notr.3~I(qB_notr+dsq[1]), data=spcoexdat[sbs,][order(spcoexdat[sbs,]$qB_notr),], type="b", col=collst_tr[2],
       xlab="", ylab="", ylim=c(0,1))
  
  segments(spcoexdat[sbs,]$qB_notr+dsq[1],
           1-spcoexdat[sbs,]$pzero_notr.1,
           spcoexdat[sbs,]$qB_notr+dsq[1],
           1-spcoexdat[sbs,]$pzero_notr.5, lwd=1.5, lend=2, col=collst_tr[2])
  
  segments(spcoexdat[sbs,]$qB_notr+dsq[1],
           1-spcoexdat[sbs,]$pzero_notr.2,
           spcoexdat[sbs,]$qB_notr+dsq[1],
           1-spcoexdat[sbs,]$pzero_notr.4, lwd=3, lend=2, col=collst_tr[2])
  
  lines(1-pzero_tr.3~I(qB_notr+dsq[2]), data=spcoexdat[sbs,][order(spcoexdat[sbs,]$qB_tr),], type="b", col=collst_tr[3])
  
  segments(spcoexdat[sbs,]$qB_notr+dsq[2],
           1-spcoexdat[sbs,]$pzero_tr.1,
           spcoexdat[sbs,]$qB_notr+dsq[2],
           1-spcoexdat[sbs,]$pzero_tr.5, lwd=1.5, lend=2, col=collst_tr[3])
  
  segments(spcoexdat[sbs,]$qB_notr+dsq[2],
           1-spcoexdat[sbs,]$pzero_tr.2,
           spcoexdat[sbs,]$qB_notr+dsq[2],
           1-spcoexdat[sbs,]$pzero_tr.4, lwd=3, lend=2, col=collst_tr[3])
  
  
  
  lines(pzero_obs.3~I(qB_notr+dsq[3]), data=spcoexdat[sbs,][order(spcoexdat[sbs,]$qB_tr),], type="b", col=collst_tr[1])
  
  segments(spcoexdat[sbs,]$qB_notr+dsq[3],
           spcoexdat[sbs,]$pzero_obs.1,
           spcoexdat[sbs,]$qB_notr+dsq[3],
           spcoexdat[sbs,]$pzero_obs.5, lwd=1.5, lend=2, col=collst_tr[1])
  
  segments(spcoexdat[sbs,]$qB_notr+dsq[3],
           spcoexdat[sbs,]$pzero_obs.2,
           spcoexdat[sbs,]$qB_notr+dsq[3],
           spcoexdat[sbs,]$pzero_obs.4, lwd=3, lend=2, col=collst_tr[1])
  
  text((max(spcoexdat[sbs,]$qB_notr)-min(spcoexdat[sbs,]$qB_notr))*0.7+min(spcoexdat[sbs,]$qB_notr),
       0.2,
       spnumnames[i], cex=1.3)
  put.fig.letter(paste(LETTERS[i+1], ".", sep=""), "topleft", offset=ofs2, cex=1.4)
  
}

mtext(expression(paste(italic("qB"), "*"[italic(mono)], ", g m"^-1, sep="")), 1, cex=1.2, outer=T, line=0.9)
