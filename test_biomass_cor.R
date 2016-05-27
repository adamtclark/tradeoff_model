dat<-read.csv("data/e120_spAbundance_01-14.csv")
soilinfo<-read.csv("data/e120_psoilN_1994.csv")
plotinfo<-pltloc

#planted species only
plantedbm<-dat
plantedbm[,5:18]<-(plotinfo[match(plantedbm$Plot, plotinfo$plot),2:15])*(plantedbm[,5:18])

#Add planted info
plantedbm_csv<-cbind(plantedbm, plotinfo[match(plantedbm$Plot, plotinfo$plot),2:15])

#1. Get column showing monomass for each species
monomass<-apply(plantedbm_csv[rowSums(plantedbm_csv[,19:32])==1,5:18], 2,
                function(x) mean(x[x!=0]))
names(monomass)<-paste(names(monomass), "_mono", sep="")
monomat<-matrix(rep(monomass, each=nrow(plantedbm_csv)), nrow=nrow(plantedbm_csv))
colnames(monomat)<-names(monomass)
plantedbm_csv<-data.frame(plantedbm_csv, monomat)

plantedbm_csv$ismono<-plantedbm_csv$NumSp==1
plantedbm_csv$realized_richness<-rowSums(plantedbm_csv[,19:32])
plantedbm_csv$soilN<-soilinfo$pNSoil940_20[match(plantedbm_csv$Plot, soilinfo$Plot)]
plantedbm_csv$soilC<-soilinfo$pCSoil940_20[match(plantedbm_csv$Plot, soilinfo$Plot)]


#2. Fit models with each species (ordered by R*)
rstardat<-read.csv("data/data_products/tableS1_paramtertable.csv")

rstarlst<-rstardat$no3_snapped[(rstardat$ine120)]
ord<-order(rstarlst)
spnames<-names(plotinfo)[2:15]
spnames_fll<-names(plantedbm)[5:18]

coefmat<-NULL
nsp<-length(rstarlst)

pdf("figures/Competition_regression_AVP.pdf", width=6, height=4) #added variable plots
for(i in 2:length(rstarlst)) {
  plotuse<-plotinfo$plot[(plotinfo[,ord[i]+1]==1)]
  plotsubs<-plantedbm$Plot%in%plotuse
  
  moddat<-plantedbm[plotsubs,ord[(1:i)]+4]
  moddat<-moddat[subs<-moddat[,i]>0,]
  
  plot<-(plantedbm[plotsubs,]$Plot)[subs]
  year<-(plantedbm[plotsubs,]$Year)[subs]
  plantedrich<-(plantedbm[plotsubs,]$NumSp)[subs]
  realized_richness<-(plantedbm_csv[plotsubs,]$realized_richness)[subs]
  soilC<-(plantedbm_csv[plotsubs,]$soilC)[subs]
  
  xvars<-(as.matrix(moddat[,-i]))
  yvar<-(as.matrix(moddat[,i]))
  
  print(spnames[ord[i]])
  
  for(j in 1:ncol(xvars)) {
    if(ncol(xvars)>1) {
      xvarout<-xvars[,j]

      subs<-1:length(xvarout)
      
      xvarout<-xvarout[subs]
      xvars2<-as.matrix(xvars[subs,-j])
      if(is.null(colnames(xvars2)[1])) {
        colnames(xvars2)<-colnames(xvars)[-j]
      }
      
      Plot<-plot[subs]
      gmoddat<-data.frame(y=yvar[subs], xvars2, Plot=Plot)
      
      eq<-paste("lm(log10(y+0.01)~(log10(soilC)+year)^2+", paste(paste("log10(", colnames(xvars2), "+0.01)", sep=""), collapse=" + "), ",",
            "data=gmoddat)")#, corr=corCompSymm(,form=~1|Plot), method='ML')", sep="")
     
      mod2<-eval(parse(text=eq))
      mod2<-step(mod2, trace = 0)
      
      resids<-residuals(mod2)
      
      modx<-update(mod2, log10(xvarout+0.01)~.)
      
      xresids<-residuals(modx)
      
      modfull<-update(mod2, .~.+log10(xvarout+0.01))
    } else {
      #subs<-(xvars>=0)
      subs<-1:length(xvars)
      
      y<-yvar[subs]
      Plot<-plot[subs]
      mod2<-lm(log10(y+0.01)~(log10(soilC)+log10(year))^2)
      resids<-residuals(mod2)
      
      xvarout<-xvars[subs]
      
      modx<-lm(log10(xvarout+0.01)~(log10(soilC)+log10(year))^2)
      
      xresids<-residuals(modx)
      
      modfull<-update(mod2, .~.+log10(xvarout+0.01))
    }
    
    mod3<-lm(resids~xresids)
    
    plot(resids~xresids, xlab=spnames[ord[1:i]][j], ylab="Residuals", main=spnames[ord[i]], cex=0.5)
    
    mux<-mean(xresids, na.rm=T)
    muy<-mean(resids, na.rm=T)
    
    abline(h=muy, lwd=2, col=1, lty=3)
    coefmod<-coef(mod3)
    
    abline(a=coefmod[1], b=coefmod[2], lwd=2, col=2)
    
    sumry<-summary(mod3)$coefficients
    slp<-sumry[2,1]-sumry[2,2]*2
    int<-muy-(mux*slp)
    abline(a=int, b=slp, col=2, lwd=1, lty=2)
    
    slp<-sumry[2,1]+sumry[2,2]*2
    int<-muy-(mux*slp)
    abline(a=int, b=slp, col=2, lwd=1, lty=2)
    
    r2<-1-sum((resids-c(predict(mod3)))^2, na.rm=T)/sum((resids-mean(resids,na.rm=T))^2, na.rm=T)
    
    xstart<-min(xresids)+(max(xresids)-min(xresids))*0.2
    ystart_d<-(max(resids)-min(resids))
    
    text(xstart, min(resids)+ystart_d*0.35, paste("R2-adj =", round(r2,3)), pos=2)
    text(xstart, min(resids)+ystart_d*0.2, paste("p <", ceiling(sumry[2,4]*1000)/1000), pos=2)
    text(xstart, min(resids)+ystart_d*0.05, paste("slp =", round(sumry[2,1],3)), pos=2)
    
    ps<-2
    coefmat<-rbind(coefmat, data.frame(slp=sumry[ps,1], slp_sd=sumry[ps,2], rs_i=sort(rstarlst)[i], rs_j=sort(rstarlst)[j], pvalue=sumry[ps,4], spi=spnames[ord][i], spj=spnames[ord][j]))
  }
}
dev.off()



#3. Plot sign of the effects
pvcutoff<-0.05
fglst<-c("C3", "C4", "F", "L")[c(3,4,2,3,1,4,3,4,2,4,1,2,3,2)]

pdf("figures/FigureS1_Competition_regression_table.pdf", width=5, height=4.5)
par(mar=c(5,5,5,5))
plot(c(0, 1), c(0, 1), type="n", axes=F, xlab="", ylab="", yaxs="i", xaxs="i")
sq<-seq(0, 1, length=(length(ord)))

coefmat$binslp<-NA
magslp<-max(abs(c(min(coefmat[coefmat$pvalue<=pvcutoff,]$slp), max(coefmat[coefmat$pvalue<=pvcutoff,]$slp))))
magslp<-ceiling(magslp*100)/100

cutlst<-cut(coefmat[coefmat$pvalue<=pvcutoff,]$slp, breaks = seq(-magslp, magslp, length=3))
coefmat[coefmat$pvalue<=pvcutoff,]$binslp<-cutlst
cutlst<-levels(cutlst)
collst<-c("black", "white")

for(j in 1:(length(ord)-1)) {
  for(i in j:length(ord)) {
    subs<-which(coefmat$spi==spnames[ord][i] & coefmat$spj==spnames[ord][j])
    if(sum(subs)>0) {
      pv<-coefmat[subs,"pvalue"]
      if(pv<=pvcutoff) {
        eff<-coefmat[subs,"slp"]
        binslp<-coefmat[subs,"binslp"]
        
        polygon(sq[c(j, j, j+1, j+1)], sq[c(i, i+1, i+1, i)-1], col=collst[binslp])
      } else {
        polygon(sq[c(j, j, j+1, j+1)], sq[c(i, i+1, i+1, i)-1], col="grey")
        polygon(sq[c(j, j, j+1, j+1)], sq[c(i, i+1, i+1, i)-1], col="black", density = 10, angle=45)
      }
    }
  }
}
axis(3, sq[-14]+mean(diff(sq))/2, paste(spnames, " (", fglst, ")", sep="")[ord][-14], las=2, lty = 0, cex.axis=0.8)
axis(2, sq[-14]+mean(diff(sq))/2, paste(spnames, " (", fglst, ")", sep="")[ord][-1], las=2, lty = 0, cex.axis=0.8)

axis(1, sq[-14]+mean(diff(sq))/2, round(rstarlst[ord],3)[-14], las=2, cex.axis=0.8)
mtext(expression(paste(widehat(italic("R"[j])), "*, mg kg"^-1)), side = 1, line=4, cex=1.2)

axis(4, sq[-14]+mean(diff(sq))/2, round(rstarlst[ord],3)[-1], las=2, cex.axis=0.8)
mtext(expression(paste(widehat(italic("R"[i])), "*, mg kg"^-1)), side = 4, line=4, cex=1.2)

dsq<-mean(diff(sq))
polygon(sq[rep(11,4)]+dsq/2*c(0,0,1,1),
        sq[rep(4,4)]-dsq/2*c(0,1,1,0), col="black")
polygon(sq[rep(11,4)]+dsq/2*c(0,0,1,1),
        sq[rep(3,4)]-dsq/2*c(0,1,1,0), col="white")
polygon(sq[rep(11,4)]+dsq/2*c(0,0,1,1),
        sq[rep(2,4)]-dsq/2*c(0,1,1,0), col="grey")
polygon(sq[rep(11,4)]+dsq/2*c(0,0,1,1),
        sq[rep(2,4)]-dsq/2*c(0,1,1,0), col="black", density=15)
text(sq[11]+dsq*1.2, sq[4]-dsq/4, "< 0")
text(sq[11]+dsq*1.2, sq[3]-dsq/4, "> 0")
text(sq[11]+dsq*1.2, sq[2]-dsq/4, "n.s.")
text(sq[11]-dsq/2.3, sq[4]+dsq*(3/4),
     expression(paste(beta[paste(italic(j), " on ", italic("i"))]))
     , pos=4, cex=1.2)
dev.off()

