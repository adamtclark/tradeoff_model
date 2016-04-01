#Load filtered tradeoff data
tradeoffdat<-read.csv("data/data_products/filtered_tradeoff_data.csv")

#names of E120 species
splst<-c("Achmi", "Amoca", "Andge", "Asctu", "Koecr", "Lesca", "Liaas",
         "Luppe", "Panvi", "Petpu", "Poapr", "Schsc", "Solri", "Sornu")

#################################################
# Plot tradeoff surface 3D
#################################################
#set up plot layout
m<-cbind(c(rep(1,3), rep(5, 3)),
         c(rep(1,3), rep(5, 3)),
         c(rep(1,3), rep(5, 3)),
         c(rep(1,3), rep(5, 3)),
         c(rep(1,3), rep(5, 3)),
         c(rep(2, 2), rep(3, 2), rep(4,2)),
         c(rep(2, 2), rep(3, 2), rep(4,2)),
         c(rep(2, 2), rep(3, 2), rep(4,2)))

par(oma=c(0.5,0.5,0.5,0.5))
layout(m)

#fit tradeoff surface
tdat_fit<-data.frame(no3=log10(tradeoffdat$no3), abv=log10(tradeoffdat$abv), ptisn=logit(tradeoffdat$ptisn))
tdat_fit$no3<-(-tdat_fit$no3) #reverse no3, becase it is inversely related to competitive hierarchy
trout<-nondirfit(tdat_fit)

x<-tdat_fit$no3
y<-tdat_fit$abv
z<-tdat_fit$ptisn

#extract fitted surface
xfit<-trout$possnap[,"no3"]
yfit<-trout$possnap[,"abv"]
zfit<-trout$possnap[,"ptisn"]
kp<-which(is.finite(x)&is.finite(y)&is.finite(z))
x<-x[kp]; y<-y[kp]; z<-z[kp]
xfit<-xfit[kp]; yfit<-yfit[kp]; zfit<-zfit[kp]

#linear regression for plotting fitted surface (note this will be a perfect fit)
fit <- lm(zfit ~ xfit + yfit)

#predict values on regular xy grid
ngrid<-20

#ranges for plotting
xrng<-c(0.2, 1.3)
yrng<-c(1.2, 2.5)
#zrng<-c(-2.5, -1.24)
zrng<-c(-5.7, -2.8)

x.pred <- seq(xrng[1], xrng[2], length.out = ngrid)
y.pred <- seq(yrng[1], yrng[2], length.out = ngrid)
xy <- expand.grid(xfit = x.pred, 
                  yfit = y.pred)

z.pred <- matrix (nrow = ngrid, ncol = ngrid, 
                  data = predict(fit, newdata = data.frame(xy)))

#fitted points for droplines to surface
fitpoints <- predict(fit)

par(mar=c(0,0,0,0))
scatter3D(z = z, x = x, y = y, col=NULL, colvar=NULL,
          xlab="R*", ylab="B*", zlab="q",
          ylim=c(yrng[1], yrng[2]), xlim=c(xrng[1], xrng[2]), zlim=c(zrng[1], zrng[2]),
          cex=c(1,1.5)[as.numeric((tradeoffdat$fg[kp]=="L")&(tradeoffdat$ine120[kp]))+1],
          pch=c(0:2, 5, 15:18)[as.numeric(as.factor(as.character(tradeoffdat$fg[kp])))+as.numeric(tradeoffdat$ine120[kp])*4],
          type="p",lwd=2, cex.axis = 1e-9,
          theta = 135, phi = 20, ticktype = "detailed",
          surf = list(x = x.pred, y = y.pred, z = z.pred,
                      facets = NA, fit = fitpoints, col=adjustcolor("black", alpha.f = 0.5), lty=0, lwd=1))

text3D(x = c(seq(xrng[1], xrng[2], length=6), rep(xrng[2], 3), rep(xrng[2], 6)), #note axis labels true no3 values
       y = c(rep(yrng[2], 6), seq(yrng[1], yrng[2], length=3), rep(yrng[1], 6)),
       z = c(rep(zrng[1], 6), rep(zrng[1], 3), seq(zrng[1], zrng[2], length=6)),
       labels = c(round(10^(-seq(xrng[1], xrng[2], length=6)), 2), round(10^seq(yrng[1], yrng[2], length=3), 0), round(ilogit(seq(zrng[1], zrng[2], length=6)),3)),
       add = TRUE, adj = 0)

scatter3D(z = z, x = x, y = y, col=NULL, colvar=NULL,
          xlab="R*", ylab="B*", zlab="q",
          ylim=c(yrng[1], yrng[2]), xlim=c(xrng[1], xrng[2]), zlim=c(zrng[1], zrng[2]),
          cex=c(1,1.5)[as.numeric((tradeoffdat$fg[kp]=="L")&(tradeoffdat$ine120[kp]))+1],
          pch=c(0:2, 5, 15:18)[as.numeric(as.factor(as.character(tradeoffdat$fg[kp])))+as.numeric(tradeoffdat$ine120[kp])*4],
          type="p",lwd=2, cex.axis = 1e-9, add=TRUE,
          theta = 135, phi = 20, ticktype = "detailed",
          surf = list(x = x.pred, y = y.pred, z = z.pred,
                      facets = NA, fit = fitpoints, col=adjustcolor("black", alpha.f = 0.5), lty=1, lwd=0.8))

#calculate mean r-squared (coef. of det.) for total surface
obs<-trout$vars[is.finite(rowSums(trout$vars)),]
pred<-trout$possnap[is.finite(rowSums(trout$vars)),]

rss_mod<-colSums((obs-pred)^2)
rss_tot<-colSums(t(colMeans(obs)-t(obs))^2)

r<-"R"
r2<-1-sum(rss_mod)/sum(rss_tot)
rd<-paste(" =", round(r2,3))
text(grconvertX(0.35, "nfc", "user"),
     grconvertY(0.85, "nfc", "user"),
     bquote(.(r[1])^2 ~ .(rd[1])))

put.fig.letter("A.", location = "topleft", offset=c(0.05, 0), cex=1.5)

#################################################
# Make residual plots
#################################################
par(mar=c(3,3,2,1))
############
#R*
############
plot(x, xfit,
     col=1,
     pch=c(0:2, 5, 15:18)[as.numeric(as.factor(as.character(tradeoffdat$fg[kp])))+as.numeric(tradeoffdat$ine120[kp])*4],
     lwd=1,
     xlab="",
     ylab="",
     cex=c(1,1.5)[as.numeric((tradeoffdat$fg[kp]=="L")&(tradeoffdat$ine120[kp]))+1],
     main="",
     cex.main=1.2,axes=F)
sq<-c(0.055, 0.1, 0.15, 0.25, 0.4, 0.55)
axis(1, -log10(sq), sq) #label axes with true no3 values
axis(2, -log10(sq), sq)
box()
abline(a=0, b=1, lty=3, lwd=1.5, col="darkgrey")

#estimate sd around trait values
qt<-rbind(xfit-tradeoffdat$sd_no3[kp]/log(10),
          xfit+tradeoffdat$sd_no3[kp]/log(10))
xv<-x[order(x)]
xvsq<-seq(min(xv), max(xv), length=100)
yv1<-(qt[1,][is.finite(qt[1,])])[order(x)]
yv2<-(qt[2,][is.finite(qt[2,])])[order(x)]
yp1<-predict(loess(yv1~xv, span = 0.5), newdata=data.frame(xv=xvsq))
yp2<-predict(loess(yv2~xv, span = 0.5), newdata=data.frame(xv=xvsq))
lines(xvsq, yp1, lty=2, col=1)
lines(xvsq, yp2, lty=2, col=1)
#label axes
mtext(expression(paste(italic("R"), "*",", mg kg"^-1, sep="")), side = 1, line=2.5, cex=0.8)
mtext(expression(paste(widehat(italic("R")), "*",", mg kg"^-1, sep="")), side = 2, line=2.2, cex=0.8)

#calculate r-squred for fit
r2=1-sum((x-xfit)^2)/sum((x-mean(x))^2)
r<-"R"
rd<-paste(" =", round(r2,3))
text((max(x)-min(x))*0.8+min(x),
     (max(xfit)-min(xfit))*0.1+min(xfit),
     bquote(.(r[1])^2 ~ .(rd[1])))
put.fig.letter("B.", location = "topleft", offset=c(0.2, -0.05), cex=1.5)

############
#B*
############
plot(y, yfit,
     pch=c(0:2, 5, 15:18)[as.numeric(as.factor(as.character(tradeoffdat$fg[kp])))+as.numeric(tradeoffdat$ine120[kp])*4],
     col=1,
     lwd=1,
     xlab="",
     ylab="",
     cex=c(1,1.5)[as.numeric((tradeoffdat$fg[kp]=="L")&(tradeoffdat$ine120[kp]))+1],
     main="",
     cex.main=1.2,
     axes=F)
sq<-c(25, 40, 60, 100, 150, 250)
axis(1, log10(sq), sq)
axis(2, log10(sq), sq)
box()

abline(a=0, b=1, lty=3, lwd=1.5, col="darkgrey")

#estimate sd around trait values
qt<-rbind(yfit-tradeoffdat$sd_abv[kp]/log(10),
          yfit+tradeoffdat$sd_abv[kp]/log(10))
xv<-y[order(y)]
xvsq<-seq(min(xv), max(xv), length=100)
yv1<-(qt[1,][is.finite(qt[1,])])[order(y)]
yv2<-(qt[2,][is.finite(qt[2,])])[order(y)]
yp1<-predict(loess(yv1~xv, span = 0.5), newdata=data.frame(xv=xvsq))
yp2<-predict(loess(yv2~xv, span = 0.5), newdata=data.frame(xv=xvsq))
lines(xvsq, yp1, lty=2, col=1)
lines(xvsq, yp2, lty=2, col=1)

#label axes
mtext(expression(paste(italic("B"), "*"[italic(mono)], ", g m"^-1, sep="")), side = 1, line=2.5, cex=0.8)
mtext(expression(paste(widehat(italic("B")), "*"[italic(mono)], ", g m"^-1, sep="")), side = 2, line=2.2, cex=0.8)

#calculate r-squred for fit
r2=1-sum((y-yfit)^2)/sum((y-mean(y))^2)
rd<-paste(" =", round(r2,3))
text((max(y)-min(y))*0.8+min(y),
     (max(yfit)-min(yfit))*0.1+min(yfit),
     bquote(.(r[1])^2 ~ .(rd[1])))
put.fig.letter("C.", location = "topleft", offset=c(0.2, -0.05), cex=1.5)

############
#q
############
plot(z, zfit,
     pch=c(0:2, 5, 15:18)[as.numeric(as.factor(as.character(tradeoffdat$fg[kp])))+as.numeric(tradeoffdat$ine120[kp])*4],
     col=1,
     lwd=1,
     xlab="",
     ylab="",
     cex=c(1,1.5)[as.numeric((tradeoffdat$fg[kp]=="L")&(tradeoffdat$ine120[kp]))+1],
     main="",
     cex.main=1.2,
     axes=F)
sq<-c(0.008, 0.01, 0.012, 0.015, 0.02)
axis(1, logit(sq), sq)
axis(2, logit(sq), sq)
box()
abline(a=0, b=1, lty=3, lwd=1.5, col="darkgrey")

#estimate sd around trait values
qt<-rbind(zfit-tradeoffdat$sd_ptisn[kp],
          zfit+tradeoffdat$sd_ptisn[kp])
xv<-z[order(z)]
xvsq<-seq(min(xv), max(xv), length=100)
yv1<-(qt[1,][is.finite(qt[1,])])[order(z)]
yv2<-(qt[2,][is.finite(qt[2,])])[order(z)]
yp1<-predict(loess(yv1~xv, span = 0.5), newdata=data.frame(xv=xvsq))
yp2<-predict(loess(yv2~xv, span = 0.5), newdata=data.frame(xv=xvsq))
lines(xvsq, yp1, lty=2, col=1)
lines(xvsq, yp2, lty=2, col=1)

#label axes
mtext(expression(paste(italic(q), ", g g"^-1, sep="")), side = 1, line=2.5, cex=0.8)
mtext(expression(paste(widehat(italic(q)), ", g g"^-1, sep="")), side = 2, line=2.2, cex=0.8)

#calculate r-squred for fit
r2=1-sum((z-zfit)^2)/sum((z-mean(z))^2)
rd<-paste(" =", round(r2,3))
text((max(z)-min(z))*0.8+min(z),
     (max(zfit)-min(zfit))*0.1+min(zfit),
     bquote(.(r[1])^2 ~ .(rd[1])))
put.fig.letter("D.", location = "topleft", offset=c(0.2, -0.05), cex=1.5)


#################################################
# Plot qB* vs. R*
#################################################
xfit<-(-xfit) #return no3 back to true values
x<-(-x)

par(mar=c(3,4,2,2))
qBfit<-yfit+log10(ilogit(zfit))
plot(xfit, qBfit,
     pch=c(0:2, 5, 15:18)[as.numeric(as.factor(as.character(tradeoffdat$fg[kp])))+as.numeric(tradeoffdat$ine120[kp])*4],
     col=1,
     lwd=1,
     cex=c(1,1.5)[as.numeric((tradeoffdat$fg[kp]=="L")&(tradeoffdat$ine120[kp]))+1],
     xlab="",
     ylab="",
     type="n",
     axes=F)
sq<-c(0.06, 0.1, 0.15, 0.25, 0.4, 0.5)
axis(1, log10(sq), sq)
sq<-c(0.4, 0.6, 1, 1.6, 2.5, 4)
axis(2, log10(sq), sq)

#show difference between Panvi and Koecr
ps1<-which(tradeoffdat[kp,]$Species=="Panicum virgatum")
polygon(c(xfit[ps1], 0, 0, xfit[ps1]), c((qBfit)[ps1], (qBfit)[ps1], -1, -1),
        col=adjustcolor("black", alpha.f = 0.2), border="darkgrey", density=20, angle=45)

ps2<-which(tradeoffdat[kp,]$Species=="Koeleria cristata")
polygon(c(xfit[ps2], 0, 0, xfit[ps2]), c((qBfit)[ps2], (qBfit)[ps2], -1, -1),
        col=adjustcolor("black", alpha.f = 0.2), border="darkgrey", density=20, angle=135)
#arrows(-0.6, (qBfit)[ps1], -0.6, (qBfit)[ps2], lwd=2, length=0.1, code=3)
#text(-0.6, 0, "Koecr - Panvi", pos=4, cex=1.1)

box()

#add points
points(xfit, qBfit,
  pch=c(0:2, 5, 15:18)[as.numeric(as.factor(as.character(tradeoffdat$fg[kp])))+as.numeric(tradeoffdat$ine120[kp])*4],
  cex=c(1,1.5)[as.numeric((tradeoffdat$fg[kp]=="L")&(tradeoffdat$ine120[kp]))+1])

#label axes
mtext(expression(paste(widehat(italic("R")), "*",", mg kg"^-1, sep="")), 1, line=2.5, cex=0.8)
mtext(expression(paste(widehat(italic(qB)), "*", ", g m"^-2, sep="")), 2, line=2.2, cex=0.8)

#label E120 species
text(xfit[tradeoffdat[kp,]$ine120], qBfit[tradeoffdat[kp,]$ine120], splst, pos=3)

#add error bars to traits
sdrs<-rbind(xfit-tradeoffdat$sd_no3[kp]/log(10),
            xfit+tradeoffdat$sd_no3[kp]/log(10))

#simulate qB*, as there is no trivial closed form solution
sdqB<-matrix(nrow=2, ncol=length(kp))
for(i in 1:length(kp)) {
  sdqB[,i]<-quantile((exp(rnorm(nrep, log(10^yfit)[i], tradeoffdat$sd_abv[kp][i]))*
    (ilogit(rnorm(nrep, zfit[i], tradeoffdat$sd_ptisn[kp][i])))),
    c(pnorm(-1, 0, 1), pnorm(1, 0, 1)))
}
sdqB<-log10(sdqB)

segments(xfit[tradeoffdat$ine120[kp]], sdqB[1,tradeoffdat$ine120[kp]], xfit[tradeoffdat$ine120[kp]], sdqB[2,tradeoffdat$ine120[kp]], lwd=1)
segments(sdrs[1,tradeoffdat$ine120[kp]],
                   (qBfit)[tradeoffdat$ine120[kp]],
                   sdrs[2,tradeoffdat$ine120[kp]],
                   (qBfit)[tradeoffdat$ine120[kp]], lwd=1)

put.fig.letter("E.", location = "topleft", offset=c(0.154, -0.0255), cex=1.5)


#Do slope test
if(dotradeofftest==TRUE) { #set to TRUE to calculate p-values for slopes
  tdat_fit_full<-tdat_fit[is.finite(rowSums(tdat_fit)),]
  pardat<-matrix(nrow=nrep, ncol=4)
  
  #bootstrap data and refit
  bootfun<-function(...) {
    smpps<-sample(1:nrow(tdat_fit_full), rep=T)
    tdat_fit_sim<-tdat_fit_full[smpps,]
    trout_sim<-nondirfit(tdat_fit_sim)
    trout_sim$pars
  }
  
  clusterExport(cl, c("tdat_fit_full", "nondirfit", "optfun"))
  #run parallel program
  clusterout<-parLapply(cl=cl, 1:nrep, fun=bootfun)
  pardat<-t(matrix(nrow=4, unlist(clusterout)))
  
  colnames(pardat)<-c("(Intercept)", "abv", "no3", "ptisn")
  
  #save output
  write.csv(pardat, "data/data_products/tradeoff_boot_pvalues.csv", row.names=F)
  
  parameter_intervals<-apply(pardat, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
  
  #p-value testing for negative relationship between all traits
  pvale<-1-sum((sign(pardat[,1])==sign(pardat[,2]))&(sign(pardat[,3])==sign(pardat[,4])))/nrow(pardat)
  print(paste("3-way tradeoff p-value:", pvale))
} else if(dotradeofftest=="saved") {
  #load saved output
  pardat<-read.csv("data/data_products/tradeoff_boot_pvalues.csv")
  parameter_intervals<-apply(pardat, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
  #p-value testing for negative relationship between all traits
  pvale<-1-sum((sign(pardat[,2])==sign(pardat[,3]))&(sign(pardat[,3])==sign(pardat[,4])))/nrow(pardat)
  print(paste("3-way tradeoff p-value:", pvale))
}



