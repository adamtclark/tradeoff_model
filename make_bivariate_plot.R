#################################################
# Bi-variate plots
#################################################
####RAW
plot(-x, y,
     col=1,
     pch=c(0:2, 5, 15:18)[as.numeric(as.factor(as.character(tradeoffdat$fg[kp])))+as.numeric(tradeoffdat$ine120[kp])*4],
     lwd=1,
     xlab="",
     ylab="",
     cex=c(1,1.5)[as.numeric((tradeoffdat$fg[kp]=="L")&(tradeoffdat$ine120[kp]))+1],
     main="",
     cex.main=1.2,axes=F)
mtext(expression(paste(italic("R"), "*",", mg kg"^-1, sep="")), side = 1, line=2.5, cex=0.8)
mtext(expression(paste(italic("B"), "*"[italic(mono)], ", g m"^-1, sep="")), side = 2, line=2.2, cex=0.8)

sqx<-c(0.055, 0.1, 0.15, 0.25, 0.4, 0.55)
sqy<-c(25, 40, 60, 100, 150, 250)
axis(1, -log10(sqx), sqx)
axis(2, log10(sqy), sqy)
box()
trouttmp<-nondirfit(data.frame(-x, y))
r2=trouttmp$rsq$rsq_est_adj
r<-"R"
rd<-paste(" =", round(r2,3))
text((max(-x)-min(-x))*0.2+min(-x),
     (max(y)-min(y))*0.1+min(y),
     bquote(.(r[1])^2 ~ .(rd[1])))

abline(a=-trouttmp$pars[1]/trouttmp$pars[3], b=-trouttmp$pars[2]/trouttmp$pars[3], lty=2, lwd=1, col="black")
put.fig.letter("A.", location = "topleft", offset=c(0.2, -0.05), cex=1.5)

plot(-x, z,
     col=1,
     pch=c(0:2, 5, 15:18)[as.numeric(as.factor(as.character(tradeoffdat$fg[kp])))+as.numeric(tradeoffdat$ine120[kp])*4],
     lwd=1,
     xlab="",
     ylab="",
     cex=c(1,1.5)[as.numeric((tradeoffdat$fg[kp]=="L")&(tradeoffdat$ine120[kp]))+1],
     main="",
     cex.main=1.2,axes=F)
mtext(expression(paste(italic("R"), "*",", mg kg"^-1, sep="")), side = 1, line=2.5, cex=0.8)
mtext(expression(paste(italic(q), ", g g"^-1, sep="")), side = 2, line=2.2, cex=0.8)

sqz<-c(0.008, 0.01, 0.012, 0.015, 0.02)
axis(1, -log10(sqx), sqx)
axis(2, logit(sqz), sqz)
box()
trouttmp<-nondirfit(data.frame(-x, z))
r2=trouttmp$rsq$rsq_est_adj
r<-"R"
rd<-paste(" =", round(r2,3))
text((max(-x)-min(-x))*0.2+min(-x),
     (max(z)-min(z))*0.1+min(z),
     bquote(.(r[1])^2 ~ .(rd[1])))
abline(a=-trouttmp$pars[1]/trouttmp$pars[2], b=-trouttmp$pars[3]/trouttmp$pars[2], lty=2, lwd=1, col="black")
put.fig.letter("B.", location = "topleft", offset=c(0.2, -0.05), cex=1.5)

plot(z, y,
     col=1,
     pch=c(0:2, 5, 15:18)[as.numeric(as.factor(as.character(tradeoffdat$fg[kp])))+as.numeric(tradeoffdat$ine120[kp])*4],
     lwd=1,
     xlab="",
     ylab="",
     cex=c(1,1.5)[as.numeric((tradeoffdat$fg[kp]=="L")&(tradeoffdat$ine120[kp]))+1],
     main="",
     cex.main=1.2,axes=F)
mtext(expression(paste(italic(q), ", g g"^-1, sep="")), side = 1, line=2.5, cex=0.8)
mtext(expression(paste(italic("B"), "*"[italic(mono)], ", g m"^-1, sep="")), side = 2, line=2.2, cex=0.8)
axis(1, logit(sqz), sqz)
axis(2, log10(sqy), sqy)
box()
trouttmp<-nondirfit(data.frame(z, y))
r2=trouttmp$rsq$rsq_est_adj
r<-"R"
rd<-paste(" =", round(r2,3))
text((max(z)-min(z))*0.2+min(z),
     (max(y)-min(y))*0.1+min(y),
     bquote(.(r[1])^2 ~ .(rd[1])))
abline(a=-trouttmp$pars[1]/trouttmp$pars[2], b=-trouttmp$pars[3]/trouttmp$pars[2], lty=2, lwd=1, col="black")
put.fig.letter("C.", location = "topleft", offset=c(0.2, -0.05), cex=1.5)



########FITTED
xfit<-(-xfit)

plot(xfit, yfit,
     col=1,
     pch=c(0:2, 5, 15:18)[as.numeric(as.factor(as.character(tradeoffdat$fg[kp])))+as.numeric(tradeoffdat$ine120[kp])*4],
     lwd=1,
     xlab="",
     ylab="",
     cex=c(1,1.5)[as.numeric((tradeoffdat$fg[kp]=="L")&(tradeoffdat$ine120[kp]))+1],
     main="",
     cex.main=1.2,axes=F)
mtext(expression(paste(widehat(italic("R")), "*",", mg kg"^-1, sep="")), side = 1, line=2.5, cex=0.8)
mtext(expression(paste(widehat(italic("B")), "*"[italic(mono)], ", g m"^-1, sep="")), side = 2, line=2.2, cex=0.8)

sqx<-c(0.055, 0.1, 0.15, 0.25, 0.4, 0.55)
sqy<-c(25, 40, 60, 100, 150, 250)
axis(1, -log10(sqx), sqx)
axis(2, log10(sqy), sqy)
box()
trouttmp<-nondirfit(data.frame(xfit, yfit))
r2=trouttmp$rsq$rsq_est_adj
r<-"R"
rd<-paste(" =", round(r2,3))
text((max(xfit)-min(xfit))*0.2+min(xfit),
     (max(yfit)-min(yfit))*0.1+min(yfit),
     bquote(.(r[1])^2 ~ .(rd[1])))
abline(a=-trouttmp$pars[1]/trouttmp$pars[3], b=-trouttmp$pars[2]/trouttmp$pars[3], lty=2, lwd=1, col="black")
put.fig.letter("D.", location = "topleft", offset=c(0.2, -0.05), cex=1.5)

plot(xfit, zfit,
     col=1,
     pch=c(0:2, 5, 15:18)[as.numeric(as.factor(as.character(tradeoffdat$fg[kp])))+as.numeric(tradeoffdat$ine120[kp])*4],
     lwd=1,
     xlab="",
     ylab="",
     cex=c(1,1.5)[as.numeric((tradeoffdat$fg[kp]=="L")&(tradeoffdat$ine120[kp]))+1],
     main="",
     cex.main=1.2,axes=F)
mtext(expression(paste(widehat(italic("R")), "*",", mg kg"^-1, sep="")), side = 1, line=2.5, cex=0.8)
mtext(expression(paste(widehat(italic(q)), ", g g"^-1, sep="")), side = 2, line=2.2, cex=0.8)

sqz<-c(0.008, 0.01, 0.012, 0.015, 0.02)
axis(1, -log10(sqx), sqx)
axis(2, logit(sqz), sqz)
box()
trouttmp<-nondirfit(data.frame(xfit, zfit))
r2=trouttmp$rsq$rsq_est_adj
r<-"R"
rd<-paste(" =", round(r2,3))
text((max(xfit)-min(xfit))*0.2+min(xfit),
     (max(zfit)-min(zfit))*0.1+min(zfit),
     bquote(.(r[1])^2 ~ .(rd[1])))
abline(a=-trouttmp$pars[1]/trouttmp$pars[2], b=-trouttmp$pars[3]/trouttmp$pars[2], lty=2, lwd=1, col="black")
put.fig.letter("E.", location = "topleft", offset=c(0.2, -0.05), cex=1.5)

plot(zfit, yfit,
     col=1,
     pch=c(0:2, 5, 15:18)[as.numeric(as.factor(as.character(tradeoffdat$fg[kp])))+as.numeric(tradeoffdat$ine120[kp])*4],
     lwd=1,
     xlab="",
     ylab="",
     cex=c(1,1.5)[as.numeric((tradeoffdat$fg[kp]=="L")&(tradeoffdat$ine120[kp]))+1],
     main="",
     cex.main=1.2,axes=F)
mtext(expression(paste(widehat(italic(q)), ", g g"^-1, sep="")), side = 1, line=2.5, cex=0.8)
mtext(expression(paste(widehat(italic("B")), "*"[italic(mono)], ", g m"^-1, sep="")), side = 2, line=2.2, cex=0.8)
axis(1, logit(sqz), sqz)
axis(2, log10(sqy), sqy)
box()
trouttmp<-nondirfit(data.frame(zfit, yfit))
r2=trouttmp$rsq$rsq_est_adj
r<-"R"
rd<-paste(" =", round(r2,3))
text((max(zfit)-min(zfit))*0.2+min(zfit),
     (max(yfit)-min(yfit))*0.1+min(yfit),
     bquote(.(r[1])^2 ~ .(rd[1])))
abline(a=-trouttmp$pars[1]/trouttmp$pars[2], b=-trouttmp$pars[3]/trouttmp$pars[2], lty=2, lwd=1, col="black")
put.fig.letter("F.", location = "topleft", offset=c(0.2, -0.05), cex=1.5)



xfit<-(-xfit)