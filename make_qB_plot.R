#Make qB vs. R* plot
yz<-log10(ilogit(z)*10^(y))
plot(-x, yz,
     col=1,
     pch=c(0:2, 5, 15:18)[as.numeric(as.factor(as.character(tradeoffdat$fg[kp])))+as.numeric(tradeoffdat$ine120[kp])*4],
     lwd=1,
     xlab="",
     ylab="",
     cex=c(1,1.5)[as.numeric((tradeoffdat$fg[kp]=="L")&(tradeoffdat$ine120[kp]))+1],
     main="",
     cex.main=1.2,axes=F)
mtext(expression(paste(italic("R"), "*",", mg kg"^-1, sep="")), side = 1, line=2.5, cex=0.8)
mtext(expression(paste(italic("qB"), "*"[italic(mono)], ", g m"^-1, sep="")), side = 2, line=2.2, cex=0.8)

sqx<-c(0.055, 0.1, 0.15, 0.25, 0.4, 0.55)
sqy<-c(0.5, 1, 2, 3, 4, 5)
axis(1, -log10(sqx), sqx)
axis(2, log10(sqy), sqy)
box()
trouttmp<-nondirfit(data.frame(-x, yz))
r2=trouttmp$rsq$rsq_est_adj
r<-"R"
rd<-paste(" =", round(r2,3))
text((max(-x)-min(-x))*0.2+min(-x),
     (max(yz)-min(yz))*0.1+min(yz),
     bquote(.(r[1])^2 ~ .(rd[1])))

abline(a=-trouttmp$pars[1]/trouttmp$pars[3], b=-trouttmp$pars[2]/trouttmp$pars[3], lty=2, lwd=1, col="black")
put.fig.letter("G.", location = "topleft", offset=c(0.08, -0.05), cex=1.5)
