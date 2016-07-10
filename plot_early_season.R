############################################################
#Early season species example
############################################################
agmondatall<-with(dattot,
                  aggregate(cbind(no3=no3, abvbiomass=abvbiomass, totaln=totaln), list(monoculture=monoculture, exp=exp, plot=plot, subplot=subplot, month=month),
                            function(x) exp(mean(log(x[is.finite(log(x))])))))

subs<-which(agmondatall$monoculture%in%c("Poa pratensis", "Achillea millefolium(lanulosa)", "Andropogon gerardi"))
agmondatall<-agmondatall[subs,]

plot(c(5.5, 8.5), c(15, 170), type="n", xlab="", ylab="", axes=F, log="y")

axis(2); axis(1, 6:8, c("June", "July", "August")); box()

subslst<-as.character(sort(unique(agmondatall$monoculture)))
adj<-seq(-0.1, 0.1, length=length(subslst))

for(i in 1:length(subslst)) {
  subs<-which(agmondatall$monoculture==subslst[i])
  mupp<-tapply(agmondatall$abvbiomass[subs], agmondatall$month[subs], function(x) exp(mean(log(x[is.finite(log(x))]))))
  sdpp<-tapply(agmondatall$abvbiomass[subs], agmondatall$month[subs], function(x) (sd(log(x[is.finite(log(x))]))/sqrt(length(x[is.finite(log(x))]))))
  
  segments(6:8+adj[i], exp(log(mupp)-sdpp*2), 6:8+adj[i], exp(log(mupp)+sdpp*2), lwd=1.51, lend=2)
  segments(6:8+adj[i], exp(log(mupp)-sdpp), 6:8+adj[i], exp(log(mupp)+sdpp), lwd=4, lend=2)
  
  lines(6:8+adj[i], mupp, lwd=1.51)
  points(6:8+adj[i], mupp, pch=c(18, 15, 15)[i], cex=c(1.8, 1.5, 1.5)[i], col="white")
  points(6:8+adj[i], mupp, pch=c(9,12,14)[i], cex=1.5, lwd=1.5)
}

mtext("Month", side = 1, line=2.5)
mtext(expression(paste(italic("B"), "*"[mono],", g m"^-2, sep="")), side = 2, line=2.5)
