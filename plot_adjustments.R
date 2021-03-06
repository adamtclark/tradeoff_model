replacelower<-(0.02)

m<-rbind(c(1,1,1,1,1,10,2,2,2,2,2,3,3,3,3,3),
          c(1,1,1,1,1,10,2,2,2,2,2,3,3,3,3,3),
          c(1,1,1,1,1,10,2,2,2,2,2,3,3,3,3,3),
          c(1,1,1,1,1,10,2,2,2,2,2,3,3,3,3,3),
          c(1,1,1,1,1,10,2,2,2,2,2,3,3,3,3,3),
          c(4,4,4,4,4,10,5,5,5,5,5,6,6,6,6,6),
          c(4,4,4,4,4,10,5,5,5,5,5,6,6,6,6,6),
          c(4,4,4,4,4,10,5,5,5,5,5,6,6,6,6,6),
          c(4,4,4,4,4,10,5,5,5,5,5,6,6,6,6,6),
          c(4,4,4,4,4,10,5,5,5,5,5,6,6,6,6,6),
          c(7,7,7,7,7,10,8,8,8,8,8,9,9,9,9,9),
          c(7,7,7,7,7,10,8,8,8,8,8,9,9,9,9,9),
          c(7,7,7,7,7,10,8,8,8,8,8,9,9,9,9,9),
          c(7,7,7,7,7,10,8,8,8,8,8,9,9,9,9,9),
          c(7,7,7,7,7,10,8,8,8,8,8,9,9,9,9,9))
layout(m)

par(mar=c(2,2,2,1), oma=c(3,3,0.5,0))
############################################################
#Andropogon and Schizachyrium
############################################################
agmondatall<-with(mondatall,
                  aggregate(cbind(no3=no3, abvbiomass=abvbiomass, totaln=totaln), list(monoculture=monoculture, exp=exp, plot=plot, subplot=subplot),
                            function(x) exp(mean(log(x[is.finite(log(x))])))))

subs<-which(agmondatall$monoculture%in%c("Andropogon gerardi", "Schizachyrium scoparium") & agmondatall$exp=="e055")
agmondatall<-agmondatall[subs,]

plot(c(0, 0.5), c(0, 0.5), type="n", xlab="", ylab="")
put.fig.letter("A.", location = "topleft", offset=c(0.05, -0.05), cex=1.5)

tmp<-agmondatall[agmondatall$monoculture=="Andropogon gerardi",]
points(tmp$totaln, tmp$no3, pch=1, lwd=1)
abline(mand<-lm(no3~totaln, tmp), lwd=1.51)

tmp<-agmondatall[agmondatall$monoculture=="Schizachyrium scoparium",]
points(tmp$totaln, tmp$no3, pch=2, lwd=1)
abline(mss<-lm(no3~totaln, tmp), lwd=1.51, lty=2)

mtext("Total Soil N, %", side = 1, line=2.5)
mtext(expression(paste(italic("R"), "*",", mg kg"^-1, sep="")), side = 2, line=2.5)


#plot fit
datout_old<-datoutlst[[2]]
includesp<-c("Andge", "Schsc")

lims<-c(10, 150)
ul<-10^(log10(lims[1])+(log10(lims[2])-log10(lims[1]))*0.6)
bl1<-10^((log10(lims[2])-log10(lims[1]))*0.125+log10(lims[1]))
bl2<-10^((log10(lims[2])-log10(lims[1]))*0.02+log10(lims[1]))

ag_new<-with(datout_altered_andge[datout_altered_andge$sp%in%includesp & datout_altered_andge$plantedsr>1,],
                   aggregate(cbind(obs, est),
                             list(sp=sp, plantedsr=plantedsr),
                             function(x) cbind(hurdlemodel(x), sd(log(x)[is.finite(log(x))]), length(x[is.finite(log(x))]))))
ag_old<-with(datout_old[datout_old$sp%in%includesp & datout_old$plantedsr>1,],
                   aggregate(cbind(obs, est),
                             list(sp=sp, plantedsr=plantedsr),
                             function(x) cbind(hurdlemodel(x), sd(log(x)[is.finite(log(x))]), length(x[is.finite(log(x))]))))

#get p-values
datout_altered_andge$obs[datout_altered_andge$obs<replacelower]<-replacelower
datout_altered_andge$est[datout_altered_andge$est<replacelower]<-replacelower

datout_altered_andge$res<-abs(log10(datout_altered_andge$obs)-log10(datout_altered_andge$est))
tmp1<-with(datout_altered_andge[datout_altered_andge$sp%in%includesp & datout_altered_andge$plantedsr>1,], aggregate(cbind(res=res), list(plt=plt), mean))

datout_old$obs[datout_old$obs<replacelower]<-replacelower
datout_old$est[datout_old$est<replacelower]<-replacelower

datout_old$res<-abs(log10(datout_old$obs)-log10(datout_old$est))
tmp2<-with(datout_old[datout_old$sp%in%includesp & datout_old$plantedsr>1,], aggregate(cbind(res=res), list(plt=plt), mean))

pv<-wilcox.test(tmp1$res, tmp2$res, paired=T)$p.value

#old
plot(ag_old$est[,1], ag_old$obs[,1], pch=as.numeric(as.factor(ag_old$sp)),
     log="xy", lwd=1.51, xlim=lims, ylim=lims,
     xlab="",
     ylab="")
put.fig.letter("B.", location = "topleft", offset=c(0.05, -0.05), cex=1.5)

segments(ag_old$est[,1], exp(log(ag_old$obs[,1])+ag_old$obs[,2]/sqrt(ag_old$obs[,3])), ag_old$est[,1], exp(log(ag_old$obs[,1])-ag_old$obs[,2]/sqrt(ag_old$obs[,3])), lend=2, lwd=1.51)
segments(exp(log(ag_old$est[,1])+ag_old$est[,2]/sqrt(ag_old$est[,3])), ag_old$obs[,1], exp(log(ag_old$est[,1])-ag_old$est[,2]/sqrt(ag_old$est[,3])), ag_old$obs[,1], lend=2, lwd=1.51)
abline(a=0, b=1, lty=3, lwd=1.51)
points(ag_old$est[,1], ag_old$obs[,1], pch=15+as.numeric(as.factor(ag_old$sp)),
       col=grey.colors(4, start=0, end=1)[as.numeric(as.factor(ag_old$plantedsr))])
text(ul, bl2, paste("MAE =", round(10^mean(tmp2$res)-1, 2)), adj=c(0,0), cex=1.2)

#adj
plot(ag_new$est[,1], ag_new$obs[,1], pch=as.numeric(as.factor(ag_new$sp)),
     log="xy", lwd=1.51, xlim=lims, ylim=lims,
     xlab="",
     ylab="")
put.fig.letter("C.", location = "topleft", offset=c(0.05, -0.05), cex=1.5)

segments(ag_new$est[,1], exp(log(ag_new$obs[,1])+ag_new$obs[,2]/sqrt(ag_new$obs[,3])), ag_new$est[,1], exp(log(ag_new$obs[,1])-ag_new$obs[,2]/sqrt(ag_new$obs[,3])), lend=2, lwd=1.51)
segments(exp(log(ag_new$est[,1])+ag_new$est[,2]/sqrt(ag_new$est[,3])), ag_new$obs[,1], exp(log(ag_new$est[,1])-ag_new$est[,2]/sqrt(ag_new$est[,3])), ag_new$obs[,1], lend=2, lwd=1.51)
abline(a=0, b=1, lty=3, lwd=1.51)
points(ag_new$est[,1], ag_new$obs[,1], pch=15+as.numeric(as.factor(ag_new$sp)),
       col=grey.colors(4, start=0, end=1)[as.numeric(as.factor(ag_new$plantedsr))])
text(ul, bl1, paste("MAE =", round(10^mean(tmp1$res)-1, 2)), adj=c(0,0), cex=1.2)
text(ul, bl2, paste("p <", ceiling(pv*1000)/1000), adj=c(0,0), cex=1.2)

############################################################
#Lupine and Amorpha
############################################################
pltlst<-pltloc$plot[pltloc$Luppe>0 | pltloc$Amoca>0 | pltloc$Lesca>0 | pltloc$Petpu>0]
subs<-which(e120dat$Plot%in%pltlst)
agpolydatall<-with(e120dat[subs,], aggregate(cbind(LP=Planted.Lupinus.perennis.g.m2., AC=Planted.Amorpha.canescens.g.m2.,
                                                   LC=Planted.Lespedeza.capitata.g.m2., PP=Planted.Petalostemum.purpureum.g.m2.),
                                             list(plot=Plot, plantedsr=NumSp),
                                             hurdlemodel))
agpolydatall[is.na(agpolydatall)]<-0
for(i in 2:4) {
  pltlst<-pltloc$plot[!as.logical(pltloc[,c("Luppe", "Amoca", "Lesca", "Petpu")[i]])]
  agpolydatall[agpolydatall$plot%in%pltlst,i+2]<-NA
}

tmpLeg<-agpolydatall$LP
tmpLeg[tmpLeg==0]<-runif(sum(tmpLeg==0), -1, 1)
matplot(tmpLeg, agpolydatall[,c("AC", "LC", "PP")],
        pch=c(0, 13, 5),
        col=1,
        cex=0.8,
        xlab="",
        ylab="",
        xlim=c(0, 85),
        ylim=c(0,150))
put.fig.letter("D.", location = "topleft", offset=c(0.05, -0.05), cex=1.5)

lpseq<-seq(0, max(agpolydatall$LP), length=100)
Lllst<-c(3,2,1)
Lcollst<-c(1,1,adjustcolor(1, alpha.f = 0.75))
for(i in 1:3) {
  moddat<-data.frame(agpolydatall[,c(3, 3+i)])
  colnames(moddat)<-c("LP", "y")
  mod<-loess(y~LP, moddat, span=1.25)
  lines(lpseq, predict(mod, newdata=data.frame(LP=lpseq)), lty=Lllst[i], lwd=1.51, col=Lcollst[i])
}

mtext(expression(paste(italic("B"), "* Lupine,", " g m"^-2, sep="")), side = 1, line=2.5)
mtext(expression(paste(italic("B"), "* Other Legumes,", " g m"^-2, sep="")), side = 2, line=2.5)

#plot fit
datout_old<-datoutlst[[2]]
includesp<-c("Luppe", "Amoca", "Lesca", "Petpu")

lims<-c(5, 100)
ul<-10^(log10(lims[1])+(log10(lims[2])-log10(lims[1]))*0.6)
bl1<-10^((log10(lims[2])-log10(lims[1]))*0.125+log10(lims[1]))
bl2<-10^((log10(lims[2])-log10(lims[1]))*0.02+log10(lims[1]))

ag_new<-with(datout_altered_luppe[datout_altered_luppe$sp%in%includesp & datout_altered_luppe$plantedsr>1,],
             aggregate(cbind(obs, est),
                       list(sp=sp, plantedsr=plantedsr),
                       function(x) cbind(hurdlemodel(x), sd(log(x)[is.finite(log(x))]), length(x[is.finite(log(x))]))))
ag_old<-with(datout_old[datout_old$sp%in%includesp & datout_old$plantedsr>1,],
             aggregate(cbind(obs, est),
                       list(sp=sp, plantedsr=plantedsr),
                       function(x) cbind(hurdlemodel(x), sd(log(x)[is.finite(log(x))]), length(x[is.finite(log(x))]))))

#get p-values
datout_altered_luppe$obs[datout_altered_luppe$obs<replacelower]<-replacelower
datout_altered_luppe$est[datout_altered_luppe$est<replacelower]<-replacelower

datout_altered_luppe$res<-abs(log10(datout_altered_luppe$obs)-log10(datout_altered_luppe$est))
tmp1<-with(datout_altered_luppe[datout_altered_luppe$sp%in%includesp & datout_altered_luppe$plantedsr>1,], aggregate(cbind(res=res), list(plt=plt), mean))

datout_old$obs[datout_old$obs<replacelower]<-replacelower
datout_old$est[datout_old$est<replacelower]<-replacelower

datout_old$res<-abs(log10(datout_old$obs)-log10(datout_old$est))
tmp2<-with(datout_old[datout_old$sp%in%includesp & datout_old$plantedsr>1,], aggregate(cbind(res=res), list(plt=plt), mean))

pv<-wilcox.test(tmp1$res, tmp2$res, paired=T)$p.value

#old
plot(ag_old$est[,1], ag_old$obs[,1], pch=c(0, 13, 6, 5)[as.numeric(as.factor(ag_old$sp))],
     log="xy", lwd=1.51, xlim=lims, ylim=lims,
     xlab="",
     ylab="")
put.fig.letter("E.", location = "topleft", offset=c(0.05, -0.05), cex=1.5)

segments(ag_old$est[,1], exp(log(ag_old$obs[,1])+ag_old$obs[,2]/sqrt(ag_old$obs[,3])), ag_old$est[,1], exp(log(ag_old$obs[,1])-ag_old$obs[,2]/sqrt(ag_old$obs[,3])), lend=2, lwd=1.51)
segments(exp(log(ag_old$est[,1])+ag_old$est[,2]/sqrt(ag_old$est[,3])), ag_old$obs[,1], exp(log(ag_old$est[,1])-ag_old$est[,2]/sqrt(ag_old$est[,3])), ag_old$obs[,1], lend=2, lwd=1.51)
abline(a=0, b=1, lty=3, lwd=1.51)

cls<-grey.colors(4, start=0, end=1)[as.numeric(as.factor(ag_old$plantedsr))]
cls[ag_old$sp=="Luppe"]<-1
points(ag_old$est[,1], ag_old$obs[,1], pch=c(15,16,25,18)[as.numeric(as.factor(ag_old$sp))],
       col=cls,
       cex=c(1,1,1,1.2)[as.numeric(as.factor(ag_new$sp))], bg=grey.colors(4, start=0, end=1)[3])
subs<-which(ag_old$sp=="Lesca")
points(ag_old$est[subs,1], ag_old$obs[subs,1], pch=13, lwd=0.8)
text(ul, bl2, paste("MAE =", round(10^mean(tmp2$res)-1, 2)), adj=c(0,0), cex=1.2)


mtext(expression(paste("Observed Abundance, g m"^-2)), side = 2, line=3, cex=1.2)
#adj
plot(ag_new$est[,1], ag_new$obs[,1], pch=c(0, 13, 6, 5)[as.numeric(as.factor(ag_new$sp))],
     log="xy", lwd=1.51, xlim=lims, ylim=lims,
     xlab="",
     ylab="")
put.fig.letter("F.", location = "topleft", offset=c(0.05, -0.05), cex=1.5)

segments(ag_new$est[,1], exp(log(ag_new$obs[,1])+ag_new$obs[,2]/sqrt(ag_new$obs[,3])), ag_new$est[,1], exp(log(ag_new$obs[,1])-ag_new$obs[,2]/sqrt(ag_new$obs[,3])), lend=2, lwd=1.51)
segments(exp(log(ag_new$est[,1])+ag_new$est[,2]/sqrt(ag_new$est[,3])), ag_new$obs[,1], exp(log(ag_new$est[,1])-ag_new$est[,2]/sqrt(ag_new$est[,3])), ag_new$obs[,1], lend=2, lwd=1.51)
abline(a=0, b=1, lty=3, lwd=1.51)

cls<-grey.colors(4, start=0, end=1)[as.numeric(as.factor(ag_old$plantedsr))]
cls[ag_new$sp=="Luppe"]<-1
points(ag_new$est[,1], ag_new$obs[,1], pch=c(15,16,25,18)[as.numeric(as.factor(ag_new$sp))],
       col=cls,
       cex=c(1,1,1,1.2)[as.numeric(as.factor(ag_new$sp))], bg=grey.colors(4, start=0, end=1)[3])
subs<-which(ag_new$sp=="Lesca")
points(ag_new$est[subs,1], ag_new$obs[subs,1], pch=13, lwd=0.8)
text(ul, bl1, paste("MAE =", round(10^mean(tmp1$res)-1, 2)), adj=c(0,0), cex=1.2)
text(ul, bl2, paste("p <", ceiling(pv*1000)/1000), adj=c(0,0), cex=1.2)


############################################################
#Early season
############################################################
agmondatall<-with(dattot,
                  aggregate(cbind(no3=no3, abvbiomass=abvbiomass, totaln=totaln), list(monoculture=monoculture, exp=exp, plot=plot, subplot=subplot, month=month),
                            function(x) exp(mean(log(x[is.finite(log(x))])))))

subs<-which(agmondatall$monoculture%in%c("Poa pratensis", "Liatris aspera"))
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

#plot fit
datout_old<-datoutlst[[2]] #no within-species variation
includesp<-c("Poapr", "Liaas")

lims<-c(0.1, 150)
ul<-10^(log10(lims[1])+(log10(lims[2])-log10(lims[1]))*0.6)
bl1<-10^((log10(lims[2])-log10(lims[1]))*0.125+log10(lims[1]))
bl2<-10^((log10(lims[2])-log10(lims[1]))*0.02+log10(lims[1]))

#get p-values
datout_altered_early$obs[datout_altered_early$obs<replacelower | !is.finite(datout_altered_early$obs)]<-replacelower
datout_altered_early$est[datout_altered_early$est<replacelower | !is.finite(datout_altered_early$est)]<-replacelower

datout_altered_early$res<-abs(log10(datout_altered_early$obs)-log10(datout_altered_early$est))
tmp1<-with(datout_altered_early[datout_altered_early$sp%in%includesp & datout_altered_early$plantedsr>1,], aggregate(cbind(res=res), list(plt=plt), mean))

datout_old$obs[datout_old$obs<replacelower | !is.finite(datout_old$obs)]<-replacelower
datout_old$est[datout_old$est<replacelower | !is.finite(datout_old$est)]<-replacelower

datout_old$res<-abs(log10(datout_old$obs)-log10(datout_old$est))
tmp2<-with(datout_old[datout_old$sp%in%includesp & datout_old$plantedsr>1,], aggregate(cbind(res=res), list(plt=plt), mean))

pv<-wilcox.test(tmp1$res, tmp2$res, paired=T)$p.value

ag_new<-with(datout_altered_early[datout_altered_early$sp%in%includesp & datout_altered_early$plantedsr>1,],
             aggregate(cbind(obs, est),
                       list(sp=sp, plantedsr=plantedsr),
                       function(x) cbind(hurdlemodel(x), sd(log(x)[is.finite(log(x))]), length(x[is.finite(log(x))]))))
ag_old<-with(datout_old[datout_old$sp%in%includesp & datout_old$plantedsr>1,],
             aggregate(cbind(obs, est),
                       list(sp=sp, plantedsr=plantedsr),
                       function(x) cbind(hurdlemodel(x), sd(log(x)[is.finite(log(x))]), length(x[is.finite(log(x))]))))

#old
plot(ag_old$est[,1], ag_old$obs[,1], pch=c(18,15,15)[as.numeric(as.factor(ag_old$sp))],
     log="xy", xlim=lims, ylim=lims,
     xlab="",
     ylab="",
     type="n")

segments(ag_old$est[,1], exp(log(ag_old$obs[,1])+ag_old$obs[,2]/sqrt(ag_old$obs[,3])), ag_old$est[,1], exp(log(ag_old$obs[,1])-ag_old$obs[,2]/sqrt(ag_old$obs[,3])), lend=2, lwd=1.51)
segments(exp(log(ag_old$est[,1])+ag_old$est[,2]/sqrt(ag_old$est[,3])), ag_old$obs[,1], exp(log(ag_old$est[,1])-ag_old$est[,2]/sqrt(ag_old$est[,3])), ag_old$obs[,1], lend=2, lwd=1.51)
points(ag_old$est[,1], ag_old$obs[,1], pch=c(18,15,15)[as.numeric(as.factor(ag_old$sp))],
       lwd=1,
       col=grey.colors(4, start=0, end=1)[as.numeric(as.factor(ag_old$plantedsr))],
       cex=c(1.8,1.4,1.4)[as.numeric(as.factor(ag_old$sp))])

abline(a=0, b=1, lty=3, lwd=1.51)
points(ag_old$est[,1], ag_old$obs[,1], pch=c(9,12,14)[as.numeric(as.factor(ag_old$sp))],
       col=c("darkgrey", 1,1,1)[as.numeric(as.factor(ag_old$plantedsr))], lwd=1,
       cex=1.4)
text(ul, bl2, paste("MAE =", round(10^mean(tmp2$res)-1, 2)), adj=c(0,0), cex=1.2)

#adj
plot(ag_new$est[,1], ag_new$obs[,1], pch=c(18,15,15)[as.numeric(as.factor(ag_new$sp))],
     log="xy", xlim=lims, ylim=lims,
     xlab="",
     ylab="",
     type="n")

segments(ag_new$est[,1], exp(log(ag_new$obs[,1])+ag_new$obs[,2]/sqrt(ag_new$obs[,3])), ag_new$est[,1], exp(log(ag_new$obs[,1])-ag_new$obs[,2]/sqrt(ag_new$obs[,3])), lend=2, lwd=1.51)
segments(exp(log(ag_new$est[,1])+ag_new$est[,2]/sqrt(ag_new$est[,3])), ag_new$obs[,1], exp(log(ag_new$est[,1])-ag_new$est[,2]/sqrt(ag_new$est[,3])), ag_new$obs[,1], lend=2, lwd=1.51)
points(ag_new$est[,1], ag_new$obs[,1], pch=c(18,15,15)[as.numeric(as.factor(ag_new$sp))],
       lwd=1,
       col=grey.colors(4, start=0, end=1)[as.numeric(as.factor(ag_new$plantedsr))],
       cex=c(1.8,1.4,1.4)[as.numeric(as.factor(ag_new$sp))])

abline(a=0, b=1, lty=3, lwd=1.51)
points(ag_new$est[,1], ag_new$obs[,1], pch=c(9,12,14)[as.numeric(as.factor(ag_new$sp))],
       col=c("darkgrey", 1,1,1)[as.numeric(as.factor(ag_new$plantedsr))], lwd=1,
       cex=1.4)
text(ul, bl1, paste("MAE =", round(10^mean(tmp1$res)-1, 2)), adj=c(0,0), cex=1.2)
text(ul, bl2, paste("p <", ceiling(pv*1000)/1000), adj=c(0,0), cex=1.2)




############################################################
#Labels
############################################################
mtext(expression(paste("Predicted Abundance, g m"^-2)), side = 1, outer = TRUE, adj = 0.78, line=2, cex=1.2)
mtext("Original Model", side = 3, outer = TRUE, adj = 0.46, line=-3.4, cex=0.9)
mtext("Augmented Model", side = 3, outer = TRUE, adj = 0.833, line=-3.4, cex=0.9)

mtext("Original Model", side = 3, outer = TRUE, adj = 0.46, line=-22.6, cex=0.9)
mtext("Augmented Model", side = 3, outer = TRUE, adj = 0.833, line=-22.6, cex=0.9)


mtext("Original Model", side = 3, outer = TRUE, adj = 0.46, line=-41.6, cex=0.9)
mtext("Augmented Model", side = 3, outer = TRUE, adj = 0.833, line=-41.6, cex=0.9)
