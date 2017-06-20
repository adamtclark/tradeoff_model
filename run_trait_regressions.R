########################################################
# Prep data for regressions
########################################################
nrep_coef_boot<-0
#nrep_coef_boot<-nrep_traits #un-comment this to get p-values for the major axis regressions


#e120 wide to long
e120dat_long<-(melt(data.table(e120dat), id.vars = names(e120dat)[c(1:4,33)],
     measure.vars = names(e120dat)[19:32], variable.name = "Species", value.name = "Biomass"))

#fix plrich
e120dat_long$NumSp<-pltloc$plantedsr[match(e120dat_long$Plot, pltloc$plot)]

#get mean within each plot
e120dat_long<-data.frame(e120dat_long[,j=list(Biomass=mean(Biomass)), by=list(Plot, NumSp, Species)])

#fix species names
e120dat_long$Species<-gsub(" g m2 ", "", gsub("  g m2 ", "", gsub("Planted ", "", gsub(".", " ", as.character(e120dat_long$Species), fixed=T))))
e120dat_long$Species[e120dat_long$Species=="Achillea millefolium lanulosa"]<-"Achillea millefolium(lanulosa)"

#remove non-planted instances
e120dat_long$Planted<-NA
e120dat_long<-e120dat_long[order(e120dat_long$Plot, e120dat_long$Species),]
for(i in 1:nrow(pltloc)) {
  sbs<-which(e120dat_long$Plot==pltloc$plot[i])
  e120dat_long$Planted[sbs]<-unname(unlist(pltloc[i,2:15]))
}
e120dat_long<-e120dat_long[e120dat_long$Planted==1,]

#merge in trait data
e120dat_long<-data.frame(e120dat_long, tradeoffdat[match(e120dat_long$Species, as.character(tradeoffdat$Species)),-1])

#merge in soil C data
e120dat_long$soilC<-e120_tn$pCSoil940_20[match(e120dat_long$Plot, e120_tn$Plot)]

########################################################
#fit models
########################################################
stepfun_use<-function(x, ...) {x}
#stepfun_use<-step #un-comment this line for step-wise regression simplification

############################
#without plrich
############################
m0<-stepfun_use(lm(log(Biomass)~log(soilC), e120dat_long[is.finite(log(e120dat_long$Biomass)),]), trace = FALSE)
summary(m0)

m1.1<-stepfun_use(lm(log(Biomass)~log(soilC)*log(abv), e120dat_long[is.finite(log(e120dat_long$Biomass)),]), trace=F)
m1.2<-stepfun_use(lm(log(Biomass)~log(soilC)*logit(ptisn), e120dat_long[is.finite(log(e120dat_long$Biomass)),]), trace=F)
m1.3<-stepfun_use(lm(log(Biomass)~log(soilC)*log(no3), e120dat_long[is.finite(log(e120dat_long$Biomass)),]), trace=F)
summary(m1.1)
summary(m1.2)
summary(m1.3)

m2<-stepfun_use(lm(log(Biomass)~log(soilC)*(log(abv)+logit(ptisn)+log(no3)), e120dat_long[is.finite(log(e120dat_long$Biomass)),]), trace=F)
summary(m2)

m3<-stepfun_use(lm(log(Biomass)~log(soilC)*(log(abv)+logit(ptisn)+log(no3))^2, e120dat_long[is.finite(log(e120dat_long$Biomass)),]), trace=F)
summary(m3)

m4<-stepfun_use(lm(log(Biomass)~log(soilC)*(log(abv)+logit(ptisn)+log(no3))^3, e120dat_long[is.finite(log(e120dat_long$Biomass)),]), trace=F)
summary(m4)

############################
#with plrich
############################
m0.PR<-stepfun_use(lm(log(Biomass)~log(soilC)*log(NumSp), e120dat_long[is.finite(log(e120dat_long$Biomass)),]), trace=F)
summary(m0.PR)

m1.1.PR<-stepfun_use(lm(log(Biomass)~log(soilC)*log(NumSp)*log(abv), e120dat_long[is.finite(log(e120dat_long$Biomass)),]), trace=F)
m1.2.PR<-stepfun_use(lm(log(Biomass)~log(soilC)*log(NumSp)*logit(ptisn), e120dat_long[is.finite(log(e120dat_long$Biomass)),]), trace=F)
m1.3.PR<-stepfun_use(lm(log(Biomass)~log(soilC)*log(NumSp)*log(no3), e120dat_long[is.finite(log(e120dat_long$Biomass)),]), trace=F)
summary(m1.2.PR)
summary(m1.3.PR)

m2.PR<-stepfun_use(lm(log(Biomass)~log(soilC)*log(NumSp)*(log(abv)+logit(ptisn)+log(no3)), e120dat_long[is.finite(log(e120dat_long$Biomass)),]), trace=F)
summary(m2.PR)

m3.PR<-stepfun_use(lm(log(Biomass)~log(soilC)*log(NumSp)*(log(abv)+logit(ptisn)+log(no3))^2, e120dat_long[is.finite(log(e120dat_long$Biomass)),]), trace=F)
summary(m3.PR)

m4.PR<-stepfun_use(lm(log(Biomass)~log(soilC)*log(NumSp)*(log(abv)+logit(ptisn)+log(no3))^3, e120dat_long[is.finite(log(e120dat_long$Biomass)),]), trace=F)
summary(m4.PR)


nabund<-nrow(e120dat_long[is.finite(log(e120dat_long$Biomass)),])

########################################################
#test fits (species level)
########################################################

############################
#no planted richness
############################
lb<-log(e120dat_long$Biomass[is.finite(log(e120dat_long$Biomass))])

suppressMessages(lm0<-lmodel2(lb~predict(m0), range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))

suppressMessages(lm1.1<-lmodel2(lb~predict(m1.1), range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm1.2<-lmodel2(lb~predict(m1.2), range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm1.3<-lmodel2(lb~predict(m1.3), range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm2<-lmodel2(lb~predict(m2), range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm3<-lmodel2(lb~predict(m3), range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm4<-lmodel2(lb~predict(m4), range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))

1-(1-lm0$rsquare)*(nabund-1)/(nabund-length(coef(m0))-1)

1-(1-lm1.1$rsquare)*(nabund-1)/(nabund-length(coef(m1.1))-1)
1-(1-lm1.2$rsquare)*(nabund-1)/(nabund-length(coef(m1.2))-1)
1-(1-lm1.3$rsquare)*(nabund-1)/(nabund-length(coef(m1.3))-1)

1-(1-lm2$rsquare)*(nabund-1)/(nabund-length(coef(m2))-1)
1-(1-lm3$rsquare)*(nabund-1)/(nabund-length(coef(m3))-1)
1-(1-lm4$rsquare)*(nabund-1)/(nabund-length(coef(m4))-1)

############################
#with planted richness
############################
suppressMessages(lm0.PR<-lmodel2(lb~predict(m0.PR), range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm1.1.PR<-lmodel2(lb~predict(m1.1.PR), range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm1.2.PR<-lmodel2(lb~predict(m1.2.PR), range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm1.3.PR<-lmodel2(lb~predict(m1.3.PR), range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm2.PR<-lmodel2(lb~predict(m2.PR), range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm3.PR<-lmodel2(lb~predict(m3.PR), range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm4.PR<-lmodel2(lb~predict(m4.PR), range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))

1-(1-lm0.PR$rsquare)*(nabund-1)/(nabund-length(coef(m0.PR))-1)

1-(1-lm1.1.PR$rsquare)*(nabund-1)/(nabund-length(coef(m1.1.PR))-1)
1-(1-lm1.2.PR$rsquare)*(nabund-1)/(nabund-length(coef(m1.2.PR))-1)
1-(1-lm1.3.PR$rsquare)*(nabund-1)/(nabund-length(coef(m1.3.PR))-1)

1-(1-lm2.PR$rsquare)*(nabund-1)/(nabund-length(coef(m2.PR))-1)
1-(1-lm3.PR$rsquare)*(nabund-1)/(nabund-length(coef(m3.PR))-1)
1-(1-lm4.PR$rsquare)*(nabund-1)/(nabund-length(coef(m4.PR))-1)

########################################################
#test fits (total biomass)
########################################################

############################
#no planted richness
############################
lb<-log(e120dat_long$Biomass[is.finite(log(e120dat_long$Biomass))])
pltnum<-e120dat_long$Plot[is.finite(log(e120dat_long$Biomass))]

lb_plttot<-tapply(lb, pltnum, function(x) log(sum(exp(x))))

pm0_plttot<-tapply(predict(m0), pltnum, function(x) log(sum(exp(x))))

pm1.1_plttot<-tapply(predict(m1.1), pltnum, function(x) log(sum(exp(x))))
pm1.2_plttot<-tapply(predict(m1.2), pltnum, function(x) log(sum(exp(x))))
pm1.3_plttot<-tapply(predict(m1.3), pltnum, function(x) log(sum(exp(x))))

pm2_plttot<-tapply(predict(m2), pltnum, function(x) log(sum(exp(x))))
pm3_plttot<-tapply(predict(m3), pltnum, function(x) log(sum(exp(x))))
pm4_plttot<-tapply(predict(m4), pltnum, function(x) log(sum(exp(x))))

suppressMessages(lm0_tot<-lmodel2(lb_plttot~pm0_plttot, range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))

suppressMessages(lm1.1_tot<-lmodel2(lb_plttot~pm1.1_plttot, range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm1.2_tot<-lmodel2(lb_plttot~pm1.2_plttot, range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm1.3_tot<-lmodel2(lb_plttot~pm1.3_plttot, range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))

suppressMessages(lm2_tot<-lmodel2(lb_plttot~pm2_plttot, range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm3_tot<-lmodel2(lb_plttot~pm3_plttot, range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm4_tot<-lmodel2(lb_plttot~pm4_plttot, range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))

ntot<-length(unique(pltnum))

1-(1-lm0_tot$rsquare)*(nabund-1)/(nabund-length(coef(m0))-1)

1-(1-lm1.1_tot$rsquare)*(nabund-1)/(nabund-length(coef(m1.1))-1)
1-(1-lm1.2_tot$rsquare)*(nabund-1)/(nabund-length(coef(m1.2))-1)
1-(1-lm1.3_tot$rsquare)*(nabund-1)/(nabund-length(coef(m1.3))-1)

1-(1-lm2_tot$rsquare)*(nabund-1)/(nabund-length(coef(m2))-1)
1-(1-lm3_tot$rsquare)*(nabund-1)/(nabund-length(coef(m3))-1)
1-(1-lm4_tot$rsquare)*(nabund-1)/(nabund-length(coef(m4))-1)

############################
#with planted richness
############################

pm0_plttot.PR<-tapply(predict(m0.PR), pltnum, function(x) log(sum(exp(x))))

pm1.1_plttot.PR<-tapply(predict(m1.1.PR), pltnum, function(x) log(sum(exp(x))))
pm1.2_plttot.PR<-tapply(predict(m1.2.PR), pltnum, function(x) log(sum(exp(x))))
pm1.3_plttot.PR<-tapply(predict(m1.3.PR), pltnum, function(x) log(sum(exp(x))))

pm2_plttot.PR<-tapply(predict(m2.PR), pltnum, function(x) log(sum(exp(x))))
pm3_plttot.PR<-tapply(predict(m3.PR), pltnum, function(x) log(sum(exp(x))))
pm4_plttot.PR<-tapply(predict(m4.PR), pltnum, function(x) log(sum(exp(x))))


suppressMessages(lm0_tot.PR<-lmodel2(lb_plttot~pm0_plttot.PR, range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))

suppressMessages(lm1.1_tot.PR<-lmodel2(lb_plttot~pm1.1_plttot.PR, range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm1.2_tot.PR<-lmodel2(lb_plttot~pm1.2_plttot.PR, range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm1.3_tot.PR<-lmodel2(lb_plttot~pm1.3_plttot.PR, range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))

suppressMessages(lm2_tot.PR<-lmodel2(lb_plttot~pm2_plttot.PR, range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm3_tot.PR<-lmodel2(lb_plttot~pm3_plttot.PR, range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))
suppressMessages(lm4_tot.PR<-lmodel2(lb_plttot~pm4_plttot.PR, range.y = "interval", range.x = "interval", nperm = nrep_coef_boot))


1-(1-lm0_tot.PR$rsquare)*(nabund-1)/(nabund-length(coef(m0.PR))-1)

1-(1-lm1.1_tot.PR$rsquare)*(nabund-1)/(nabund-length(coef(m1.1.PR))-1)
1-(1-lm1.2_tot.PR$rsquare)*(nabund-1)/(nabund-length(coef(m1.2.PR))-1)
1-(1-lm1.3_tot.PR$rsquare)*(nabund-1)/(nabund-length(coef(m1.3.PR))-1)

1-(1-lm2_tot.PR$rsquare)*(nabund-1)/(nabund-length(coef(m2.PR))-1)
1-(1-lm3_tot.PR$rsquare)*(nabund-1)/(nabund-length(coef(m3.PR))-1)
1-(1-lm4_tot.PR$rsquare)*(nabund-1)/(nabund-length(coef(m4.PR))-1)

########################################################
#Calculate partial R2 for each variable
########################################################

modlst<-list(m1.1, m1.2, m1.3, m2, m3, m4,
             m0.PR, m1.1.PR, m1.2.PR, m1.3.PR, m2.PR, m3.PR, m4.PR)
modlst_names<-c("m1.1", "m1.2", "m1.3", "m2", "m3", "m4",
                "m0.PR", "m1.1.PR", "m1.2.PR", "m1.3.PR", "m2.PR", "m3.PR", "m4.PR")
coeflst<-c("abv", "ptisn", "no3", "soilC", "NumSp")

pr2dat<-NULL
for(i in 1:length(modlst)) {
  for(j in 1:length(coeflst)) {
    mtmp<-modlst[[i]]
    
    tmp<-names(coef(mtmp))
    rmlst<-grep(coeflst[j], tmp)
    if(length(rmlst)>0) {
      mc<-paste("update(mtmp, .~.-", paste(tmp[rmlst], collapse="-"), ")")
      mnew<-eval(parse(text = mc))
      
      SSfull<-sum((resid(summary(mtmp))-lb)^2)
      SSred<-sum((resid(summary(mnew))-lb)^2)
      
      pr2<-(SSfull-SSred)/SSred
      
      pr2dat<-rbind(pr2dat, data.frame(mod=modlst_names[i], coef=coeflst[j], pr2=pr2))
    }
  }
}

#remove N from null model
SSfull<-sum((resid(summary(m0.PR))-lb)^2)
SSred<-sum((resid(summary(m0))-lb)^2)
pr2_N<-(SSfull-SSred)/SSred

#remove C from null model
SSred<-sum((resid(summary(update(m0.PR, .~.-log(soilC)-log(soilC):log(NumSp))))-lb)^2)
pr2_C<-(SSfull-SSred)/SSred

########################################################
#Create summary table
########################################################

#model names
regression_sumtable<-data.frame(Model=rep(c("B*mono", "q", "R*", "all traits", "2-way interactions", "3-way interactions"), each=2), Type=rep(c("Base Model", "Planted Richness"), 6))

#number of parameters
regression_sumtable$npar<-c(length(coef(m1.1)),
                            length(coef(m1.1.PR)),
                            length(coef(m1.2)),
                            length(coef(m1.2.PR)),
                            length(coef(m1.3)),
                            length(coef(m1.3.PR)),
                            length(coef(m2)),
                            length(coef(m2.PR)),
                            length(coef(m3)),
                            length(coef(m3.PR)),
                            length(coef(m4)),
                            length(coef(m4.PR)))

#p-value of ANOVA
regression_sumtable$pvalue<-c(anova(m0, m1.1)[2,6],
  anova(m0.PR, m1.1.PR)[2,6],
  anova(m0, m1.2)[2,6],
  anova(m0.PR, m1.2.PR)[2,6],
  anova(m0, m1.3)[2,6],
  anova(m0.PR, m1.3.PR)[2,6],
  anova(m1.1, m2)[2,6],
  anova(m1.1.PR, m2.PR)[2,6],
  anova(m2, m3)[2,6],
  anova(m2.PR, m3.PR)[2,6],
  anova(m3, m4)[2,6],
  anova(m3.PR, m4.PR)[2,6])
  

#partial r2 for all coefficients
tmp<-pr2dat[pr2dat$coef=="abv",]
regression_sumtable$partialr2_Bmono<-
  tmp$pr2[match(c("m1.1", "m1.1.PR", "m1.2", "m1.2.PR", "m1.3", "m1.3.PR",
    "m2", "m2.PR", "m3", "m3.PR", "m4", "m4.PR"), tmp$mod)]

tmp<-pr2dat[pr2dat$coef=="ptisn",]
regression_sumtable$partialr2_q<-
  tmp$pr2[match(c("m1.1", "m1.1.PR", "m1.2", "m1.2.PR", "m1.3", "m1.3.PR",
                  "m2", "m2.PR", "m3", "m3.PR", "m4", "m4.PR"), tmp$mod)]

tmp<-pr2dat[pr2dat$coef=="no3",]
regression_sumtable$partialr2_R<-
  tmp$pr2[match(c("m1.1", "m1.1.PR", "m1.2", "m1.2.PR", "m1.3", "m1.3.PR",
                  "m2", "m2.PR", "m3", "m3.PR", "m4", "m4.PR"), tmp$mod)]

tmp<-pr2dat[pr2dat$coef=="soilC",]
regression_sumtable$partialr2_soilC<-
  tmp$pr2[match(c("m1.1", "m1.1.PR", "m1.2", "m1.2.PR", "m1.3", "m1.3.PR",
                  "m2", "m2.PR", "m3", "m3.PR", "m4", "m4.PR"), tmp$mod)]

tmp<-pr2dat[pr2dat$coef=="NumSp",]
regression_sumtable$partialr2_PlRich<-
  tmp$pr2[match(c("m1.1", "m1.1.PR", "m1.2", "m1.2.PR", "m1.3", "m1.3.PR",
                  "m2", "m2.PR", "m3", "m3.PR", "m4", "m4.PR"), tmp$mod)]
                            

#p-values for base model
#NA corresponds to "not significant" (i.e. not significantly better fitting than a simpler model)
regression_sumtable$pvalue_base<-c(lm1.1$regression.results[4,5],
                                  lm1.1.PR$regression.results[4,5],
                                  lm1.2$regression.results[4,5],
                                  lm1.2.PR$regression.results[4,5],
                                  lm1.3$regression.results[4,5],
                                  lm1.3.PR$regression.results[4,5],
                                  lm2$regression.results[4,5],
                                  lm2.PR$regression.results[4,5],
                                  lm3$regression.results[4,5],
                                  lm3.PR$regression.results[4,5],
                                  lm4$regression.results[4,5],
                                  lm4.PR$regression.results[4,5])

#r2 and adj. r2 for base model
regression_sumtable$r2_base<-c(lm1.1$rsquare,
                                   lm1.1.PR$rsquare,
                                   lm1.2$rsquare,
                                   lm1.2.PR$rsquare,
                                   lm1.3$rsquare,
                                   lm1.3.PR$rsquare,
                                   lm2$rsquare,
                                   lm2.PR$rsquare,
                                   lm3$rsquare,
                                   lm3.PR$rsquare,
                                   lm4$rsquare,
                                   lm4.PR$rsquare)

regression_sumtable$r2_adj_base<-1-(1-regression_sumtable$r2_base)*(nabund-1)/(nabund-regression_sumtable$npar-1)

#p-values for planted richness model
regression_sumtable$pvalue_tot<-c(lm1.1_tot$regression.results[4,5],
                              lm1.1_tot.PR$regression.results[4,5],
                              lm1.2_tot$regression.results[4,5],
                              lm1.2_tot.PR$regression.results[4,5],
                              lm1.3_tot$regression.results[4,5],
                              lm1.3_tot.PR$regression.results[4,5],
                              lm2_tot$regression.results[4,5],
                              lm2_tot.PR$regression.results[4,5],
                              lm3_tot$regression.results[4,5],
                              lm3_tot.PR$regression.results[4,5],
                              lm4_tot$regression.results[4,5],
                              lm4_tot.PR$regression.results[4,5])

#r2 and adj. r2 for planted richness model
regression_sumtable$r2_tot<-c(lm1.1_tot$rsquare,
                               lm1.1_tot.PR$rsquare,
                               lm1.2_tot$rsquare,
                               lm1.2_tot.PR$rsquare,
                               lm1.3_tot$rsquare,
                               lm1.3_tot.PR$rsquare,
                               lm2_tot$rsquare,
                               lm2_tot.PR$rsquare,
                               lm3_tot$rsquare,
                               lm3_tot.PR$rsquare,
                               lm4_tot$rsquare,
                               lm4_tot.PR$rsquare)

regression_sumtable$r2_adj_tot<-1-(1-regression_sumtable$r2_tot)*(nabund-1)/(nabund-regression_sumtable$npar-1)

#drop extra decimal places. Note p is always <= value
regression_sumtable[,c(4,10,13)]<-ceiling(regression_sumtable[,c(4,10,13)]*1000)*0.001
regression_sumtable[,c(5:9,11:12,14:15)]<-round(regression_sumtable[,c(5:9,11:12,14:15)],3)


write.csv(regression_sumtable, "data/data_products/table1_regression_summaries.csv", row.names=F)
