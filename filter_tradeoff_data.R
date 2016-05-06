############################################################
#Cleaning and fitting species trait data for e120 polyculture model
############################################################
mondatall<-read.csv("data/total_monoculture_data.csv")
dattot<-mondatall

#CLEANING RULES:
#0: Merge experiment numbers
#1: Remove high total soil N
#2: Remove high leaf litter in non-burned plots (E26)
#3: Remove low aboveground biomass
#4: Remove years with incomplete sampling
#5: Remove woody and sedge
#6: Remove young plots

#FITTING RULES:
#1: Always remove duplicates
#2: Random effects structure: (1|exp/plot)

##############################
#CLEANING
##############################
##############
#Apply rule 0: Merge experiment numbers
##############
#E70 is a subset of E55; E249 is nested within E120
dattot$exp<-as.character(dattot$exp)
dattot$exp[dattot$exp=="e070"]<-"e055"
dattot$exp[dattot$exp=="e249"]<-"e120"

##############
#Apply rule 1: Remove high total soil N
##############
#1a. remove fertilized plot in E26
dattot<-dattot[!(dattot$exp=="e026"&dattot$plot==7),]

#1b. in E26 and E55, remove blocks with mean total soil N
#more than two standard deviations above E120
tn26<-with(dattot[dattot$exp=="e026",], tapply(totaln, plot, function(x) mean(x, na.rm=T)))
tn55<-with(dattot[dattot$exp=="e055",], tapply(totaln, plot, function(x) mean(x, na.rm=T)))

e120_tn<-read.csv("data/e120_psoilN_1994.csv")
e120<-unique(data.frame(plot=e120_tn$Plot, totaln=e120_tn$pNSoil940_20))
e120<-e120[is.finite(e120$totaln),]
mtn<-mean(e120$totaln, na.rm=T)+sd(e120$totaln, na.rm=T)

dattot<-dattot[!(dattot$exp=="e026"&dattot$plot%in%as.numeric(names(tn26[tn26>mtn]))),]
dattot<-dattot[!(dattot$exp=="e055"&dattot$plot%in%as.numeric(names(tn55[tn55>mtn]))),]

##############
#Apply rule 2: Remove high leaf litter in non-burned plots (E26)
##############
mlit<-mean(dattot[dattot$exp=="e120",]$litbiomass, na.rm=T)+sd(dattot[dattot$exp=="e120",]$litbiomass, na.rm=T)
dattot<-dattot[!((is.na(dattot$litbiomass)|dattot$litbiomass>mlit)&dattot$exp=="e026"),]

##############
#Apply rule 3: Remove low aboveground or belowground biomass
##############
dattot<-dattot[((dattot$abvbiomass>10)),]

#Remove e120 poa
dattot<-dattot[!(dattot$monoculture=="Poa pratensis"&dattot$exp=="e120"),]

##############
#Apply rule 4: Remove years with incomplete sampling
##############
#4a: remove years with only partial sampling
dattot<-dattot[!is.na(dattot$exp) & is.finite(dattot$year),]
dattot<-dattot[!(dattot$exp=="e026" & dattot$year%in%c(1989, 1991, 1993)),]

#4b: remove all plots that were sampled after weeding ceased
#E123 had grazing treatments established after 1997 sampling included here
dattot<-dattot[!((dattot$exp=="e123")&(dattot$year>1997)),]

##############
#Apply rule 5: Remove woody and sedge
##############
dattot$fg<-as.character(dattot$fg)
dattot<-dattot[dattot$fg!="W" & dattot$fg!="S",]

##############
#Apply rule 6: Remove "young" plots
##############
#2a: remove all plots that are less than 3-years old
dattot<-dattot[dattot$age>=3,]

write.csv(dattot, "data/data_products/filtered_dattot.csv", row.names=F)

##############################
#FITTING
##############################
##############
#Make species-level tradeoff data
##############
datlowN<-dattot

require(lme4)
mno3<-lmer(log(no3)~-1+monoculture+(1|exp/plot/subplot), unique(datlowN[is.finite(log(datlowN$no3)),c("exp", "plot", "subplot", "year", "month", "monoculture", "no3")]))
no3<-summary(mno3)$coefficients[,1]
sd_no3<-summary(mno3)$coefficients[,2]

mabv<-lmer(log(abvbiomass)~-1+monoculture+(1|exp/plot/subplot), unique(datlowN[is.finite(log(datlowN$abvbiomass))&datlowN$age>3,c("exp", "plot", "subplot", "year", "month", "monoculture", "abvbiomass")]))
abv<-summary(mabv)$coefficients[,1]
sd_abv<-summary(mabv)$coefficients[,2]

#none of the experiments with ptisn have subplots
ptisnmod<-lmer(logit(atissN)~-1+monoculture+(1|exp/plot), data=unique(dattot[is.finite(logit(dattot$atissN)),c("exp", "plot", "subplot", "year", "month", "monoculture", "atissN")]))
ptisn<-summary(ptisnmod)$coefficients[,1]
sd_ptisn<-summary(ptisnmod)$coefficients[,2]

tradeoffdat<-data.frame(unique(dattot[,c("monoculture", "fg")]))
colnames(tradeoffdat)<-c("Species", "fg")
tradeoffdat$Species<-as.character(tradeoffdat$Species)
tradeoffdat$fg<-as.character(tradeoffdat$fg)
tradeoffdat<-tradeoffdat[order(as.character(tradeoffdat$Species)),]

tradeoffdat$no3<-no3[match(paste("monoculture", tradeoffdat$Species, sep=""), names(no3))]
tradeoffdat$sd_no3<-sd_no3[match(paste("monoculture", tradeoffdat$Species, sep=""), names(sd_no3))]

tradeoffdat$abv<-abv[match(paste("monoculture", tradeoffdat$Species, sep=""), names(abv))]
tradeoffdat$sd_abv<-sd_abv[match(paste("monoculture", tradeoffdat$Species, sep=""), names(sd_abv))]

tradeoffdat$ptisn<-ptisn[match(paste("monoculture", tradeoffdat$Species, sep=""), names(ptisn))]
tradeoffdat$sd_ptisn<-sd_ptisn[match(paste("monoculture", tradeoffdat$Species, sep=""), names(sd_ptisn))]

tradeoffdat$no3<-exp(tradeoffdat$no3)
tradeoffdat$abv<-exp(tradeoffdat$abv)
tradeoffdat$ptisn<-ilogit(tradeoffdat$ptisn)

tradeoffdat$ine120<-FALSE
tradeoffdat$ine120[tradeoffdat$Species%in%dattot$monoculture[dattot$exp=="e120"]]<-TRUE
tradeoffdat$ine120[tradeoffdat$Species=="Poa pratensis"]<-TRUE

write.csv(tradeoffdat, "data/data_products/filtered_tradeoff_data.csv", row.names=F)

