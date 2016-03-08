##########################
#load data file
##########################
e120dat<-read.csv("data/e120_spAbundance_01-14.csv") #observed e120 data

#location of biomass for each planted species
bmpos<-c(5:18) #total biomass
e120dat<-cbind(e120dat, matrix(nrow=nrow(e120dat), ncol=length(splst)))
plbmpos<-c(19:32) #empty list for total biomass, only in plots where it was planted
colnames(e120dat)[plbmpos]<-paste("Planted.", colnames(e120dat)[bmpos], sep="")

#load data on which species were planted in which plots
plotlst<-read.csv("data/e120_plantedpos.csv")

#make simple matrix of plot locations for each species
tmp<-strsplit(as.character(plotlst$Species), " ", fixed=T)
pltloc<-data.frame(matrix(nrow=length(tmp), ncol=length(splst)+1, data=0))
colnames(pltloc)<-c("plot", splst)

for(i in 1:length(tmp)) {
  pltloc[i,1]<-plotlst$Plot[i]
  ps<-match(tmp[[i]], splst)
  pltloc[i,1+ps[is.finite(ps)]]<-1
}

#achieved species richness in each plot
pltloc$srich<-rowSums(pltloc[,-1])

#planted species richness in each plot
pltloc$plantedsr<-plotlst$Number.of.species[match(pltloc$plot, plotlst$Plot)]

#if not species were planted in plot or if plot is a monoculture, set planted richness to 0 or to 1
pltloc$plantedsr[pltloc$srich==0]<-0
pltloc$plantedsr[pltloc$srich==1]<-1

#update e120dat to match pltloc
for(i in 1:nrow(e120dat)) {
  e120dat[i,plbmpos]<-0
  ps<-as.logical(pltloc[match(e120dat$Plot[i], pltloc$plot),2:15])
  e120dat[i,plbmpos[ps]]<-e120dat[i,bmpos[ps]]
}

pCstart<-e120_tn
pCstart<-pCstart[order(pCstart$Plot),]

pCmedian<-median(pCstart$pCSoil940_20)
agPC<-pCstart[pCstart$Plot%in%pltloc$plot[pltloc$srich==1],]
agPC$monoculture<-apply(pltloc[pltloc$srich==1,2:15], 1, function(x) which(x==1))
pCmed<-tapply(agPC$pCSoil940_20, agPC$monoculture, function(x) exp(mean(log(x[is.finite(log(x))]), na.rm=T)))

if(centermeans) {
  #Re-center monoculture means on E120 values
  e120dat$sr<-pltloc$srich[match(e120dat$Plot, pltloc$plot)]
  e120monodat<-e120dat[e120dat$sr==1,]
  e120monodat$monoculture<-NA
  e120monodat$monoculturebiomass<-NA
  for(i in 1:length(splst)) {
    sbs<-pltloc$plot[pltloc$sr==1&pltloc[,i+1]==1]
    sbs<-which(e120monodat$Plot%in%sbs)
    e120monodat$monoculture[sbs]<-splst[i]
    e120monodat$monoculturebiomass[sbs]<-e120monodat[sbs,plbmpos[i]]
  }
  
  e120monodat<-e120monodat[e120monodat$Plot%in%pltloc$plot[pltloc$plantedsr==1],]
  abv<-with(aggregate(cbind(obs=e120monodat$monoculturebiomass), by=list(plot=e120monodat$Plot, mono=e120monodat$monoculture), hurdlemodel),
            aggregate(cbind(obs=obs), list(mono=mono), hurdlemodel))$obs
  e120_abv<-abv
  
  if(adjustS) {
    #if biomass will be adjusted to match plot-level soil fertility,
    #scale centered values such that monocultures will still match observations
    for(i in 1:length(splst)) {
      lap<-pltloc$plot[as.logical((pltloc$plantedsr==1)&(pltloc[,i+1]==1))]
      slw<-(pCstart$pCSoil940_20[pCstart$Plot%in%lap])/pCmedian
      
      wi<-(1/prod(slw))^(1/length(slw))
      tre120$abv[i]<-abv[i]*wi
    }
  }
}

if(usetr==2) {
  #fit trait data to tradeoff
  tottr<-nondirfit(vardf = data.frame(no3=log(trdat$no3), ptisn=logit(trdat$ptisn), abv=log(trdat$abv)))#, slp=trdat$tmp_slope))
  #extract positions snapped to tradeoff surface
  ndf<-tottr$possnap
  tre120$no3<-exp(ndf[,"no3"])[trdat$ine120]
  
  if(!centermeans) {
    tre120$abv<-exp(ndf[,"abv"])[trdat$ine120]
  }
  
  #if means were centered, update ptisn to maintain same B*q as before
  tmp<-exp(ndf[,"abv"])[trdat$ine120]*ilogit(ndf[,"ptisn"])[trdat$ine120]
  tre120$ptisn<-tmp/tre120$abv
}

#remove Amorpha canescens records from 2-8 species plots, as it was not planted there
pltloc$Amoca[(pltloc$plantedsr>1)&(pltloc$plantedsr<16)]<-0

#calculate achieved species richness in each plot
pltloc$srich<-rowSums(pltloc[,2:15])

#use all monocultures as monocutlures
pltloc$plantedsr[pltloc$srich==1]<-1
plts<-pltloc$plot[pltloc$srich>0]

#make empty list to save outputs
datout<-NULL
ssout<-NULL


