usetr<-2

#load tradeoff data
trdat<-read.csv("data/data_products/filtered_tradeoff_data.csv")

tre120<-trdat[trdat$ine120==TRUE,] #make list for species in E120

splst<-c("Achmi", "Amoca", "Andge", "Asctu", "Koecr", "Lesca", "Liaas",
         "Luppe", "Panvi", "Petpu", "Poapr", "Schsc", "Solri", "Sornu")
fglst<-as.character(tre120$fg)

#Load e120 data
source("load_e120_data.R")

#make predictions for each plot
for(i in 1:length(plts)) {
  #which species are in this plot?
  ps<-as.logical(pltloc[match(plts[i], pltloc$plot),2:15])
  #how many species are in this plot?
  nsp<-sum(ps)
  #load trait data for these species
  no3lst_dat<-tre120[which(ps),c("no3", "sd_no3")]
  pNi_dat<-tre120[which(ps),c("ptisn", "sd_ptisn")]
  abmi_dat<-tre120[which(ps),c("abv", "sd_abv")]
  
  #fix poa pretensis in non-monoculture plots
  if(!exists("lowerdlim")) {
    lowerdlim<-1
  }
  if(nsp>lowerdlim) {
    if(usetr==2) {
      abmi_dat[splst[ps]=="Poapr",1]<-exp(ndf[,"abv"])[trdat$Species=="Poa pratensis"]
      pNi_dat[splst[ps]=="Poapr",1]<-exp(ndf[,"ptisn"])[trdat$Species=="Poa pratensis"]
    } else {
      abmi_dat[splst[ps]=="Poapr",1]<-trdat$abv[trdat$Species=="Poa pratensis"]
      pNi_dat[splst[ps]=="Poapr",1]<-trdat$ptisn[trdat$Species=="Poa pratensis"]
    }
  }
  
  #how many years of observations are there in this plot?
  niter<-sum(e120dat$Plot==plts[i])
  #which years are these?
  yrs<-sort(unique(e120dat$Year[e120dat$Plot==plts[i]]))
  
  #matrix for storing prediction output for this plot
  esti<-matrix(nrow=nrep_traits,ncol=NROW(abmi_dat))
  esti_sd<-matrix(nrow=nrep_traits,ncol=NROW(abmi_dat))
  
  if(adjustS) {
    dS<-(pCstart$pCSoil940_20[which(plts[i]==pCstart$Plot)])/(pCmedian)
    
    abmi_dat$abv<-(abmi_dat$abv*dS)
  }
  
  plrich<-pltloc$plantedsr[which(plts[i]==pltloc$plot)]
  if(nrep_traits>1) { #Simulate error distribiutions if nrep>1
    #export parameters for parallel function
    clusterExport(cl, c("niter", "no3lst_dat", "pNi_dat", "abmi_dat",
                        "colSums_safe", "colMeans_safe", "colSD_safe", "nsp", "ilogit", "logit",
                        "ps", "plrich",
                        "alter_whichspecies", "alter_whichno3", "alter_whichrichness", "alter_probability"))
    
    #run parallel program for predicting community biomass
    clusterout<-t(matrix(nrow=2*nsp, unlist(parLapply(cl=cl, 1:nrep_traits, fun=repsmp_altered))))
    
    #extract results
    esti<-clusterout[,1:nsp]
    esti_sd<-clusterout[,(nsp+1):ncol(clusterout)]
    
    esti[!is.finite(esti)]<-0
    esti_sd[(!is.finite(esti_sd))|(esti_sd==0)]<-NA
    
    #apply hurdle model to take mean across years
    pzero<-colSums_safe(esti==0)/NROW(esti)
    est<-exp(colMeans_safe(log(esti)))*(1-pzero)
    est_sd<-colMeans_safe(esti_sd)
    est_sd[!is.finite(est_sd)]<-0
  } else { #Just run once
    clusterout<-repsmp_single()
    
    est<-esti<-clusterout[1:nsp]
    est_sd<-clusterout[(nsp+1):length(clusterout)]
  }
  
  #extract observed biomass from e120
  tmpobs<-as.matrix(e120dat[which(e120dat$Plot==plts[i]),plbmpos[which(ps)]])
  pzero<-colSums_safe(tmpobs==0)/NROW(tmpobs)
  if(nsp>1) {
    obs<-unname(exp(colMeans_safe(log(tmpobs))))*(1-pzero)
  } else {
    obs<-unname(exp(colMeans_safe(log(tmpobs))))
  }
  obs_sd<-unname(colSD_safe(log(tmpobs)))
  obs_sd[!is.finite(obs_sd)]<-0
  
  if(bootr2) {
    #extract R2 estimates
    ssout<-rbind(ssout, data.frame(plot=plts[i],
                                   sp=splst[ps],
                                   obs=unname(obs),
                                   est=c(t(esti)),
                                   iter=rep(1:nrep_traits, each=length(obs))))
  }
  
  #save mean results
  tmp<-data.frame(sp=as.character(splst[ps]), est=est, est_sd=est_sd, obs=obs, obs_sd=obs_sd, sr=nsp, plt=plts[i])
  tmp$sp<-as.character(tmp$sp)
  row.names(tmp)<-NULL
  datout<-rbind(datout, tmp)
  
  #print progress
  if(i/10 == round(i/10)) {
    print(round(i/length(plts), 2))
  }
}

#record simulation output
datout$plantedsr<-pltloc$plantedsr[match(datout$plt, pltloc$plot)]
datout_altered<-datout
ssout_altered<-ssout