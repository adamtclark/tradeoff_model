##########################
#get community biomass
##########################
#R version of .C function
getbmest<-function(no3lst, pNi, abmi) {
  abm_esti<-numeric(length(abmi))
  
  subsi<-no3lst==min(no3lst, na.rm=T)
  abm_esti[subsi]<-abmi[subsi]
  
  pNuptake<-pNi
  
  ord<-order(no3lst)
  ord<-ord[!(ord%in%which(subsi))]
  
  if(sum(ord)>0) {
    for(j in 1:length(ord)) {
      abm_esti[ord[j]]<-abmi[ord[j]]-sum((pNuptake*abm_esti)/pNi[ord[j]])
      abm_esti[!is.finite(abm_esti)]<-0
      if(abm_esti[ord[j]]<0) {
        abm_esti[ord[j]]<-0
      }
    }
  }
  return(abm_esti)
}

#parallelized .C function
repsmp<-function(...) {
  dyn.load("getbmest.so")
  abvest<-matrix(ncol=NROW(abmi_dat), nrow=niter)
  for(j in 1:niter) {
    no3lst<-exp(rnorm(nsp, log(no3lst_dat[,1]), no3lst_dat[,2]))
    pNi<-ilogit(rnorm(nsp, logit(pNi_dat[,1]), pNi_dat[,2]))
    abmi<-exp(rnorm(nsp, log(abmi_dat[,1]), abmi_dat[,2]))
    
    abvest[j,]<-.C("getbmest", no3lst=as.double(no3lst), pNi=as.double(pNi), abmi=as.double(abmi), plabmi=as.integer(length(abmi)), abm_esti=as.double(numeric(length = length(abmi))))$abm_esti
    
  }
  pzero<-colSums_safe(abvest==0)/NROW(abvest)
  esti<-exp(colMeans_safe(log(abvest)))*(1-pzero)
  esti_sd<-colSD_safe(log(abvest))
  
  return(c(esti, esti_sd))
}

#run once in .C framework
repsmp_single<-function(...) {
  dyn.load("getbmest.so")
  no3lst<-no3lst_dat[,1]
  pNi<-pNi_dat[,1]
  abmi<-abmi_dat[,1]
  
  esti<-.C("getbmest", no3lst=as.double(no3lst), pNi=as.double(pNi), abmi=as.double(abmi), plabmi=as.integer(length(abmi)), abm_esti=as.double(numeric(length = length(abmi))))$abm_esti
  esti_sd<-rep(NA, length(esti))
  
  return(c(esti, esti_sd))
}

##########################
#logit and ilogit functions
##########################
logit<-function(x) {
  suppressWarnings(res<-(-log(1/x-1)))
  res[!is.finite(res)]<-NA
  res
}

ilogit<-function(x) {
  1/(1+exp(-x))
}

##########################
#matrix functions for 1+D matrices
##########################
colMeans_safe<-function(x) {
  x[!is.finite(x)]<-NA
  if((NCOL(x)==1)|(NROW(x)==1)) {
    mean(x, na.rm=T)
  } else {
    colMeans(x, na.rm=T)
  }
}

colSums_safe<-function(x) {
  if((NCOL(x)==1)|(NROW(x)==1)) {
    sum(x, na.rm=T)
  } else {
    colSums(x, na.rm=T)
  }
}

colSD_safe<-function(x) {
  x[!is.finite(x)]<-NA
  if((NCOL(x)==1)|(NROW(x)==1)) {
    sd(x, na.rm=T)
  } else {
    apply(x, 2, function(x) sd(x, na.rm=T))
  }
}

##########################
#hurdle model function
##########################
hurdlemodel<-function(x) {
  exp(mean(log(x[is.finite(log(x))]), na.rm=T))*(sum(is.finite(x) & x!=0)/length(x))
}

##########################
#major axis regression
##########################
nondirfit<-function(vardf, doscale=TRUE) {
  #vardf is a dataframe with columns of variables to be fit nodirectionally
  
  var<-vardf
  if(ncol(var)<=1 | nrow(var)<=1) {
    stop("too few variables to regress")
  }
  
  if(is.null(colnames(var))) {
    if(ncol(var)>3)
      stop("too many columns - don't even think about it")
    colnames(var)<-letters(1:ncol(var))
  }
  
  #scale data
  if(doscale) {
    sdmat<-apply(var, 2, function(x) sd(x, na.rm=T))
    mnmat<-colMeans(var, na.rm=T)
  } else {
    sdmat<-rep(1, ncol(var)); names(sdmat)=colnames(var)
    mnmat<-rep(0, ncol(var)); names(mnmat)=colnames(var)
  }
  
  varsc<-data.frame(t((t(var)-mnmat)/sdmat))
  
  vn<-names(var)
  
  #fit linear model to get initial guess of parameters
  expr<-paste("lm(", vn[1], "~", paste(vn[-1], collapse="+"), ", data=varsc)", sep="")    
  modc<-coef(eval(parse(text=expr)))
  modc[vn[1]]<-(-1)
  
  #optimize nondirectionally
  optc<-optim(par = modc, fn = optfun, gr = NULL, extraparms=list(varsc=varsc, vn=vn))$par
  optc<-optc/optc[2]
  
  #Calcuate predicted (and unscaled) values
  hatv<-matrix(nrow=nrow(varsc), ncol=ncol(varsc)); colnames(hatv)<-vn
  hatvsc<-matrix(nrow=nrow(varsc), ncol=ncol(varsc)); colnames(hatvsc)<-vn
  for(i in 1:ncol(hatv)) {
    pos<-which(names(optc)==vn[i])
    hatvsc[,i]<-(optc[-c(1, pos)]%*%t(varsc[-i])+optc[1])/(-optc[pos])
    hatv[,i]<-(hatvsc[,i]*sdmat[i])+mnmat[i]
  }
  hatv<-data.frame(hatv)
  hatvsc<-data.frame(hatvsc)
  
  #project variables onto fit line or plane
  #given ax+by+cz+...+int = 0, and point {x0, y0, z0...}
  #k = (ax0+by0+cz0+...+int)/(a^2+b^2+c^2...)
  #xnew = x0-a*k
  
  k<-(optc[c(1, match(vn, names(optc)))]%*%t(cbind(1, varsc)))/sum(optc[-1]^2)
  possnapsc<-t(t(varsc)-optc[match(vn, names(optc))]*t(matrix(ncol=length(vn), nrow=length(k), k)))
  possnap<-array(dim=dim(possnapsc))
  
  for(i in 1:ncol(possnapsc)) {
    possnap[,i]<-possnapsc[,i]*sdmat[i]+mnmat[i]
  }
  
  #Transform parms into non-scaled space
  optc_notsc<-numeric(length(optc)); names(optc_notsc)<-names(optc)
  optc_notsc[1]<-optc[1]
  for(i in 2:length(optc)) {
    scpos<-which(names(sdmat)==names(optc[i]))
    optc_notsc[i]<-optc[i]/sdmat[scpos]
    optc_notsc[1]<-optc_notsc[1]-mnmat[scpos]/sdmat[scpos]*optc[i]
  }
  
  colnames(possnap)<-colnames(var)
  
  #Get estimate of model fit
  SSres<-sum((possnapsc-varsc)^2, na.rm=T)
  SStot<-sum((t(varsc)-colMeans(varsc, na.rm=T))^2, na.rm=T) #NB - ybar ~ 0
  rsq_est<-(1-SSres/SStot)
  p<-ncol(varsc); n<-nrow(varsc)
  rsq_est_adj<-rsq_est-(1-rsq_est)*p/(n-p-1)
  
  return(list(possnap=possnap, pred=hatv, vars=var, pars=optc_notsc, rsq=list(rsq_est=rsq_est, rsq_est_adj=rsq_est_adj),
              scl=list(sdmat=sdmat, mnmat=mnmat, possnapsc=possnapsc, predsc=hatvsc, varsc=varsc, parssc=optc)))
  #regular and scaled lists
  #possnap is values snapped to tradeoff
  #hatv is predicted values
  #pars is pars for nonparametric regression (in scaled space!!)
  #scl is scaling parameters for xscal = (x-mean(x))/sd(x)
}


optfun<-function(modc, extraparms) {
  k<-(modc[c(1, match(extraparms$vn, names(modc)))]%*%t(cbind(1, extraparms$varsc)))/sum(modc[-1]^2)
  possnapsc<-t(t(extraparms$varsc)-modc[match(extraparms$vn, names(modc))]*t(matrix(ncol=length(extraparms$vn), nrow=length(k), k)))
  
  diffobs<-sum((rowSums((possnapsc-extraparms$varsc)^2)), na.rm=T)
  
  return(diffobs)
}

get_rsq_nondir<-function(trout, colnums) {
  possnapsc<-trout$scl$possnapsc[,colnums]
  varsc<-trout$scl$varsc[,colnums]
  
  SSres<-sum((possnapsc-varsc)^2, na.rm=T)
  SStot<-sum((t(varsc)-colMeans(as.matrix(varsc),na.r=T))^2, na.rm=T) #NB - ybar ~ 0
  rsq_est<-(1-SSres/SStot)
  p<-ncol(as.matrix(varsc)); n<-nrow(as.matrix(varsc))
  rsq_est_adj<-rsq_est-(1-rsq_est)*p/(n-p-1)
  
  return(list(rst_est=rsq_est, rsq_est_adj=rsq_est_adj))
}

##########################
#plot figure legend letters
##########################
#from https://waterprogramming.wordpress.com/2015/12/
put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, 
                           offset=c(0, 0), ...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.015,0.98),
                     topcenter = c(0.5525,0.98),
                     topright = c(0.985, 0.98),
                     bottomleft = c(0.015, 0.02), 
                     bottomcenter = c(0.5525, 0.02), 
                     bottomright = c(0.985, 0.02),
                     c(0.015, 0.98) )
  } else {
    coords <- c(x,y)
  }
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  text(labels=label[1], x=this.x, y=this.y, xpd=NA, ...)
}

