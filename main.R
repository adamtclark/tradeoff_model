#!/usr/bin/env Rscript
setwd("~/Dropbox/Active Work/Projects/001_rescomp_thesis/tradeoff_model/")
rm(list=ls())

############################################################
# set up for simulations
############################################################
#load packages
require(parallel) #parallelize calculations
require(lmodel2) #fit bivariate major axis regressions
require(mvtnorm) #calculate multivariate normal distribution
require(plot3D) #make tradeoff plots
require(lme4) #fit trait values from data
require(data.table) #for reforming data

#Load R functions
source("get_filtered_estimate_functions.R")

#build data for analyses
source("filter_tradeoff_data.R")

#lists for outputs
datoutlst<-NULL #list for output means
ssoutlst<-NULL #list for output total data

#set preferences
bootr2<-TRUE #get bootstrapped estimates for R2?
dotradeofftest<-"saved" #get p-values for tradeoff slopes. Value "saved" loads an existing file

#simulation options
centermeans<-TRUE #center to true E120 monoculture means?
nrep<-1000 #number of iterations for nonparametric analyses
adjustS<-TRUE #use 1994 total soil C to adjust among-plot variability?
nrep_traits<-1000 #number of iterations for testing within-species trait variation. If 1, then mean values are used

#Load C code
system("R CMD SHLIB getbmest.c")
dyn.load("getbmest.so")
if(!exists("cl") & nrep_traits>1) {
  cl <- makeCluster(mc <- getOption("cl.cores", detectCores())) #cluters for simulations
}

############################################################
#Make 3D tradeoff plot
############################################################
pdf("figures/Figure1_tradeoff_surface_3D.pdf", width=6, height=6)
  source("make_3d_plot.R")
dev.off()

############################################################
#Plot bi-variate relationships
############################################################
pdf("figures/FigureS2_BivariatePlots.pdf",width=6, height=6)
  par(mfcol=c(3,2), mar=c(3,3,2,2), oma=c(1,1,0,0))
  source("make_bivariate_plot.R")
dev.off()

############################################################
# run simulations
############################################################
source("run_simulations.R")
save.image("data/data_products/simulated_results.RData") #save output for long simulations
#load("data/data_products/simulated_results.RData")

############################################################
# plot outputs
############################################################
  #Get plots of prediction fits and CD
  pdf("figures/Figure2_fit_figure.pdf", width=8, height=12)
    source("aggregate_data.R") #observed vs. fitted
    source("get_rsquared_intervals.R") #MAE fit for biomass by diversity level
    source("get_richness_metrics.R") #observed vs. sampled diversity
  
    #make figure labels
    mtext("Planted Richness", 1, line=-44.5, cex=1.5, outer=T)
    mtext(expression(paste("Observed Biomass, g m"^"-2", sep="")), 2, line=0.4, cex=1.5, outer=T, adj=0.18)
    mtext(expression(paste("Predicted Biomass, g m"^"-2", sep="")), 1, line=2, cex=1.5, outer=T, adj=0.5)
  
    mtext(expression(paste("MAE, fold change", sep="")), 2, line=0.4, cex=1.5, outer=T, adj=0.68)
    mtext("Sample Richness", 2, line=0.85, cex=1.5, outer=T, adj=0.965)
  
    mtext("With Snapping", 4, line=-1, cex=1, outer=T, adj=0.081)
    mtext("Without Snapping", 4, line=-1, cex=1, outer=T, adj=0.33)
    
    adj<-0.19
    mtext("Species Abundance", 1, line=-40, cex=1, outer=T, adj=adj)
    mtext("Net Primary Productivity", 1, line=-40, cex=1, outer=T, adj=1-adj)
    
    mtext("Species Abundance", 1, line=-65.5, cex=1, outer=T, adj=adj)
    mtext("Net Primary Productivity", 1, line=-65.5, cex=1, outer=T, adj=1-adj)
  
    mtext("Without Intraspecific Variation", 1, line=-86, cex=1, outer=T, adj=adj-0.04)
    mtext("With Intraspecific Variation", 1, line=-86, cex=1, outer=T, adj=1-adj+0.04)
  dev.off()

############################################################
# simulate from tradeoff surface
############################################################
pdf("figures/Figure3_simulated_community.pdf", width=8, height=4)
  source("simulate_communities.R")
dev.off()
save.image("data/data_products/simulated_results_simulated.RData") #save output for long simulations

############################################################
# output supplementary tables
############################################################
source("make_sup_tables.R")

############################################################
#Get interaction coefficients
############################################################
source("test_biomass_cor.R")

############################################################
# run altered simulations
############################################################
#refresh cluster
stopCluster(cl)
cl <- makeCluster(mc <- getOption("cl.cores", detectCores())) #cluters for simulations
datout_alteredlst<-NULL #list for saving results

#Andropogon gerardi
alter_whichspecies<-which(splst=="Andge") #alter R* for which species?
alter_whichno3<-which(splst%in%c("Schsc")) #alter it to lower than which species?
alter_whichrichness<-c(8,16) #for which richness treatments?
alter_probability<-0.5 #with what probability?

source("run_simulations_altered.R")
datout_alteredlst[[1]]<-datout_altered

#Lupinus perennis
alter_whichspecies<-which(splst%in%c("Luppe"))
alter_whichno3<-which(splst=="Amoca")
alter_whichrichness<-c(1,2,4,8,16)
alter_probability<-0.8

source("run_simulations_altered.R")
datout_alteredlst[[2]]<-datout_altered

#make plots
datout_altered_andge<-datout_alteredlst[[1]]
datout_altered_luppe<-datout_alteredlst[[2]]

pdf("figures/Figure4_Augmented_models.pdf", width=8, height=4)
  source("plot_adjustments.R")
dev.off()

pdf("figures/FigureS3_seasonality_effects.pdf", width=6, height=10)  
  source("plot_early_season.R")
dev.off()

############################################################
#close cluster
############################################################
if(exists("cl")) {
  stopCluster(cl)
}

save.image("data/data_products/simulated_results_altered.RData") #save output for long simulations