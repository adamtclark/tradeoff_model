#!/usr/bin/env Rscript
#setwd("~/Dropbox/ActiveWork/Projects/001_rescomp_thesis/tradeoff_model/")
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
require(RColorBrewer) #for mixing palettes

#Load R functions
source("get_filtered_estimate_functions.R")

#build data for analyses
source("filter_tradeoff_data.R")

#lists for outputs
datoutlst<-NULL #list for output means
ssoutlst<-NULL #list for output total data

#set preferences
bootr2<-TRUE #get bootstrapped estimates for R2?
dotradeofftest<-"saved" #get p-values for tradeoff slopes. Value "saved" loads an existing file; value "TRUE" re-runs

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
pdf("figures/Figure2_tradeoff_surface_3D.pdf", width=6, height=6, colormodel = "cmyk", useDingbats = FALSE)
  source("make_3d_plot.R")
dev.off()

############################################################
#Plot bi-variate relationships
############################################################
pdf("figures/FigureS2_bivariate_plots.pdf",width=6, height=7, colormodel = "cmyk", useDingbats = FALSE)
  m2<-cbind(c(1,2,3,7), c(4,5,6,7))
  layout(m2)
  par(mar=c(3,3,2,2), oma=c(1,1,0,0))
  
  source("make_bivariate_plot.R")
  source("make_qB_plot.R")
dev.off()

############################################################
# run simulations
############################################################
source("run_simulations.R")
save.image("data/data_products/simulated_results_coex.RData") #save output for long simulations
#load("data/data_products/simulated_results_coex.RData") #save output for long simulations

#run trait regressions (this creates Table1)
source("run_trait_regressions.R")

############################################################
# plot outputs
############################################################
#plot coexistence
pdf("figures/Figure3_coexistence_predictions.pdf", width=6, height=5, colormodel = "cmyk", useDingbats = FALSE)
  source("plot_coexistence.R")
dev.off()


#Get plots of prediction fits and CD
pdf("figures/Figure4_fit_figure.pdf", width=8, height=9, colormodel = "cmyk", useDingbats = FALSE)
  source("aggregate_data.R") #observed vs. fitted
  source("get_rsquared_intervals.R") #MAE fit for biomass by diversity level
  
  #make figure labels
  mtext("Planted Richness", 1, line=-43, cex=1.5, outer=T)
  mtext(expression(paste("Observed Biomass, g m"^"-2", sep="")), 2, line=0.4, cex=1.5, outer=T, adj=0.27)
  mtext(expression(paste("Predicted Biomass, g m"^"-2", sep="")), 1, line=2, cex=1.5, outer=T, adj=0.5)

  mtext(expression(paste("MAE, fold change", sep="")), 2, line=0.4, cex=1.5, outer=T, adj=0.94)
  
  adj<-0.19
  mtext("Species Abundance", 1, line=-39.2, cex=1, outer=T, adj=adj)
  mtext("Total Plot-Level Aboveground Biomass", 1, line=-39.2, cex=1, outer=T, adj=1-adj)
  
  mtext("Species Abundance", 1, line=-64, cex=1, outer=T, adj=adj)
  mtext("Total Plot-Level Aboveground Biomass", 1, line=-64, cex=1, outer=T, adj=1-adj)
dev.off()

pdf("figures/FigureS3_intraspecific_varaition.pdf", width=8, height=4, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(1,2), oma=c(1,1,1,0), mar=c(3,3,2,1))
  source("get_richness_metrics.R") #observed vs. sampled diversity
  mtext("Sample Richness", 2, line=-0.5, cex=1.5, outer=T)
  mtext("Planted Richness", 1, line=-0.5, cex=1.5, outer=T)
  
  mtext("Without Intraspecific Variation", 3, line=-0.5, cex=1, outer=T, adj=0.15)
  mtext("With Intraspecific Variation", 3, line=-0.5, cex=1, outer=T, adj=0.85)
dev.off()

############################################################
# simulate from tradeoff surface
############################################################
pdf("figures/Figure5_simulated_community.pdf", width=8, height=4, colormodel = "cmyk", useDingbats = FALSE)
  source("simulate_communities.R")
dev.off()
save.image("data/data_products/simulated_results_simulated_coex.RData") #save output for long simulations
#load("data/data_products/simulated_results_simulated_coex.RData") #load output for long simulations

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

pdf("figures/FigureS4_augmented_models.pdf", width=10, height=8, colormodel = "cmyk", useDingbats = FALSE)
  source("plot_adjustments.R")
dev.off()

############################################################
#close cluster
############################################################
if(exists("cl")) {
  stopCluster(cl)
}

save.image("data/data_products/simulated_results_altered_coex.RData") #save output for long simulations
#load("data/data_products/simulated_results_altered_coex.RData") #load output for long simulations
