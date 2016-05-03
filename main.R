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
nrep<-100#20000 #number of iterations for nonparametric analyses
adjustS<-TRUE #use 1994 total soil C to adjust among-plot variability?
nrep_traits<-1 #number of iterations for testing within-species trait variation. If 1, then mean values are used

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
pdf("figures/FigureS1_BivariatePlots.pdf",width=6, height=6)
  par(mfcol=c(3,2), mar=c(3,3,2,2), oma=c(1,1,0,0))
  source("make_bivariate_plot.R")
dev.off()

############################################################
# run simulations
############################################################
source("run_simulations.R")
#save.image("data/data_products/simulated_results.RData") #save output for long simulations
#load("data/data_products/simulated_results.RData")

############################################################
# plot outputs
############################################################
#Get plots of prediction fits and CD
pdf("figures/Figure2_fit_figure.pdf", width=8, height=8)
  source("aggregate_data.R")
  source("get_rsquared_intervals.R")
dev.off()

#Fit by diversity level
pdf("figures/FigureS2_species_fits.pdf", width=6, height=10)  
  par(mfrow=c(7,2), mar=c(2,2,2,2), oma=c(3,3,0,0))
  source("plot_species_level.R")
dev.off()

pdf("figures/FigureS3_richness_fits.pdf", width=8, height=8)
  source("richness_level_plots.R")
dev.off()

############################################################
# simulate from tradeoff surface
############################################################
pdf("figures/Figure3_simulated_community.pdf", width=8, height=4)
  source("simulate_communities.R")
dev.off()

#save.image("data/data_products/simulated_results_simulated.RData") #save output for long simulations


############################################################
# output supplementary tables
############################################################
source("make_sup_tables.R")


