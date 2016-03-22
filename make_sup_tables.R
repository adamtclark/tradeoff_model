#Table S1

tableS1<-data.frame(tradeoffdat[,c("Species", "no3", "sd_no3", "abv", "sd_abv", "ptisn", "sd_ptisn", "ine120")])
tableS1$e120_abv<-NA
tableS1$e120_abv[tableS1$ine120]<-e120_abv

tableS1$no3_snapped<-10^(trout$possnap[,"no3"])
tableS1$abv_snapped<-10^(trout$possnap[,"abv"])
tableS1$ptisn_snapped<-ilogit(trout$possnap[,"ptisn"])

tableS1<-tableS1[is.finite(tableS1$no3)&is.finite(tableS1$abv)&is.finite(tableS1$ptisn),]

write.csv(tableS1, "data/data_products/tableS1_paramtertable.csv", row.names=F)

#Table S2
tmp<-dattot[dattot$monoculture%in%tableS1$Species,]
index<-with(tmp, paste(monoculture, exp, plot, subplot))
tmp<-tmp[!duplicated(index),]
tableS2<-table(tmp$monoculture, tmp$exp)

write.csv(tableS2, "data/data_products/tableS2_replicates.csv", row.names=F)

#Table S3
tableS3<-rbind(parameter_intervals, pvalue=colSums(pardat<0)/nrow(pardat))
tableS3<-data.frame(tableS3, threeway=c(NA, NA, NA, pvale))

write.csv(tableS3, "data/data_products/tableS3_tradeoffparameters.csv")
