#!/bin/R

# plot PCA correction figures for splcing data

#Libraries
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(sva)
library(cowplot)


# Master directory
dir = Sys.getenv('RARE_DIS_DIR')

# get PCA data
load("PCA_splicing.RData")

# Get metadata
metadata = paste(dir,"/data/metadata/2018_06_12_Rare_Disease_Metadata.tsv",sep="")
metadata = read_tsv(metadata) %>% filter(in_freeze=="yes") %>% filter(source_of_RNA=="Blood")


meta_pivus=paste(dir,"/data/metadata/PIVUS_RNASequencingInfo.csv", sep="")
meta_pivus=read_csv(meta_pivus) %>% filter(Age=="70") %>% filter(RunOK=="Yes")
meta_pivus2=paste(dir,"/data/metadata/Pivus_expCovariatesAll_unscaled.txt", sep="")
meta_pivus2=read_tsv(meta_pivus2)%>% filter(Age==0)

# get useful info for pca in data frame
data_pca=data.frame(sample=rownames(res.comp$completeObs), 
	institution=c(metadata$institution, rep("PIVUS", 65)), 
	batch=as.factor(as.integer(as.factor(c(metadata$batch,meta_pivus$RUN)))), 
	concentration=c(metadata$bioanalyzer_cDNA_conc_pM, meta_pivus$RNA_conc_pM), 
	RIN=c(metadata$RIN,meta_pivus$RIN))

# get pca results for splicing data with main pcs regressed out
pca_regressed=prcomp(fitsv$residuals)

# handle covariates info 
covariates=data.frame(
	concentration=as.numeric(c(metadata$bioanalyzer_cDNA_conc_pM, meta_pivus$RNA_conc_pM)),
	RIN=c(metadata$RIN,meta_pivus$RIN),
	age=c(metadata$age, meta_pivus$Age), 
	sex=c(as.integer(as.factor(metadata$sex)) -1, meta_pivus2$Sex),
	affected_status=as.numeric(as.factor(c(metadata$affected_status,rep("Control",65)))),
	read_length=c(as.integer(as.factor(metadata$read_length)),rep(1,65)))

batch=as.factor(as.integer(as.factor(c(metadata$batch,meta_pivus$RUN))))
batch_mat=model.matrix(~batch)
batch_mat=as.data.frame(batch_mat[,2:ncol(batch_mat)])
batch_mat$batch1=c(ifelse(metadata$batch==1,1,0), rep(0,65))


institution=as.factor(c(metadata$institution, rep("PIVUS", 65)))
inst_mat=model.matrix(~institution)
inst_mat=as.data.frame(inst_mat[,2:ncol(inst_mat)])
inst_mat$institutionCGS=c(ifelse(metadata$institution=="CGS",1,0), rep(0,65))
colnames(inst_mat)=c("CHEO","PIVUS","UDN","CGS")

covariates=cbind(covariates,batch_mat,inst_mat)

# get correlation matrix between pcs and covariates 
cor.pc.cov=cor(pcs_selec[,1:10],covariates,use="complete.obs")

cor.pc.cov.m=melt(cor.pc.cov)
cor.pc.cov.m$Var2=factor(cor.pc.cov.m$Var2,levels=c("UDN", "CGS", "CHEO","PIVUS","batch1", "batch2", "batch3", "batch4", "batch5", "batch6", "batch7", "batch8", "batch9", "batch10", "batch11", "batch12", "concentration", "RIN", "read_length","age","sex", "affected_status"))


save.image(file = paste(dir,"/data/FigureS9.in.RData",sep=""))

