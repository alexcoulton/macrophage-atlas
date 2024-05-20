source('~/work/ucl/scripts/misc/functions.R')
library(rhdf5)
library(Seurat)
library(gridExtra)
library(dplyr)
library(Matrix)

############################
#LOAD DATA 
############################

krishna.dat = readRDS('~/work/ucl/bigdata/ucl.projects/macrophage/raw.study.data/krishna/ccRCC_6pat_Seurat')
krishna.dat = UpdateSeuratObject(krishna.dat)


krishna.dat
