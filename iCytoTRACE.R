#Setting up first time installation #more installation on their website https://cytotrace.stanford.edu

setwd("~/Documents/work/PhD /R")
install.packages("devtools")
library("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("sva")
library(sva)
devtools::install_local("~/Documents/work/PhD /R/CytoTRACE_0.3.3.tar")
install.packages("reticulate")
library(reticulate)
py_install("pandas")
#installation of a base version of python https://www.python.org/downloads/

#create a new conda environment
conda_create("cytotrace")
Sys.setenv(RETICULATE_PYTHON="/Users/ruth/Library/r-miniconda/envs/cytotrace/bin/python") 

# install python packages in the conda cytotrace environment with conda 
conda_install("cytotrace", "numpy")
# install python packages via pip (not available on conda)
py_install("scanoramaCT",pip=TRUE)
use_condaenv("cytotrace")
#import packages
numpy <- import("numpy")
scanoramaCT <- import("scanoramaCT")
#intall cytotrace
library(CytoTRACE)
#istall seurat
install.packages("Seurat")
library(Seurat)

#loading of packages for the second time
library(devtools)
library(sva)
library(reticulate)
py_install("scanoramaCT",pip=TRUE)
use_condaenv("cytotrace")
numpy <- import("numpy")
scanoramaCT <- import("scanoramaCT")
library(CytoTRACE)

#Read in data
brain <- readRDS("seurat_object_with_annotation.rds")
brain[["annotation"]] <- Idents(object = brain)
brain <- SetIdent(brain, value = brain@meta.data$orig.ident)
onlycells <- subset(brain, idents = "cells")
onlynuclei <- subset(brain, idents = "nuclei")

tsneembedding_brain <- Embeddings(object = brain, reduction = "tsne")
brain[["annotation"]] <- Idents(object = brain)

# Separate nuclei and cells and only retain neuronal clusters
brain <- SetIdent(brain, value = brain@meta.data$seurat_clusters)
onlycells <- SetIdent(onlycells, value = onlycells@meta.data$seurat_clusters)
onlynuclei <- SetIdent(onlynuclei, value = onlynuclei@meta.data$seurat_clusters)
cellsneurons <- subset(onlycells, idents = c('4','10','54','57','70','73','82','80','45','67'), invert = TRUE)
nucleineurons <- subset(onlynuclei, idents = c('4','10','54','57','70','73','82','45'), invert = TRUE)
cellsraw <- GetAssayData(object = cellsneurons, slot = "counts")
nucleiraw <- GetAssayData(object = nucleineurons, slot = "counts")
cellsraw_matrix <- as.matrix(cellsraw)
nucleiraw_matrix <- as.matrix(nucleiraw)

#read in data matrix nuclei
datasets <- list(nucleiraw_matrix, cellsraw_matrix)
results <- iCytoTRACE(datasets)
plotCytoTRACE(results, emb = tsneembedding_brain)

plotCytoGenes(results, numOfGenes = 10)





