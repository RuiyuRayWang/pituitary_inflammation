#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(Matrix)

source("rscripts/utilities.R")
source("rscripts/seurat_to_anndata_f.R")

SeuratToAnndata(f_name = "cells_postprocessed", f_path = "data/processed/", out_dir = "data/processed/")
