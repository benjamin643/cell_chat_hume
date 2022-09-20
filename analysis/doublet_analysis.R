#!/usr/bin/env -S Rscript --vanilla
#
# Set the log file and job name (%j is replaced with job id)
#SBATCH -o /Moore/home/steinhartb/HuLPS_PBMC/logs/091922.%j.log
#SBATCH -e /Moore/home/steinhartb/HuLPS_PBMC/errs/091922.%j.log

#SBATCH -J 091922_QC
#
# Set RAM, nodes, and cores per node
#SBATCH --mem 70G
#SBATCH -n 1
#SBATCH -c 96



#####################################
# PRIMERO SC Data 
# READ DATA INTO SEURAT
# PLOT QUALITY CONTROL METRICS TO MAKE
# DECISIONS ON FILTERING STRATEGIES
#####################################

# Load Libraries
# Be sure to install all of these
library(Seurat)
library(dplyr)
library(Matrix)
library(openxlsx)
library(ggplot2)
# library(scater)
# library(loomR)
library(scds)
library(scran)


options(future.globals.maxSize= 800000000000)

#### sample name (for naming files)
project_name <- "PBMC_batches"


#####################################
# READ IN DATA
#####################################


pbmc1 <- "/Moore/data/mould/PBMC_batch1/"
sample_names_gex = list.files(pbmc1, pattern =  "_GEX$")


pbmc2 <- "/Moore/data/mould/PBMC_batch2/"
sample_names_pbmc = list.files(pbmc2, pattern =  "_PBMC$")

print(sample_names_pbmc)

sample_names = c(sample_names_pbmc,sample_names_gex)
print(sample_names)

print(length(sample_names))
# 
# 
# ####  create a list where you read in the 10X data for each sample
# 
# my_dir = "/Volumes/Mould/data/HuLPS_PBMC"
# dat <- list()
# for(sample_name in sample_names) {
#   data_dir <- paste0(my_dir, "/", sample_name, "/filtered_feature_bc_matrix")
#   dat[[sample_name]] <- Read10X(data.dir = data_dir)
# }

my_dir = "/Moore/home/steinhartb/HuLPS_PBMC"
my_dir_gex = "/Moore/data/mould/PBMC_batch1"
my_dir_pbmc = "/Moore/data/mould/PBMC_batch2"



dat <- list()
for(sample_name in sample_names_gex) {
  data_dir <- paste0(my_dir_gex, "/", sample_name, "/outs/filtered_feature_bc_matrix")
  dat[[sample_name]] <- Read10X(data.dir = data_dir)
}

for(sample_name in sample_names_pbmc) {
  data_dir <- paste0(my_dir_pbmc, "/", sample_name, "/outs/filtered_feature_bc_matrix")
  dat[[sample_name]] <- Read10X(data.dir = data_dir)
}



# #####################################
# # CALCULATE QC METRICS
# #####################################
# 
# #### READ DATA INTO SEURAT
# # HERE WE REQUIRE GENES TO HAVE >1 COUNT IN AT LEAST 3 CELLS
# # AND CELLS TO HAVE AT LEAST 100 FEATURES OR ELSE THEY ARE FILTERED OUT
# # YOU CAN ADJUST THESE, BUT THIS IS A PRETTY LIBERAL PLACE TO START
# 
sDat <- list()
for(sample_name in sample_names) {
  sDat[[sample_name]] <- CreateSeuratObject(counts = dat[[sample_name]], min.cells = 3, min.features = 100, project = sample_name)
}

# #### CALCULATE SOME QC STATISTICS
# # % MITOCHONDRIAL READS AND % RIBOSOMAL READS
# # DOUBLET DETECTION METRICS

for(sample_name in names(sDat)) {
  mito.genes <- grep(pattern = "^MT-|^MRPS|^MRPL", x = rownames(x = sDat[[sample_name]]), value = TRUE)
  percent.mito <- Matrix::colSums(x=GetAssayData(sDat[[sample_name]], slot='counts')[mito.genes,])/Matrix::colSums(x=GetAssayData(sDat[[sample_name]], slot='counts'))
  
  ribo.genes <- grep(pattern = "^RPS|^RPL", x = rownames(x = sDat[[sample_name]]), value = TRUE)
  percent.ribo <- Matrix::colSums(x=GetAssayData(sDat[[sample_name]], slot='counts')[ribo.genes,])/Matrix::colSums(x=GetAssayData(sDat[[sample_name]], slot='counts'))
  
  sDat[[sample_name]][['percent.mito']] <- percent.mito
  sDat[[sample_name]][['percent.ribo']] <- percent.ribo
  
  # DOUBLET SCORING
  sce <-SingleCellExperiment(list(counts=GetAssayData(sDat[[sample_name]],
                                                      assay="RNA",
                                                      slot="counts")))
  
  # Annotate doublet using co-expression based doublet scoring:
  sce = cxds(sce)
  
  # Annotate doublet using binary classification based doublet scoring:
  sce = bcds(sce)
  
  # Combine both annotations into a hybrid annotation
  sce = cxds_bcds_hybrid(sce)
  
  # Doublet scores are now available via colData:
  sDat[[sample_name]]@meta.data$cxds_score <- colData(sce)$cxds_score
  sDat[[sample_name]]@meta.data$bcds_score <- colData(sce)$bcds_score
  sDat[[sample_name]]@meta.data$hybrid_score <- colData(sce)$hybrid_score
  
}

saveRDS(sDat, paste0(my_dir, "/data/seurat_data.RDS"))

# merge data preQC meta data into a single data frame
sDat_merge <- merge(sDat[[names(sDat)[1]]], y = sDat[names(sDat)[-1]], add.cell.ids = names(sDat), project = project_name)
sDat_merge@meta.data$sampleID <- sDat_merge@meta.data$orig.ident

# Make some QC plots by sample
pdf(paste0(my_dir,"results/all_samples_pre_QC_VlnPlot_no_dots.pdf"), width=11, height=8.5)
VlnPlot(object=sDat_merge, features="nCount_RNA", group.by="sampleID", pt.size=0)+ theme(legend.position="none")
VlnPlot(object=sDat_merge, features="nCount_RNA", group.by="sampleID", log=T, pt.size=0)+ theme(legend.position="none")
VlnPlot(object=sDat_merge, features="nFeature_RNA", group.by="sampleID", pt.size=0)+ theme(legend.position="none")
VlnPlot(object=sDat_merge, features="nFeature_RNA", group.by="sampleID", log=T, pt.size=0)+ theme(legend.position="none")
VlnPlot(object=sDat_merge, features="percent.mito", group.by="sampleID", pt.size=0)+ theme(legend.position="none")
VlnPlot(object=sDat_merge, features="percent.ribo", group.by="sampleID", pt.size=0)+ theme(legend.position="none")
VlnPlot(object=sDat_merge, features="cxds_score", group.by="sampleID", pt.size=0)+ theme(legend.position="none")
VlnPlot(object=sDat_merge, features="bcds_score", group.by="sampleID", pt.size=0)+ theme(legend.position="none")
VlnPlot(object=sDat_merge, features="hybrid_score", group.by="sampleID", pt.size=0)+ theme(legend.position="none")
dev.off()

