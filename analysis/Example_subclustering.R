################################################################
# Mould Recluster Macrophages
# This is an example of subsetting a seurat object to cell of 
# Inteerest and re-clustering
#################################################################
# This is Kara's data which is located here: /Moore/home/moorec/mould/seurat/joint_integration_new/

# Load Libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(openxlsx)
library(ggplot2)
library(here)



options(future.globals.maxSize = 100*1000 * 1024^2)

# Read in data
integrated_hume <- readRDS( here('data','integrated_all_samples.RDS'))

# Cell type info provided by Kara
# integrated_hume@meta.data$overall_celltype <- ifelse(integrated_hume@meta.data$seurat_cluster ==17, 'B cell',
#                                                  ifelse(integrated_hume@meta.data$seurat_cluster%in% c(11, 13), 'Cycling',
#                                                         ifelse(integrated_hume@meta.data$seurat_cluster%in% c(12, 15), 'DC',
#                                                                ifelse(integrated_hume@meta.data$seurat_cluster ==14, 'Epi and Mast',
#                                                                       ifelse(integrated_hume@meta.data$seurat_cluster ==18, 'Migratory DC',
#                                                                              ifelse(integrated_hume@meta.data$seurat_cluster ==16, 'pDC',
#                                                                                     ifelse(integrated_hume@meta.data$seurat_cluster%in% c(0,2,3,7,9,10), 'RAM',
#                                                                                            ifelse(integrated_hume@meta.data$seurat_cluster%in% c(1,5,8), 'RecAM',
#                                                                                                   'T cell'
#                                                         ))))))))

#17	B cells
#11	Cycling1
#13	Cycling2
#15	DC1
#12	DC2
#14	Epi and Mast
#18	Migratory DC
#16	pDC
#0	RAM
#2	RAM
#3	RAM
#7	RAM
#9	RAM
#10	RAM
#1	RecAM
#5	RecAM
#8	RecAM
#4	T cells
#6	T cells

# Check that it was coded right
# table(integrated_hume@meta.data$seurat_cluster, integrated_hume@meta.data$overall_celltype)

# Question: have t cell and b cell clusters but would need the rest to be able to do all clusters 
# Save cell type info
# saveRDS(integrated_hume, 'integrated_hume_all_samples.RDS')



include <- c(1,10)

# Subset data

recluster <- integrated_hume[,integrated_hume@meta.data$seurat_clusters %in% include]


# Save overall cluster as reclusting will write over the seurat_clusters variable
recluster@meta.data$overall_cluster <- recluster@meta.data$seurat_clusters


# remove(integrated_hume)

# Choose number of PCAs to use in clustering - we used 25 for our original clustering
# so I stuck with this setting
nPCAs = 25
res = 0.5

DefaultAssay(recluster) <- "integrated"

recluster <- FindNeighbors(object = recluster, dims = 1:nPCAs)

recluster <- RunUMAP(object = recluster, reduction = "pca", dims = 1:nPCAs, n.neighbors=50)

recluster <- FindClusters(recluster, reduction="pca", resolution=res, algorithm=3) # 

##### UMAP Plots ######

pdf(here('Tcell_recluster_UMAP.pdf'))

p1 = DimPlot(recluster, reduction='umap', label=T)
print(p1)

p1 <- DimPlot(object = recluster, reduction = "umap", group.by = "LibraryID")
print(p1)

p1 <- DimPlot(object = recluster, reduction = "umap", group.by = "tissue")
print(p1)

# p1 <- DimPlot(object = recluster, reduction = "umap", group.by = "subject")
# print(p1)

p1 = FeaturePlot(recluster, c("nCount_RNA"), cols=c("gray", "blue"), pt.size=1, min.cutoff="q5", max.cutoff="q95")
print(p1)
p1 = FeaturePlot(recluster, c("nFeature_RNA"), cols=c("gray", "blue"), pt.size=1, min.cutoff="q5", max.cutoff="q95")
print(p1)

p1 = FeaturePlot(recluster, c("percent.mito"), cols=c("gray", "blue"), pt.size=1, min.cutoff="q5", max.cutoff="q95")
print(p1)

p1 = FeaturePlot(recluster, c("percent.ribo"), cols=c("gray", "blue"), pt.size=1, min.cutoff="q5", max.cutoff="q95")
print(p1)

dev.off()


##### QC Plots by clusters

# MAKE QC PLOTS BY CLUSTER
pdf(paste0("tcell","_cluster_QC_plots.pdf"))
VlnPlot(recluster, c("nCount_RNA"), pt.size=0)
VlnPlot(recluster, c("nFeature_RNA"), pt.size=0)
VlnPlot(recluster, c("percent.mito"), pt.size=0)
VlnPlot(recluster, c("percent.ribo"), pt.size=0)
# VlnPlot(recluster, c("bcds_score"), pt.size=0)
# VlnPlot(recluster, c("cxds_score"), pt.size=0)
# VlnPlot(recluster, c("hybrid_score"), pt.size=0)
dev.off()




# # Follow Seurat clustering steps - make sure you are using the integrated_hume assay
# DefaultAssay(recluster) <- "integrated"
# 
# pca_list = c(20, 25, 30)
# 
# for (nPCAs in pca_list){
# 
# recluster <- FindNeighbors(object = recluster, dims = 1:nPCAs)
# 
# recluster <- RunUMAP(object = recluster, reduction = "pca", dims = 1:nPCAs, n.neighbors=50)
# 
# recluster <- FindClusters(recluster, reduction="pca", resolution=0.4, algorithm=3)
# 
# # Make a plot / UMAP to display results
# # Here we color by cluster, sample, subject, time, some QC stats
# pdf(here('results',paste0("PCA_",as.character(nPCAs),"Tcell_recluster_UMAP_sctrans.pdf")))
# p1 = DimPlot(recluster, reduction='umap', label=T)
# print(p1)
# 
# p1 <- DimPlot(object = recluster, reduction = "umap", group.by = "LibraryID")
# print(p1)
# 
# p1 <- DimPlot(object = recluster, reduction = "umap", group.by = "tissue")
# print(p1)
# 
# # p1 <- DimPlot(object = recluster, reduction = "umap", group.by = "subject")
# # print(p1)
# 
# p1 = FeaturePlot(recluster, c("nCount_RNA"), cols=c("gray", "blue"), pt.size=1, min.cutoff="q5", max.cutoff="q95")
# print(p1)
# p1 = FeaturePlot(recluster, c("nFeature_RNA"), cols=c("gray", "blue"), pt.size=1, min.cutoff="q5", max.cutoff="q95")
# print(p1)
# 
# p1 = FeaturePlot(recluster, c("percent.mito"), cols=c("gray", "blue"), pt.size=1, min.cutoff="q5", max.cutoff="q95")
# print(p1)
# 
# p1 = FeaturePlot(recluster, c("percent.ribo"), cols=c("gray", "blue"), pt.size=1, min.cutoff="q5", max.cutoff="q95")
# print(p1)
# 
# dev.off()
# 
# }
# 
# # ## after looking at results with different pcas, 25 seems reasonable
# 
# resolution_list = c(.3, .4, .5)
# 
# nPCAs = 25
# for (res in resolution_list){
# 
# 
#   # Follow Seurat clustering steps - make sure you are using the integrated_hume assay
#   DefaultAssay(recluster) <- "integrated"
# 
#   recluster <- FindNeighbors(object = recluster, dims = 1:nPCAs)
# 
#   recluster <- RunUMAP(object = recluster, reduction = "pca", dims = 1:nPCAs, n.neighbors=50)
# 
#   recluster <- FindClusters(recluster, reduction="pca", resolution=res, algorithm=3)
# 
#   # Make a plot / UMAP to display results
#   # Here we color by cluster, sample, subject, time, some QC stats
#   pdf(here('results',paste0("RES_",as.character(res),"_Tcell_recluster_UMAP_sctrans.pdf")))
#   p1 = DimPlot(recluster, reduction='umap', label=T)
#   print(p1)
# 
#   p1 <- DimPlot(object = recluster, reduction = "umap", group.by = "LibraryID")
#   print(p1)
# 
#   p1 <- DimPlot(object = recluster, reduction = "umap", group.by = "tissue")
#   print(p1)
# 
#   # p1 <- DimPlot(object = recluster, reduction = "umap", group.by = "subject")
#   # print(p1)
# 
#   p1 = FeaturePlot(recluster, c("nCount_RNA"), cols=c("gray", "blue"), pt.size=1, min.cutoff="q5", max.cutoff="q95")
#   print(p1)
#   p1 = FeaturePlot(recluster, c("nFeature_RNA"), cols=c("gray", "blue"), pt.size=1, min.cutoff="q5", max.cutoff="q95")
#   print(p1)
# 
#   p1 = FeaturePlot(recluster, c("percent.mito"), cols=c("gray", "blue"), pt.size=1, min.cutoff="q5", max.cutoff="q95")
#   print(p1)
# 
#   p1 = FeaturePlot(recluster, c("percent.ribo"), cols=c("gray", "blue"), pt.size=1, min.cutoff="q5", max.cutoff="q95")
#   print(p1)
# 
#   dev.off()
# 
# }


# 
# nPCAs = 25
# res = 0.4
# 
# DefaultAssay(recluster) <- "integrated"
# 
# recluster <- FindNeighbors(object = recluster, dims = 1:nPCAs)
# 
# recluster <- RunUMAP(object = recluster, reduction = "pca", dims = 1:nPCAs, n.neighbors=50)
# 
# recluster <- FindClusters(recluster, reduction="pca", resolution=res, algorithm=3) # 
DimPlot(recluster, reduction='umap', label=T)

# #========================================#
# # Freq and Proportion of cell type in each sample #
# #========================================#
write.xlsx(list(cluster_freq = table(recluster@meta.data$seurat_clusters),
                cluster_pct = 100*prop.table(table(recluster@meta.data$seurat_clusters)),
                overall_cluster_freq = table(recluster@meta.data$overall_cluster),
                overall_cluster_pct = 100*prop.table(table(recluster@meta.data$overall_cluster)),
                sample_freq = table(recluster@meta.data$seurat_clusters, recluster@meta.data$LibraryID),
                sample_pct = 100*prop.table(table(recluster@meta.data$seurat_clusters, recluster@meta.data$LibraryID),2),
                tissue_freq = table(recluster@meta.data$seurat_clusters, recluster@meta.data$tissue),
                tissue_pct = 100*prop.table(table(recluster@meta.data$seurat_clusters, recluster@meta.data$tissue),2)),
           here('results','joint_sample_Tcell_recluster_freq_pct.xlsx'))

# write.xlsx(list(as.data.frame(table(recluster@meta.data$seurat_clusters)) %>% rename(cluster = Var1, N = Freq),
#                 cluster_pct = as.data.frame(100*prop.table(table(recluster@meta.data$seurat_clusters))) %>% rename(cluster = Var1, Percent = Freq),
#                 overall_cluster_freq = as.data.frame(table(recluster@meta.data$overall_cluster)) %>% rename(cluster = Var1, N = Freq),
#                 overall_cluster_pct = as.data.frame(100*prop.table(table(recluster@meta.data$overall_cluster))) %>% rename(cluster = Var1, Percent = Freq),
#                 sample_freq = table(recluster@meta.data$seurat_clusters, recluster@meta.data$LibraryID),
#                 sample_pct = 100*prop.table(table(recluster@meta.data$seurat_clusters, recluster@meta.data$LibraryID),2),
#                 tissue_freq = table(recluster@meta.data$seurat_clusters, recluster@meta.data$tissue),
#                 tissue_pct = 100*prop.table(table(recluster@meta.data$seurat_clusters, recluster@meta.data$tissue),2)),
#            here('results','joint_sample_Tcell_recluster_freq_pct.xlsx'))



# Save the dataset
saveRDS(recluster, here('data','Tcell_recluster.RDS'))



#########################################################################
# Create CLOUPE file for Kara to work with
#########################################################################

md <- recluster@meta.data

md$UMAP_1 <- recluster@reductions$umap@cell.embeddings[,1]
md$UMAP_2 <- recluster@reductions$umap@cell.embeddings[,2]


# x <-strsplit(rownames(md), '_')
# md$Barcode <- unlist(lapply(x, function(x) x[1]))

write.csv(md[,c('Barcode', 'UMAP_1', 'UMAP_2')], 'tcell_recluster_umap_coordinates.csv', row.names = FALSE)

md$Tcell_recluster <- md$seurat_clusters

write.csv(md[,c('Barcode', 'Tcell_recluster')], 'tcell_recluster_cloupe_categories.csv', row.names = FALSE)

new_clusters = left_join(md[,c('Barcode', 'Tcell_recluster')], md[,c('Barcode', 'UMAP_1', 'UMAP_2')]) %>% select(Barcode, UMAP_1, UMAP_2, Tcell_recluster)
write.csv(new_clusters, 'tcell_recluster_cloupe_all.csv', row.names = FALSE)

##################################################
# Marker finding for each cluster
##################################################
###Find SCT markers

DefaultAssay(recluster) <- "SCT"

nClusters <- length(unique(recluster@meta.data$seurat_clusters))

# Marker finding ignoring clustering of cells in samples
markers.list <- list()

recluster2 = PrepSCTFindMarkers(recluster, assay = "SCT") # Question the FindMarkers function required I do this first; not sure if right
#Given a merged object with multiple SCT models, this function uses minimum of the median UMI (calculated using the raw UMI counts) of 
#individual objects to reverse the individual SCT regression model using minimum of median UMI as the sequencing depth covariate. 
#The counts slot of the SCT assay is replaced with recorrected counts and the data slot is replaced with log1p of recorrected counts.

for(i in 1:nClusters) {
  temp <- FindMarkers(recluster2, ident.1=(i-1), min.pct=0.1, logfc.threshold=0.25,  only.pos=T, assay = 'SCT') 
  markers.list[[i]] <- temp[temp$p_val_adj < 0.05,]
}

names(markers.list) <- paste0('cl_', seq(0, nClusters-1,1))
write.xlsx(markers.list, here('results','tcell_recluster_sct_combined_markers_list.xlsx'), rowNames=T)

combined <- markers.list


# Get Sample Specific Info for each gene in the combined list
# Meaning repeat the analysis per sample
# This is the annoyng code I mentioned, where we check to make sure there are enough cells per sample
markers.list <- list()

sample_freq = table(recluster2@meta.data$seurat_clusters, recluster2@meta.data$orig.ident)
orig.idents <- as.character(colnames(sample_freq))

for (i in 1:nClusters){
  low_count <- sample_freq[i,]<3
  if(sum(low_count)==0){
    temp <- FindConservedMarkers(recluster2, ident.1=i-1, min.pct=0, logfc.threshold=0,  
                                 only.pos=F, assay = 'SCT', features = rownames(combined[[i]]),
                                 grouping.var = 'orig.ident')
    
  }else{
    temp <- FindConservedMarkers(recluster2[,!(recluster2@meta.data$orig.ident %in% orig.idents[low_count])], ident.1=i-1, min.pct=0, logfc.threshold=0,  
                                 only.pos=F, assay = 'SCT', features = rownames(combined[[i]]),
                                 grouping.var = 'orig.ident')
  }
  
  markers.list[[i]] <- temp
  
  if(dim(markers.list[[i]])[1]!=0){
    markers.list[[i]]$n_samples_significant <- rowSums(markers.list[[i]][,grep("p_val_adj", colnames(markers.list[[i]]))]<0.05)
    
  }}

names(markers.list) <- paste0('cl_', seq(0,nClusters-1,1))

write.xlsx(markers.list, here('results','tcell_recluster_sct_sample_specific_markers_list.xlsx'), rowNames=T)


## Quick doublet check based on binary classification 

bcds_score <- read_csv("~/Desktop/bcds_score.csv")

bcds = left_join(recluster@meta.data, bcds_score)

bcds = bcds %>% mutate(bcds_bin = case_when(
  bcds_score == "<= 0.5" ~ 0,
  TRUE ~ 1
))

library(scales)

pdf(here("results", "doublet_binary_plot.pdf"))
ggplot(bcds, aes(seurat_clusters, bcds_bin)) +
  stat_summary(fun=mean, geom="bar", fill="grey70") +
  stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=0.2) +
  scale_y_continuous(labels=percent_format(), limits=c(0,1)) +
  theme_bw()
dev.off()
