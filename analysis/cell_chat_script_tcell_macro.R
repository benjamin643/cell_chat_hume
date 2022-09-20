#!/usr/bin/env -S Rscript --vanilla
#
# Set the log file and job name (%j is replaced with job id)
#SBATCH -o /Moore/home/steinhartb/practice/logs/practice.%j.log
#SBATCH -e /Moore/home/steinhartb/practice/errs/practice.%j.log

#SBATCH -J practice
#
# Set RAM, nodes, and cores per node
#SBATCH --mem 700G
#SBATCH -n 1
#SBATCH -c 96

library(CellChat)
library(patchwork)
library(tidyverse)
library(here)
library(Seurat)
library(SeuratObject)
library(NMF)
library(ggalluvial)
options(stringsAsFactors = FALSE, warn=-1, message = -1)

# Part I: Data input & processing and initialization of CellChat object

integrated_tcell <- readRDS('/Moore/data/hume/hume_dat_12_2021/Tcell_recluster.RDS')
integrated_macro <- readRDS('/Moore/data/hume/hume_dat_12_2021/integrated_M4.RDS')

n_macro_clusters = integrated_macro@meta.data$seurat_clusters %>% unique() %>% length()
n_tcell_clusters = integrated_tcell@meta.data$seurat_clusters %>% unique() %>% length()

## convert from from factor to numeric 

integrated_tcell@meta.data$seurat_clusters = as.numeric(integrated_tcell@meta.data$seurat_clusters)
integrated_macro@meta.data$seurat_clusters = as.numeric(integrated_macro@meta.data$seurat_clusters)

integrated_tcell@meta.data$seurat_clusters = integrated_tcell@meta.data$seurat_clusters + n_macro_clusters #clusters now 7:12
integrated_macro@meta.data$seurat_clusters = integrated_macro@meta.data$seurat_clusters # clusters now 1:6

pbmc.combined <- merge(integrated_tcell, y = integrated_macro, add.cell.ids = c("tcell", "macro"), project = "PBMC12K")
pbmc.combined

head(colnames(pbmc.combined))


data.input = GetAssayData(pbmc.combined@assays$SCT)# normalized data matrix ; Question: I think this is our normalized data via SCT transform
meta = pbmc.combined@meta.data # a dataframe with rownames containing cell mata data
# cell.use = rownames(meta)[meta$condition == "LS"] # extract the cell names from disease data 
# Question: integrated_snn)res.0.8 is all NAs


meta$labels = factor(meta$seurat_clusters, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), labels = c("mac:1", "mac:2","mac:3","mac:4","mac:5","mac:6","tcell:1","tcell:2","tcell:3","tcell:4","tcell:5","tcell:6"))


meta = meta %>% mutate(cluster_labels = case_when(
  seurat_clusters == "1" ~ "mac:1",
  seurat_clusters == "2" ~ "mac:2",
  seurat_clusters == "3" ~ "mac:3",
  seurat_clusters == "4" ~ "mac:4",
  seurat_clusters == "5" ~ "mac:5",
  seurat_clusters == "6" ~ "mac:6",
  seurat_clusters == "7" ~ "tcell:1",
  seurat_clusters == "8" ~ "tcell:2",
  seurat_clusters == "9" ~ "tcell:3",
  seurat_clusters == "10" ~ "tcell:4",
  seurat_clusters == "11" ~ "tcell:5",
  seurat_clusters == "12" ~ "tcell:6",
))


## Create a CellChat object

samples = meta$orig.ident %>% unique()
# look at interaction between macrophage subgroups and 1 tcell 

for(i in 1:length(samples)){
  if (i>1){
    
    integrated_tcell <- readRDS('/Moore/data/hume/hume_dat_12_2021/Tcell_recluster.RDS')
    integrated_macro <- readRDS('/Moore/data/hume/hume_dat_12_2021/integrated_M4.RDS')
    
    n_macro_clusters = integrated_macro@meta.data$seurat_clusters %>% unique() %>% length()
    n_tcell_clusters = integrated_tcell@meta.data$seurat_clusters %>% unique() %>% length()
    
    ## convert from from factor to numeric 
    
    integrated_tcell@meta.data$seurat_clusters = as.numeric(integrated_tcell@meta.data$seurat_clusters)
    integrated_macro@meta.data$seurat_clusters = as.numeric(integrated_macro@meta.data$seurat_clusters)
    
    integrated_tcell@meta.data$seurat_clusters = integrated_tcell@meta.data$seurat_clusters + n_macro_clusters #clusters now 7:12
    integrated_macro@meta.data$seurat_clusters = integrated_macro@meta.data$seurat_clusters # clusters now 1:6
    
    pbmc.combined <- merge(integrated_tcell, y = integrated_macro, add.cell.ids = c("tcell", "macro"), project = "PBMC12K")
    pbmc.combined
    
    head(colnames(pbmc.combined))
    
    
    data.input = GetAssayData(pbmc.combined@assays$SCT)# normalized data matrix ; Question: I think this is our normalized data via SCT transform
    meta = pbmc.combined@meta.data # a dataframe with rownames containing cell mata data
    # cell.use = rownames(meta)[meta$condition == "LS"] # extract the cell names from disease data 
    # Question: integrated_snn)res.0.8 is all NAs
    
    
    meta$labels = factor(meta$seurat_clusters, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), labels = c("mac:1", "mac:2","mac:3","mac:4","mac:5","mac:6","tcell:1","tcell:2","tcell:3","tcell:4","tcell:5","tcell:6"))
    
    
    meta = meta %>% mutate(cluster_labels = case_when(
      seurat_clusters == "1" ~ "mac:1",
      seurat_clusters == "2" ~ "mac:2",
      seurat_clusters == "3" ~ "mac:3",
      seurat_clusters == "4" ~ "mac:4",
      seurat_clusters == "5" ~ "mac:5",
      seurat_clusters == "6" ~ "mac:6",
      seurat_clusters == "7" ~ "tcell:1",
      seurat_clusters == "8" ~ "tcell:2",
      seurat_clusters == "9" ~ "tcell:3",
      seurat_clusters == "10" ~ "tcell:4",
      seurat_clusters == "11" ~ "tcell:5",
      seurat_clusters == "12" ~ "tcell:6",
    ))
    
  }
 

cell.use = rownames(meta)[meta$orig.ident ==  samples[i]] # extract the cell names from disease data

# mainDir = "/Users/steinhartb/Desktop/cell_chat_hume/results"
mainDir = "/Moore/home/steinhartb/cell-chat/results"

subDir = samples[i] 
dir.create(file.path(mainDir, subDir))
data.input = data.input[, cell.use]
meta = meta[cell.use, ]

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cluster_labels") # takes one group out : Question 

## Set the ligand-receptor interaction database

# CellChatDB in human contains 1,939 validated molecular interactions, including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions and 16.5% of cell-cell contact interactions.

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use

## Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Part II: Inference of cell-cell communication network

## Compute the communication probability and infer cellular communication network

cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10) # Question; arbitrary here? Could change? 

## Extract the inferred cellular communication network as a data frame

## Infer the cell-cell communication at a signaling pathway level

cellchat <- computeCommunProbPathway(cellchat)

## Calculate the aggregated cell-cell communication network 

cellchat <- aggregateNet(cellchat)

# Part III: Visualization of cell-cell communication network

## Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram

cellchat@netP$pathways # All the signaling pathways showing significant communications 

signal_pathways =cellchat@netP$pathways

try(
for (signal in signal_pathways){
  pathways.show <- signal 
  pdf(paste0('/Moore/home/steinhartb/cell-chat/results/',subDir ,'/signal_pathway_', signal,'.pdf'))
  
  # Circle plot
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  # Chord diagram
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  # Heatmap
  # netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds") # question ; heatmap not working 
  dev.off()
}, silent = TRUE)

### Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair ; choosing to show top 3 pathways 

pathways.show.all <- cellchat@netP$pathways

try(
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show.all, geneLR.return = FALSE), silent = TRUE)

pdf(paste0('/Moore/home/steinhartb/cell-chat/results/',subDir ,'/signal_ligand_receptor_pairs.pdf'))

try(
for (i in 1:nrow(pairLR.CXCL)){
  
  LR.show <- pairLR.CXCL[i,] # show one ligand-receptor pair
  # Circle plot
  netVisual_individual(cellchat, signaling = pathways.show.all, pairLR.use = LR.show, layout = "circle")
  # Chord diagram
  netVisual_individual(cellchat, signaling = pathways.show.all, pairLR.use = LR.show, layout = "chord")
}, silent = TRUE)

dev.off()


### Automatically save the plots of the all inferred network for quick exploration

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
pdf(paste0('/Moore/home/steinhartb/cell-chat/results/',subDir ,'/L-R_contributions.pdf'))
try(
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  # netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy") # question; these done get added to proper results folder
  gg = netAnalysis_contribution(cellchat, signaling = pathways.show.all[i]) 
  gg = gg + ggtitle(paste0("Contribution of each L-R Pair: ",pathways.show.all[i])) +
    theme(plot.title = element_text(hjust = 0.5))
  print(gg)
}, silent = TRUE)

dev.off()

# The dot color and size represent the calculated communication probability and p-values. p-values are computed from one-sided permutation test.

## Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
### Bubble plot

pdf(paste0('/Moore/home/steinhartb/cell-chat/results/',subDir ,'/cell_com_mult_receptors_or_pathways.pdf'))

try(
for (i in 1:n_macro_clusters){
  p1 = netVisual_bubble(cellchat, sources.use = i, targets.use = c(7:12), remove.isolate = FALSE, title.name = paste0("Cell Communication mediated by multiple pathways: macro", i))
  print(p1)
}, silent = TRUE)
dev.off()



### Chord diagram

# netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(7:12), lab.cex = 0.5,legend.pos.y = 30)

# pdf(here(paste0('results/',subDir),paste0('all_LR_interactions_from_macro_cluster.pdf')))

pdf(paste0('/Moore/home/steinhartb/cell-chat/results/',subDir ,'/all_LR_interactions_from_macro_cluster.pdf'))

try(
for (i in 1:n_macro_clusters){
  p1 = netVisual_bubble(cellchat, sources.use = i, targets.use = c(7:12), remove.isolate = FALSE, title.name = paste0("Cell Communication mediated by multiple pathways:macro", i))
  print(p1)
}, silent = TRUE)

dev.off()

## Plot the signaling gene expression distribution using violin/dot plot

pdf(paste0('/Moore/home/steinhartb/cell-chat/results/',subDir ,'/gene_expression_signaling_genes_to_LR_pairs.pdf'))

try(
for (path in pathways.show.all){
  p1 = plotGeneExpression(cellchat, signaling = path)  + plot_annotation(title = paste0("Signaling Pathway: ", path) ,theme = theme(plot.title = element_text(hjust = 0.5)))
  print(p1)
}, silent = TRUE)
dev.off()

# Part IV: Systems analysis of cell-cell communication network

### Compute and visualize the network centrality scores

# Compute the network centrality scores


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways # Question; weird error message that makes me nervous. 

try(
suppressWarnings({ pdf(paste0('/Moore/home/steinhartb/cell-chat/results/',subDir ,'/network_centrality_scores.pdf'))
  
  p1 = netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all, width = 8, height = 2.5, font.size = 10)
  print(p1)
  dev.off() }), silent = TRUE)



### Identify signals contributing most to outgoing or incoming signaling of certain cell groups

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways

pdf(paste0('/Moore/home/steinhartb/cell-chat/results/',subDir ,'/outgoing_and_incoming_signaling_strength.pdf'))

try(
  suppressWarnings({ 
h1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")

prac1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

print(ht1)

print(ht2)

dev.off()}), silent = TRUE)


### Identify and visualize outgoing communication pattern of secreting cells
try(
suppressWarnings({pdf(paste0('/Moore/home/steinhartb/cell-chat/results/',subDir ,'/outgoing_pattern_analysis.pdf'))
  
  outs = selectK(cellchat, pattern = "outgoing")
  
  print(outs)
  
  dev.off()
}), silent = TRUE)

try(
suppressWarnings({pdf(paste0('/Moore/home/steinhartb/cell-chat/results/',subDir ,'/incoming_pattern_analysis.pdf'))
  
  p3 = selectK(cellchat, pattern = "incoming")
  print(p3)
  
  dev.off()
}), silent = TRUE)


}
