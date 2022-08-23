#install the following packages
library("BiocManager")
library("multtest")
library(rlang)
library(Rcpp)
library(Seurat)
library(patchwork)
library(dplyr)
library(sctransform)
###### Pre-processing steps for cells ###### 

###Read in data
setwd("~/first draft/cellranger output/cells")
data = Read10X(data.dir = "filtered_feature_bc_matrix", unique.features = T)

###Filter and normalize dataset using scTransform
cells <- CreateSeuratObject(counts = data, project = "cells", min.cells = 10)
cells[["percent.mt"]] <- PercentageFeatureSet(cells, pattern = "^mt-")

#first subset
cells <- subset(cells, subset = nFeature_RNA > 800 & nFeature_RNA < 4000 & percent.mt < 5)
#then normalize
cells <- SCTransform(cells, variable.features.n = 2000, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = TRUE)

#run PCA
DefaultAssay(cells) <- "SCT"
cells <- RunPCA(cells, verbose = FALSE, npcs = 150)
cells <- RunUMAP(cells, dims = 1:150, verbose = FALSE)
cells <- RunTSNE(cells, dims = 1:150, verbose = FALSE)

###### Pre-processing steps for nuclei ###### 
###Read in data
setwd("~/first draft/cellranger output/nuclei")
data = Read10X(data.dir = "filtered_feature_bc_matrix", unique.features = T)

###Filter and normalize dataset using scTransform
nuclei <- CreateSeuratObject(counts = data, project = "nuclei", min.cells = 10)
nuclei[["percent.mt"]] <- PercentageFeatureSet(nuclei, pattern = "^mt-")

#first subset
nuclei <- subset(nuclei, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 5)
#then normalize
nuclei <- SCTransform(nuclei, variable.features.n = 2000, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = TRUE)

#run PCA
DefaultAssay(nuclei) <- "SCT"
nuclei <- RunPCA(nuclei, verbose = FALSE, npcs = 150)
nuclei <- RunUMAP(nuclei, dims = 1:150, verbose = FALSE)
nuclei <- RunTSNE(nuclei, dims = 1:150, verbose = FALSE)

saveRDS(nuclei, file = "../nucleisseuratv3.rds")
saveRDS(cells, file = "../cellsseuratv3.rds")

###### Data integration and pre-processing steps ###### 
#Read in data
nuclei <- readRDS("nucleisseuratv3.rds")
cells <- readRDS("cellsseuratv3.rds")
#merge datasets
brain  <- merge(cells, y = nuclei, add.cell.ids = c("C","N"), project = "brain")
table(brain$orig.ident)
#merging and filtering is done separately for each dataset
brain.list <- SplitObject(brain, split.by = "orig.ident")
brain.list <- brain.list[c("cells","nuclei")]
for (i in 1:length(brain.list)) {
  brain.list[[i]] <- SCTransform(brain.list[[i]], verbose = FALSE)
}
brain.features <- SelectIntegrationFeatures(object.list = brain.list, nfeatures = 3000)
brain.list <- PrepSCTIntegration(object.list = brain.list, anchor.features = brain.features)
brain.anchors <- FindIntegrationAnchors(object.list = brain.list, normalization.method = "SCT", anchor.features = brain.features)
brain.integrated <- IntegrateData(anchorset = brain.anchors, normalization.method = "SCT", verbose = FALSE   )
brain.integrated <- RunPCA(brain.integrated, npcs = 100)
brain.integrated <- RunUMAP(brain.integrated, dims = 1:100)
brain.integrated <- RunTSNE(brain.integrated, dims = 1:100)
brain.integrated <- FindNeighbors(brain.integrated, dims = 1:100)
brain.integrated <- FindClusters(brain.integrated, resolution = 2)

saveRDS(brain.integrated, file = "../integrationseuratv3.rds")

#Saved as RDS file and used as input for the scclusteval package 

###### Cluster-stability analysis ###### 

# Install packages
library(scclusteval)
library(tidyverse)
library(patchwork)
library(Seurat)
library(scclusteval)
library(ComplexHeatmap)
options(max.print = 30)
devtools::install_github("crazyhottommy/scclusteval")


##Cells: Analysis 
setwd("~/first draft/scclustaleval/cellsseuratv3")
subsample_idents<- readRDS("gather_subsample.rds")
fullsample_idents<- readRDS("gather_full_sample.rds")

message(sprintf('PCs: %s',paste0(unique(subsample_idents$pc),collapse = ', ')))
message(sprintf('K: %s',paste0(unique(subsample_idents$k_param),collapse = ', ')))
message(sprintf('Resolution: %s',paste0(unique(subsample_idents$resolution),collapse = ', ')))
subsample_idents_list<- subsample_idents %>% 
  group_by(pc, resolution, k_param) %>% 
  nest()
message('Computing cluster stability ...')
stable_clusters<- subsample_idents_list %>%
  mutate(stable_cluster = map(data, ~ AssignStableCluster(.x$original_ident,
                                                          .x$recluster_ident,
                                                          jaccard_cutoff = 0.8,
                                                          method = "jaccard_percent", 
                                                          percent_cutoff = 0.8)))
message("done.")
ParameterSetScatterPlot(stable_clusters = stable_clusters,
                        fullsample_idents = fullsample_idents,
                        x_var = "k_param",
                        y_var = "percentage",
                        facet_rows = "resolution",
                        facet_cols = "pc") +
  ggtitle("Percentage of cells in stable clusters")+ylab('Percentage of cells in stable clusters')

# Now, get the total number of clusters vs percentage of "stable" cells --> select the parameters that maximize the number of cell populations
df <- dplyr::left_join(stable_clusters, fullsample_idents) %>% 
  dplyr::ungroup() %>% dplyr::mutate(total = map_dbl(stable_cluster, ~length(.x$stable_cluster))) %>% dplyr::mutate(stable = map_dbl(stable_cluster, ~.x$number_of_stable_cluster)) %>% dplyr::mutate(percentage = map2_dbl(original_ident_full, stable_cluster, function(x, y) CalculatePercentCellInStable(x,y$stable_cluster))) %>% dplyr::select(-data, -stable_cluster, -original_ident_full) %>% dplyr::mutate_if(is.character, function(x) as.factor(as.numeric(x))) %>% tidyr::gather(total:stable, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              key = "category", value = "number")#ggplot(dff,aes(x = percentage, y = number))+geom_point()+ylab('N stable clusters')+xlab('Percentage of cells assigned')+geom_smooth()+facet_grid(k_param~pc)
ggplot(df[df$category=='stable',],aes(x = percentage, y = number))+geom_point()+ylab('N stable clusters')+xlab('Percentage of cells assigned to stable clusters')+geom_smooth()

df%>%arrange(-percentage)
# solution maximizing the number of stable clusters:
m = df[df$category == 'stable',]
m$pc = as.integer(as.character(m$pc))
m$resolution = as.numeric(as.character(m$resolution))
m$percentage = as.numeric(as.character(m$percentage))
m$k_param = as.integer(as.character(m$k_param))
m$k_param
library(psych)
m = m[,-(ncol(m)-1)] 
colnames(m)[ncol(m)] = 'n_stable_cl'
pairs.panels(m, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = F,  # show density plots
             ellipses = F, # show correlation ellipses
)
most_stable = df[df$category == 'stable',]%>%arrange(-number)%>%filter(row_number()==1)
most_stable
most_perc = df%>%filter(category == 'stable')%>%arrange(-percentage)%>%filter(row_number() == 1)

most_stable
most_perc

#save the stable cluster identities
subsample_idents_list %>% ungroup() %>% mutate(id = row_number()) %>%
  filter(pc == 150, resolution == 1.4, k_param == 10)
subsample_idents_list$data[[52]]

AssignHighestJaccard(subsample_idents_list$data[[52]]$original_ident, 
                     subsample_idents_list$data[[52]]$recluster_ident)

JaccardRainCloudPlot(subsample_idents_list$data[[52]]$original_ident,
                     subsample_idents_list$data[[52]]$recluster_ident) + 
  geom_hline(yintercept = c(0.6, 0.75), linetype = 2) +
  xlab("cluster id w/ k=10 res=1.4 pc=150") 

#if you want to save the jaccard indices
jaccard <- AssignHighestJaccard(subsample_idents_list$data[[52]]$original_ident, 
                                subsample_idents_list$data[[52]]$recluster_ident)
df <- data.frame(jaccard)
write.csv(df,file = "jaccard_most_stable_clusters.csv")


#save the stable clusters with the most % of cells
subsample_idents_list %>% ungroup() %>% mutate(id = row_number()) %>%
  filter(pc == 150, resolution == 0.6, k_param == 16)
subsample_idents_list$data[[132]]
AssignHighestJaccard(subsample_idents_list$data[[132]]$original_ident, 
                     subsample_idents_list$data[[132]]$recluster_ident)

JaccardRainCloudPlot(subsample_idents_list$data[[132]]$original_ident,
                     subsample_idents_list$data[[132]]$recluster_ident) + 
  geom_hline(yintercept = c(0.6, 0.75), linetype = 2) +
  xlab("cluster id w/ k=16 res=0.6 pc=150") 

#if you want to save the jaccard indices
jaccard <- AssignHighestJaccard(subsample_idents_list$data[[132]]$original_ident, 
                                subsample_idents_list$data[[132]]$recluster_ident)
df <- data.frame(jaccard)
write.csv(df,file = "jaccard_most_%_cells.csv")

#Save the cluster identities
cells_clusterids_most_perc_cells <- fullsample_idents %>% mutate(id = row_number()) %>% filter(pc == 150, resolution == 0.6, k_param == 16)
cells_clusterids_most_stable_clusters <- fullsample_idents %>% mutate(id = row_number()) %>% filter(pc == 150, resolution == 1.4, k_param == 10)

cells_clusterids_most_perc_cells_list <- cells_clusterids_most_perc_cells$original_ident_full
cells_clusterids_most_stable_clusters_list <- cells_clusterids_most_stable_clusters$original_ident_full


##Nuclei: analysis 
setwd("~/first draft/scclustaleval/nucleisseuratv3")
subsample_idents<- readRDS("gather_subsample.rds")
fullsample_idents<- readRDS("gather_full_sample.rds")

message(sprintf('PCs: %s',paste0(unique(subsample_idents$pc),collapse = ', ')))
message(sprintf('K: %s',paste0(unique(subsample_idents$k_param),collapse = ', ')))
message(sprintf('Resolution: %s',paste0(unique(subsample_idents$resolution),collapse = ', ')))
subsample_idents_list<- subsample_idents %>% 
  group_by(pc, resolution, k_param) %>% 
  nest()
message('Computing cluster stability ...')
stable_clusters<- subsample_idents_list %>%
  mutate(stable_cluster = map(data, ~ AssignStableCluster(.x$original_ident,
                                                          .x$recluster_ident,
                                                          jaccard_cutoff = 0.8,
                                                          method = "jaccard_percent", 
                                                          percent_cutoff = 0.8)))
message("done.")
ParameterSetScatterPlot(stable_clusters = stable_clusters,
                        fullsample_idents = fullsample_idents,
                        x_var = "k_param",
                        y_var = "percentage",
                        facet_rows = "resolution",
                        facet_cols = "pc") +
  ggtitle("Percentage of cells in stable clusters")+ylab('Percentage of cells in stable clusters')

# e.g. pc=50, res=1.5, k=100, total = 22, stable = 14
df <- dplyr::left_join(stable_clusters, fullsample_idents) %>% 
  dplyr::ungroup() %>% dplyr::mutate(total = map_dbl(stable_cluster, ~length(.x$stable_cluster))) %>% dplyr::mutate(stable = map_dbl(stable_cluster, ~.x$number_of_stable_cluster)) %>% dplyr::mutate(percentage = map2_dbl(original_ident_full, stable_cluster, function(x, y) CalculatePercentCellInStable(x,y$stable_cluster))) %>% dplyr::select(-data, -stable_cluster, -original_ident_full) %>% dplyr::mutate_if(is.character, function(x) as.factor(as.numeric(x))) %>% tidyr::gather(total:stable, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              key = "category", value = "number")#ggplot(dff,aes(x = percentage, y = number))+geom_point()+ylab('N stable clusters')+xlab('Percentage of cells assigned')+geom_smooth()+facet_grid(k_param~pc)
ggplot(df[df$category=='stable',],aes(x = percentage, y = number))+geom_point()+ylab('N stable clusters')+xlab('Percentage of cells assigned to stable clusters')+geom_smooth()

df%>%arrange(-percentage)

most_stable = df[df$category == 'stable',]%>%arrange(-number)%>%filter(row_number()==1)
most_stable
most_perc = df%>%filter(category == 'stable')%>%arrange(-percentage)%>%filter(row_number() == 1)

most_stable
most_perc


#save the stable cluster identities
subsample_idents_list %>% ungroup() %>% mutate(id = row_number()) %>%
  filter(pc == 50, resolution == 3.0, k_param == 10)

#printing the entire tibble:
print(subsample_idents_list, n = nrow(subsample_idents_list))

# => 62 50 pca   3.0   res     10 k     <tibble [20 x 3]>


subsample_idents_list$data[[62]]
AssignHighestJaccard(subsample_idents_list$data[[62]]$original_ident, 
                     subsample_idents_list$data[[62]]$recluster_ident)

JaccardRainCloudPlot(subsample_idents_list$data[[62]]$original_ident,
                     subsample_idents_list$data[[62]]$recluster_ident) + 
  geom_hline(yintercept = c(0.6, 0.75), linetype = 2) +
  xlab("cluster id w/ k=10 res=3 pc=50") 


#if you want to save the jaccard indices
jaccard <- AssignHighestJaccard(subsample_idents_list$data[[62]]$original_ident, 
                                subsample_idents_list$data[[62]]$recluster_ident)
df <- data.frame(jaccard)
write.csv(df,file = "jaccard_most_stable_clusters.csv")

#save the stable clusters with the most % of cells
subsample_idents_list %>% ungroup() %>% mutate(id = row_number()) %>%
  filter(pc == 150, resolution == 0.6, k_param == 30)

subsample_idents_list$data[[196]]
AssignHighestJaccard(subsample_idents_list$data[[196]]$original_ident, 
                     subsample_idents_list$data[[196]]$recluster_ident)

JaccardRainCloudPlot(subsample_idents_list$data[[196]]$original_ident,
                     subsample_idents_list$data[[196]]$recluster_ident) + 
  geom_hline(yintercept = c(0.6, 0.75), linetype = 2) +
  xlab("cluster id w/ k=30 res=0.6 pc=150") 

#if you want to save the jaccard indices
jaccard <- AssignHighestJaccard(subsample_idents_list$data[[196]]$original_ident, 
                                subsample_idents_list$data[[196]]$recluster_ident)
df <- data.frame(jaccard)
write.csv(df,file = "jaccard_most_%_nuclei.csv")

#Save the cluster identities

nuclei_clusterids_most_perc_nuclei <- fullsample_idents %>% mutate(id = row_number()) %>% filter(pc == 150, resolution == 0.6, k_param == 30)
nuclei_clusterids_most_perc_nuclei_list <- nuclei_clusterids_most_perc_nuclei$original_ident_full

#printing the entire tibble
print(fullsample_idents, n = nrow(fullsample_idents))
nuclei_clusterids_most_stable_clusters <- fullsample_idents$original_ident_full[[62]] 
nuclei_clusterids_most_stable_clusters <- data.frame(nuclei_clusterids_most_stable_clusters)
nuclei_clusterids_most_stable_clusters_list <- list(nuclei_clusterids_most_stable_clusters_list)



####Integrated: Analysis
setwd("~/first draft/scclustaleval/integrationseuratv3")
subsample_idents<- readRDS("gather_subsample.rds")
fullsample_idents<- readRDS("gather_full_sample.rds")

message(sprintf('PCs: %s',paste0(unique(subsample_idents$pc),collapse = ', ')))
message(sprintf('K: %s',paste0(unique(subsample_idents$k_param),collapse = ', ')))
message(sprintf('Resolution: %s',paste0(unique(subsample_idents$resolution),collapse = ', ')))
subsample_idents_list<- subsample_idents %>% 
  group_by(pc, resolution, k_param) %>% 
  nest()
message('Computing cluster stability ...')
stable_clusters<- subsample_idents_list %>%
  mutate(stable_cluster = map(data, ~ AssignStableCluster(.x$original_ident,
                                                          .x$recluster_ident,
                                                          jaccard_cutoff = 0.8,
                                                          method = "jaccard_percent", 
                                                          percent_cutoff = 0.8)))
message("done.")
ParameterSetScatterPlot(stable_clusters = stable_clusters,
                        fullsample_idents = fullsample_idents,
                        x_var = "k_param",
                        y_var = "percentage",
                        facet_rows = "resolution",
                        facet_cols = "pc") +
  ggtitle("Percentage of cells in stable clusters")+ylab('Percentage of cells in stable clusters')

# Now, get the total number of clusters vs percentage of "stable" cells --> select the parameters that maximize the number of cell populations
df <- dplyr::left_join(stable_clusters, fullsample_idents) %>% 
  dplyr::ungroup() %>% dplyr::mutate(total = map_dbl(stable_cluster, ~length(.x$stable_cluster))) %>% dplyr::mutate(stable = map_dbl(stable_cluster, ~.x$number_of_stable_cluster)) %>% dplyr::mutate(percentage = map2_dbl(original_ident_full, stable_cluster, function(x, y) CalculatePercentCellInStable(x,y$stable_cluster))) %>% dplyr::select(-data, -stable_cluster, -original_ident_full) %>% dplyr::mutate_if(is.character, function(x) as.factor(as.numeric(x))) %>% tidyr::gather(total:stable, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              key = "category", value = "number")#ggplot(dff,aes(x = percentage, y = number))+geom_point()+ylab('N stable clusters')+xlab('Percentage of cells assigned')+geom_smooth()+facet_grid(k_param~pc)
ggplot(df[df$category=='stable',],aes(x = percentage, y = number))+geom_point()+ylab('N stable clusters')+xlab('Percentage of cells assigned to stable clusters')+geom_smooth()

df%>%arrange(-percentage)

# solution maximizing the number of stable clusters:
m = df[df$category == 'stable',]
m$pc = as.integer(as.character(m$pc))
m$resolution = as.numeric(as.character(m$resolution))
m$percentage = as.numeric(as.character(m$percentage))
m$k_param = as.integer(as.character(m$k_param))
library(psych)
m = m[,-(ncol(m)-1)] 
colnames(m)[ncol(m)] = 'n_stable_cl'
pairs.panels(m, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = F,  # show density plots
             ellipses = F, # show correlation ellipses
)
most_stable = df[df$category == 'stable',]%>%arrange(-number)%>%filter(row_number()==1)
most_perc = df%>%filter(category == 'stable')%>%arrange(-percentage)%>%filter(row_number() == 1)

most_stable
most_perc


#save the stable cluster identities
subsample_idents_list %>% ungroup() %>% mutate(id = row_number()) %>%
  filter(pc == 150, resolution == 2, k_param == 10)

#printing the entire tibble:
print(subsample_idents_list, n = nrow(subsample_idents_list))

# => 60 150 pca   2  res     10 k     <tibble [20 x 3]>


subsample_idents_list$data[[60]]

AssignHighestJaccard(subsample_idents_list$data[[60]]$original_ident, 
                     subsample_idents_list$data[[60]]$recluster_ident)

JaccardRainCloudPlot(subsample_idents_list$data[[60]]$original_ident,
                     subsample_idents_list$data[[60]]$recluster_ident) + 
  geom_hline(yintercept = c(0.6, 0.75), linetype = 2) +
  xlab("cluster id w/ k=10 res=2 pc=150") 

#if you want to save the jaccard indices
jaccard <- AssignHighestJaccard(subsample_idents_list$data[[60]]$original_ident, 
                                subsample_idents_list$data[[60]]$recluster_ident)
df <- data.frame(jaccard)
write.csv(df,file = "jaccard_most_stable_clusters.csv")

#save the cluster identities
print(fullsample_idents, n = nrow(fullsample_idents))
brain_clusterids_most_stable_clusters <- fullsample_idents$original_ident_full[[60]] 
brain_clusterids_most_stable_clusters <- data.frame(brain_clusterids_most_stable_clusters)
brain_clusterids_most_stable_clusters_list <- list(brain_clusterids_most_stable_clusters)

#save the stable clusters with the most % of cells
subsample_idents_list %>% ungroup() %>% mutate(id = row_number()) %>%
  filter(pc == 30, resolution == 0.6, k_param == 30)

subsample_idents_list$data[[193]]
AssignHighestJaccard(subsample_idents_list$data[[193]]$original_ident, 
                     subsample_idents_list$data[[193]]$recluster_ident)

JaccardRainCloudPlot(subsample_idents_list$data[[193]]$original_ident,
                     subsample_idents_list$data[[193]]$recluster_ident) + 
  geom_hline(yintercept = c(0.6, 0.75), linetype = 2) +
  xlab("cluster id w/ k=30 res=0.6 pc=30") 

#if you want to save the jaccard indices
jaccard <- AssignHighestJaccard(subsample_idents_list$data[[193]]$original_ident, 
                                subsample_idents_list$data[[193]]$recluster_ident)
df <- data.frame(jaccard)
write.csv(df,file = "jaccard_most_%_brain.csv")

#Save the cluster identities
brain_clusterids_most_perc_cells <- fullsample_idents %>% mutate(id = row_number()) %>% filter(pc == 30, resolution == 0.6, k_param == 30)
brain_clusterids_most_perc_cells_list <- brain_clusterids_most_perc_cells$original_ident_full

########################################################################
#This is the end of the workflow from the scclusteval package but you can now use the cluster identities in seurat
########## cells ########## 
cells <- readRDS("cellsseuratv3.rds")
########## Cluster with optimal parameters for the most stable clusters
DefaultAssay(cells) <- "SCT"
cells <- RunPCA(cells, verbose = FALSE, npcs = 150)
cells <- RunUMAP(cells, dims = 1:150, verbose = FALSE)
cells <- RunTSNE(cells, dims = 1:150, verbose = FALSE)
cells <- FindNeighbors(cells, k.param = 10, dims = 1:150, verbose = FALSE)
cells <- FindClusters(cells, resolution = 1.4, verbose = FALSE)
##Retrieve the cell identities from the scclusteval package
Idents(cells) = cells_clusterids_most_stable_clusters_list

##Check if all cells are assigned to clusters
table(is.na(Idents(cells)))

##Check your total cluster number and save the umaps/tsne plots as pdf (A4 landscape - with and without labels)
DimPlot(cells, reduction = "umap", label=TRUE)
DimPlot(cells, reduction = "tsne", label=TRUE)

DimPlot(cells, reduction = "umap", label=FALSE)
DimPlot(cells, reduction = "tsne", label=FALSE)

##Save some parameters
FeaturePlot(cells, features="percent.mt", reduction ="tsne", label=FALSE)
FeaturePlot(cells, features="percent.mt", reduction ="umap", label=FALSE)

##Export these plots as pdf A4 -landscape + adjust image height to 4
FeaturePlot(cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),reduction ="tsne", ncol = 3)
FeaturePlot(cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),reduction ="umap", ncol = 3)

###For DE analysis it is best to use the RNA assay
DefaultAssay(cells) <- "RNA"

###This assay should also be normalized and scaled
cells <- NormalizeData(cells, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(cells)
cells <- ScaleData(cells, features = all.genes)

cells.markers <- FindAllMarkers(cells, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

cells.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10

DoHeatmap(cells, features = top10$gene, label = FALSE) + NoLegend()

saveRDS(cells, file = "cells_most_stable_clusters.rds")

###Save marker genes as csv

df <- data.frame(cells.markers)
write.csv(df,file = "cells_most_stable_clusters_markers.csv")

########## Cluster with optimal parameters for the highest number of stable clusters
DefaultAssay(cells) <- "SCT"
cells <- RunPCA(cells, verbose = FALSE, npcs = 150)
cells <- RunUMAP(cells, dims = 1:150, verbose = FALSE)
cells <- RunTSNE(cells, dims = 1:150, verbose = FALSE)
cells <- FindNeighbors(cells, k.param = 16, dims = 1:150, verbose = FALSE)
cells <- FindClusters(cells, resolution = 0.6, verbose = FALSE)
##Retrieve the cell identities from the scclusteval package
Idents(cells) = cells_clusterids_most_perc_cells_list

##Check if all cells are assigned to clusters
table(is.na(Idents(cells)))

setwd("~/first draft/seurat/cells/most_perc_cells")
##Check your total cluster number and save the umaps/tsne plots as pdf (A4 landscape - with and without labels)
DimPlot(cells, reduction = "umap", label=TRUE)
DimPlot(cells, reduction = "tsne", label=TRUE)

DimPlot(cells, reduction = "umap", label=FALSE)
DimPlot(cells, reduction = "tsne", label=FALSE)

##Save some parameters
FeaturePlot(cells, features="percent.mt", reduction ="tsne", label=FALSE)
FeaturePlot(cells, features="percent.mt", reduction ="umap", label=FALSE)

##Export these plots as pdf A4 -landscape + adjust image height to 4
FeaturePlot(cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),reduction ="tsne", ncol = 3)
FeaturePlot(cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),reduction ="umap", ncol = 3)

###For DE analysis it is best to use the RNA assay
DefaultAssay(cells) <- "RNA"

###This assay should also be normalized and scaled
cells <- NormalizeData(cells, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(cells)
cells <- ScaleData(cells, features = all.genes)

cells.markers <- FindAllMarkers(cells, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

cells.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10

DoHeatmap(cells, features = top10$gene, label = FALSE) + NoLegend()

saveRDS(cells, file = "cells_most_perc_cells.rds")

###Save marker genes as csv

df <- data.frame(cells.markers)
write.csv(df,file = "cells_most_perc_cells_markers.csv")

########## We repeated the same for the nuclei.

########## Preprocessing of integration object ########## 

########## Cluster with optimal parameters for the highest number of stable clusters
brain <- readRDS("integrationseuratv3.rds")
DefaultAssay(brain) <- "integrated"
brain <- RunPCA(brain, npcs = 150)
brain <- RunUMAP(brain, dims = 1:150)
brain <- RunTSNE(brain, dims = 1:150)
brain <- FindNeighbors(brain, k.param = 10, dims = 1:150)
brain <- FindClusters(brain, resolution = 2)

Idents(brain) = brain_clusterids_most_stable_clusters_list

##Check if all cells are assigned to clusters
table(is.na(Idents(brain)))

##Check your total cluster number and save the umaps/tsne plots as pdf (A4 landscape - with and without labels)
DimPlot(brain, reduction = "umap", label=TRUE)
DimPlot(brain, reduction = "tsne", label=TRUE)

DimPlot(brain, reduction = "umap", label=FALSE)
DimPlot(brain, reduction = "tsne", label=FALSE)

##Save some parameters
FeaturePlot(brain, features="percent.mt", reduction ="tsne", label=FALSE)
FeaturePlot(brain, features="percent.mt", reduction ="umap", label=FALSE)

##Export these plots as pdf A4 -landscape + adjust image height to 4
FeaturePlot(brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),reduction ="tsne", ncol = 3)
FeaturePlot(brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),reduction ="umap", ncol = 3)

###For DE analysis it is best to use the RNA assay
DefaultAssay(brain) <- "RNA"

###This assay should also be normalized and scaled
brain <- NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(brain)
memory.limit()
brain <- ScaleData(brain, features = all.genes)

brain.markers <- FindAllMarkers(brain, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

brain.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- brain.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10

DoHeatmap(brain, features = top10$gene, label = FALSE) + NoLegend()

saveRDS(brain, file = "brain_most_stable_clusters.rds")

###Save marker genes as csv

df <- data.frame(brain.markers)
write.csv(df,file = "brain_most_stable_clusters_markers.csv")

###Prepare plots that show the data integration
setwd("~/first draft/seurat/integration/most_stable")
brain <- readRDS("brain_most_stable_clusters.rds")

DimPlot(brain, reduction = "umap", label=FALSE, split.by = "orig.ident")
DimPlot(brain, reduction = "tsne", label=FALSE, split.by = "orig.ident")

DimPlot(brain, reduction = "umap", label=FALSE, group.by = "orig.ident")
DimPlot(brain, reduction = "tsne", label=FALSE, group.by = "orig.ident")

########## Cluster with optimal parameters for the highest percentage of cells assigned to stable clusters

DefaultAssay(brain) <- "integrated"
brain <- RunPCA(brain, verbose = FALSE, npcs = 30)
brain <- RunUMAP(brain, dims = 1:30, verbose = FALSE)
brain <- RunTSNE(brain, dims = 1:30, verbose = FALSE)
brain <- FindNeighbors(brain, k.param = 30, dims = 1:30, verbose = FALSE)
brain <- FindClusters(brain, resolution = 0.6, verbose = FALSE)
Idents(brain) = brain_clusterids_most_perc_cells_list

##Check if all cells are assigned to clusters
table(is.na(Idents(brain)))

setwd("~/first draft/seurat/integration/most_perc_cells/clusters")
##Check your total cluster number and save the umaps/tsne plots as pdf (A4 landscape - with and without labels)
DimPlot(brain, reduction = "umap", label=TRUE)
DimPlot(brain, reduction = "tsne", label=TRUE)

DimPlot(brain, reduction = "umap", label=FALSE)
DimPlot(brain, reduction = "tsne", label=FALSE)

##Save some parameters
FeaturePlot(brain, features="percent.mt", reduction ="tsne", label=FALSE)
FeaturePlot(brain, features="percent.mt", reduction ="umap", label=FALSE)

##Export these plots as pdf A4 -landscape + adjust image height to 4
FeaturePlot(brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),reduction ="tsne", ncol = 3)
FeaturePlot(brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),reduction ="umap", ncol = 3)

###For DE analysis it is best to use the RNA assay
DefaultAssay(brain) <- "RNA"

###This assay should also be normalized and scaled
brain <- NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(brain)
brain <- ScaleData(brain, features = all.genes)
brain.markers <- FindAllMarkers(brain, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

brain.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- brain.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10

DoHeatmap(brain, features = top10$gene, label = FALSE) + NoLegend()

saveRDS(brain, file = "brain_most_perc_cells.rds")

###Save marker genes as csv

df <- data.frame(brain.markers)
write.csv(df,file = "brain_most_perc_nuclei_markers.csv")


###Annotate stable clusters 
brain_most_stable_clusters <- readRDS("brain_most_stable_clusters.rds")
annotation <- c("not separated","unstable","unstable","IGL2-GLUT/DOP","GLIA1","OGL2-DOP","unstable","unstable","not separated","DOP2","FBL","SUB","not separated","unstable","unstable","IGL2-GLUT/DOP","unstable","not separated","unstable","OGL1","unstable","GLUT1","unstable","unstable","unstable","DOP3","SERT","GLUT3","DOP1","IGL3","ACH2","unstable","PREC","PEP-Fmrfa3","unstable","PEP-Fmrfa1","unstable","unstable","VL","unstable","unstable","unstable","unstable","GABA","unstable","unstable","IGL1-OA","unstable","not separated","unstable","GLUT2","OA","ACH3","ACH1","EC","unstable","IGL4-L11","GLIA2","not separated","GLUT4","TBA1","OGL3-OA","CCAP","unstable","PEP-APWG","unstable","TBA2","TBA8","unstable","unstable","HC","unstable","TBA3","unstable","TBA4","TBA5","PEP-Burs","TBA6","unstable","unstable","unstable","unstable","GLIA3","unstable","unstable","TBA7","unstable")
names(annotation) <- levels(brain_most_stable_clusters)
brain_most_stable_clusters_annotated <- RenameIdents(brain_most_stable_clusters, annotation)
saveRDS(brain_most_stable_clusters_annotated, file = "brain_most_stable_clusters_w_annotation.rds")
DimPlot(brain_most_stable_clusters_annotated, reduction = "tsne", label = TRUE, pt.size = 0.5)

###Find all markers for seurat clusters
markers <- FindAllMarkers(brain_most_stable_clusters,only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
###Find all markers for stable clusters
stable.markers <- FindAllMarkers(stable_clusters, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)


###Prepare data for GEO submission
#CELL LEVEL METADATA
#Extract cell barcodes per cell type from seurat object
#Read in annotated seurat object
brain <- readRDS("brain_most_stable_clusters_w_annotation.rds")
brain_annotated_cells <- brain@active.ident
brain_annotated_cells
write.csv(brain_annotated_cells,'cell_level_metadata.csv')

#COUNT MATRIX
#Extract raw count matrix
countmatrix <- brain@assays$RNA@counts
write.csv(countmatrix,'countmatrix.csv')

#EXTRACT EMBEDDING
so <- readRDS("seurat_object_w_annotation.rds")
tsneembedding <- so@reductions$tsne@cell.embeddings
write.csv(tsneembedding,'tsneembedding.csv')


#Proportion of cells in each cluster
table(Idents(brain))
prop.table(table(Idents(brain)))
library(Seurat)
#Proportion per cluster from each replicate
table(Idents(brain), brain$orig.ident)
prop <- prop.table(table(Idents(brain), brain$orig.ident))
prop <- as.data.frame(prop)
colnames(prop) <- c("Cluster", "Batch", "Frequency")
ggplot(prop, aes(x=Cluster, y=Frequency, fill=Batch)) + geom_bar(stat="identity") + coord_flip()
#Saved as portrait A4

#do this also for all clusters and add this to the other plot
brain <- readRDS("brain_most_stable_clusters_w_annotation.rds")
brain@active.ident <- brain@meta.data$seurat_clusters
prop <- table(Idents(brain), brain$orig.ident)
prop <- prop.table(table(Idents(brain), brain$orig.ident))
prop <- as.data.frame(prop)
colnames(prop) <- c("Cluster", "Batch", "Cells")
g2 <- ggplot(prop, aes(x=Cluster, y=Cells, fill=Batch)) + geom_bar(stat="identity") + coord_flip()
g2 <- g2 + theme(legend.position = "none")
g2

#Exact numbers per cluster from each replicate
exact <- table(Idents(brain), brain$orig.ident)
ex <- as.data.frame(exact)
colnames(ex) <- c("Cluster", "Batch", "Cells")
g4 <- ggplot(ex, aes(x=Cluster, y=Cells, fill=Batch)) + geom_bar(stat="identity") + coord_flip()
g4 <- g4 + theme(legend.position="top",legend.justification="right")      
g4

#Make boxplots for jaccard indices 
jaccardindices <- as.data.frame(jaccard_most_stable_clusters)
colnames(jaccard_most_stable_clusters) <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86")
df <- jaccardindices %>% select("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86") %>% pivot_longer(., cols = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86"), names_to = "Cluster", values_to = "Jaccard")
df$Cluster <- as.numeric(df$Cluster)
df <- df[order(df$Cluster), ]
df <- as.data.frame(df)
g2 <- ggplot(prop, aes(x=Cluster, y=Frequency, fill=Batch)) + geom_bar(stat="identity") + coord_flip()


g1 <- ggplot(df, aes(x=factor(Cluster),y=Jaccard)) + 
  geom_violin(fill = "grey90") + 
  geom_boxplot(width = .2, outlier.shape = NA, coef = 0)
g1 <- g1 + geom_hline(yintercept = c(0.6)) + coord_flip() + xlab("Cluster")
g3 <- g2 + g1 + g4


