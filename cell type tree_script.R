#more information on gitlab
https://git.embl.de/musser/profiling-cellular-diversity-in-sponges-informs-animal-cell-type-and-nervous-system-evolution/-/blob/master/single_cell_analysis/cell_type_clustering.R
install.packages('Seurat')
library(Seurat)
install.packages('ape')
library(ape)
install.packages('phytools')
library(phytools)
install.packages('ggtree')
library(ggtree)  
install.packages('treeio')
library(treeio)
install.packages('tidytree')
library(tidytree)  
install.packages('Rtools')
install.packages('phytools')
install.packages('tmvnsim')
library(tmvnsim)
library(phytools)
#Read in data
brain <- readRDS("seurat_object_with_annotation.rds")
stableclusters <- readRDS("stable_clusters.rds")

#average gene expression per cell type
DefaultAssay(stableclusters) <- "SCT"
cluster.averages <- AverageExpression(stableclusters)
types_average_norm <- as.matrix(cluster.averages$SCT)
allgenes <- rownames(stableclusters)
allgenes_matrix <- t(types_average_norm[allgenes,])
#calculate the euclidean distance matrix
allgenes_matrix_dist <- dist(allgenes_matrix, method = "euclidean")
#neighbor-joining tree estimation with the ape package
allgenes_matrix_tree <- nj(allgenes_matrix_dist)
#plot phylogenetic tree; unrooted vs midpoint root
plot.phylo(allgenes_matrix_tree, type = "unrooted")
plot.phylo(midpoint.root(allgenes_matrix_tree))
#midpoint root tree
allgenes_matrix_tree <- midpoint.root(allgenes_matrix_tree)
#perform bootstrap analysis with phytools
allgenes_matrix_tree_boot <- boot.phylo(allgenes_matrix_tree, allgenes_matrix, FUN = function(xx) nj(dist(xx)), B = 10000)
?boot.phylo()
#visualize percentages on the nodes
allgenes_matrix_tree$node.label <- allgenes_matrix_tree_boot/100
allgenes_matrix_tree$node.label <- round(allgenes_matrix_tree$node.label)
#plot tree with bootstrap values
allgenes_matrix_tree
plot.phylo(allgenes_matrix_tree)
#save as text file and open in figtree
write.tree(allgenes_matrix_tree)

#reconstruct tree based only on HVG
DefaultAssay(stableclusters) <- "RNA"
variable_genes <- FindVariableFeatures(stableclusters, selection.method = "vst", nfeatures = 2000)
plot1 <- VariableFeaturePlot(variable_genes)
#average gene expression per cell type
#subset matrix for hvg
hvg_id <- variable_genes$RNA@var.features
hvg <- subset(variable_genes, features = hvg_id)
hvg
#average gene expression per cell type
DefaultAssay(hvg) <- "SCT"
hvgid <- rownames(hvg$SCT)
cluster.averages <- AverageExpression(hvg)
types_average_norm <- as.matrix(cluster.averages$SCT)
hvg_matrix <- t(types_average_norm[hvgid,])
#calculate the euclidean distance matrix
hvggenes_matrix_dist <- dist(hvg_matrix, method = "euclidean")
#neighbor-joining tree estimation with the ape package
hvggenes_matrix_tree <- nj(hvggenes_matrix_dist)
#plot phylogenetic tree; unrooted vs midpoint root
plot.phylo(hvggenes_matrix_tree, type = "unrooted")
plot.phylo(midpoint.root(hvggenes_matrix_tree))
#midpoint root tree
hvggenes_matrix_tree <- midpoint.root(hvggenes_matrix_tree)
#perform bootstrap analysis with phytools
hvggenes_matrix_tree_boot <- boot.phylo(hvggenes_matrix_tree, hvg_matrix, FUN = function(xx) nj(dist(xx)), B = 10000)
?boot.phylo()
#visualize percentages on the nodes
hvggenes_matrix_tree$node.label <- hvggenes_matrix_tree_boot/100
hvggenes_matrix_tree$node.label <- round(hvggenes_matrix_tree$node.label)
#plot tree with bootstrap values
hvggenes_matrix_tree
plot.phylo(hvggenes_matrix_tree)
#save as text file and open in figtree
write.tree(hvggenes_matrix_tree)

####create cell type tree for all TF
tf <- read_excel("~//phylogenetic trees/tf.xlsx", +col_names = FALSE)
TF <- tf$...1
alltf <- subset(stableclusters, features = TF)

#average gene expression per cell type
DefaultAssay(alltf) <- "SCT"
cluster.averages <- AverageExpression(alltf)
types_average_norm <- as.matrix(cluster.averages$SCT)
alltfgeneids <- rownames(alltf)
alltf_matrix <- t(types_average_norm[alltfgeneids,])
#calculate the euclidean distance matrix
alltf_matrix_dist <- dist(alltf_matrix, method = "euclidean")
#neighbor-joining tree estimation with the ape package
alltf_matrix_tree <- nj(alltf_matrix_dist)
#plot phylogenetic tree; unrooted vs midpoint root
plot.phylo(alltf_matrix_tree, type = "unrooted")
plot.phylo(midpoint.root(alltf_matrix_tree))
#midpoint root tree
alltf_matrix_tree <- midpoint.root(alltf_matrix_tree)
#perform bootstrap analysis with phytools
alltf_matrix_tree_boot <- boot.phylo(alltf_matrix_tree, alltf_matrix, FUN = function(xx) nj(dist(xx)), B = 10000)
?boot.phylo()
#visualize percentages on the nodes
alltf_matrix_tree$node.label <- alltf_matrix_tree_boot/100
alltf_matrix_tree$node.label <- round(alltf_matrix_tree$node.label)
#plot tree with bootstrap values
alltf_matrix_tree
plot.phylo(alltf_matrix_tree)
#save as text file and open in figtree
write.tree(alltf_matrix_tree)

#reconstruct tree based only on HVG TF
DefaultAssay(alltf) <- "RNA"
alltf <- FindVariableFeatures(alltf, selection.method = "vst", nfeatures = 200)

#subset matrix for hvg
hvtf_id <- alltf$RNA@var.features
hvtf <- subset(alltf, features = hvtf_id)
hvtf

#average gene expression per cell type
DefaultAssay(hvtf) <- "SCT"
hvtfid <- rownames(hvtf$SCT)
cluster.averages <- AverageExpression(hvtf)
types_average_norm <- as.matrix(cluster.averages$SCT)
hvtf_matrix <- t(types_average_norm[hvtfid,])

#calculate the euclidean distance matrix
hvtf_matrix_dist <- dist(hvtf_matrix, method = "euclidean")
#neighbor-joining tree estimation with the ape package
hvtf_matrix_tree <- nj(hvtf_matrix_dist)
#plot phylogenetic tree; unrooted vs midpoint root
plot.phylo(hvtf_matrix_tree, type = "unrooted")
plot.phylo(midpoint.root(hvtf_matrix_tree))
#midpoint root tree
hvtf_matrix_tree <- midpoint.root(hvtf_matrix_tree)
#perform bootstrap analysis with phytools
hvtf_matrix_tree_boot <- boot.phylo(hvtf_matrix_tree, hvtf_matrix, FUN = function(xx) nj(dist(xx)), B = 10000)
?boot.phylo()
#visualize percentages on the nodes
hvtf_matrix_tree$node.label <- hvtf_matrix_tree_boot/100
hvtf_matrix_tree$node.label <- round(hvtf_matrix_tree$node.label)
#plot tree with bootstrap values
hvtf_matrix_tree
plot.phylo(hvtf_matrix_tree)
#save as text file and open in figtree
write.tree(hvtf_matrix_tree)

DoHeatmap(alltf, features = hvtf_id) + NoLegend()



#create tree based on 390 tf markers
DefaultAssay(alltf) <- "RNA"
TFs.markers <- FindAllMarkers(alltf, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
tfmarkers <- TFs.markers$gene
write.csv(tfmarkers, file = "tfmarkers.csv")

#subset matrix for tfmarkers
tf_markers <- subset(alltf, features = tfmarkers)
?subset()

#average gene expression per cell type
DefaultAssay(tf_markers) <- "SCT"
tf_geneids <- rownames(tf_markers)
cluster.averages <- AverageExpression(tf_markers)
types_average_norm <- as.matrix(cluster.averages$SCT)
tfmarkers_matrix <- t(types_average_norm[tf_geneids,])

#calculate the euclidean distance matrix
tfmarkers_matrix_dist <- dist(tfmarkers_matrix, method = "euclidean")
#neighbor-joining tree estimation with the ape package
tfmarkers_matrix_tree <- nj(tfmarkers_matrix_dist)
#plot phylogenetic tree; unrooted vs midpoint root
plot.phylo(tfmarkers_matrix_tree, type = "unrooted")
plot.phylo(midpoint.root(tfmarkers_matrix_tree))
#midpoint root tree
tfmarkers_matrix_tree <- midpoint.root(tfmarkers_matrix_tree)
#perform bootstrap analysis with phytools
tfmarkers_matrix_tree_boot <- boot.phylo(tfmarkers_matrix_tree, tfmarkers_matrix, FUN = function(xx) nj(dist(xx)), B = 10000)
?boot.phylo()
#visualize percentages on the nodes
tfmarkers_matrix_tree$node.label <- tfmarkers_matrix_tree_boot/100
tfmarkers_matrix_tree$node.label <- round(tfmarkers_matrix_tree$node.label)
#plot tree with bootstrap values
tfmarkers_matrix_tree
plot.phylo(tfmarkers_matrix_tree)
#save as text file and open in figtree
write.tree(tfmarkers_matrix_tree)

#looking for branch tf markers
a.markers <- FindMarkers(alltf, ident.1 = c("GLUT3","GLUT2","TBA4","OGL3-OA"), min.pct = 0.50, only.pos = TRUE,logfc.threshold = 0.25)
head(a.markers, n = 5)
b.markers <- FindMarkers(alltf, ident.1 = c("IGL2-GLUT/DOP","PEP-Burs","OGL1"), min.pct = 0.50, only.pos = TRUE, logfc.threshold = 0.25)
head(b.markers, n = 5)
c.markers <- FindMarkers(alltf, ident.1 = c("IGL3","TBA2","TBA7"), min.pct = 0.50, only.pos = TRUE, logfc.threshold = 0.25)
head(c.markers, n = 5)
d.markers <- FindMarkers(alltf, ident.1 = c("DOP1","OA","OGL2-DOP","DOP2","DOP3"), min.pct = 0.50, only.pos = TRUE, logfc.threshold = 0.25)
head(d.markers, n = 5)
e.markers <- FindMarkers(alltf, ident.1 = c("IGL1-OA","ACH1","ACH3"), min.pct = 0.50, only.pos = TRUE, logfc.threshold = 0.25)
head(e.markers, n = 5)
f.markers <- FindMarkers(alltf, ident.1 = c("PEP-Fmrfa1","ACH2"), min.pct = 0.50, only.pos = TRUE, logfc.threshold = 0.25)
head(f.markers, n = 5)
g.markers <- FindMarkers(alltf, ident.1 = c("GABA","GLUT1","SUB","SERT","FBL"), min.pct = 0.50, only.pos = TRUE, logfc.threshold = 0.25)
head(g.markers, n = 5)
h.markers <- FindMarkers(alltf, ident.1 = c("GLUT4","VL","TBA5","CCAP","TBA3","TBA1"), min.pct = 0.50, only.pos = TRUE, logfc.threshold = 0.25)
head(h.markers, n = 5)
i.markers <- FindMarkers(alltf, ident.1 = c("TBA8","PREC"), min.pct = 0.50, only.pos = TRUE, logfc.threshold = 0.25)
head(i.markers, n = 5)
j.markers <- FindMarkers(alltf, ident.1 = c("GLIA1","GLIA2","GLIA3"), min.pct = 0.50, only.pos = TRUE, logfc.threshold = 0.25)
head(j.markers, n = 5)
k.markers <- FindMarkers(alltf, ident.1 = c("HC","EC"), min.pct = 0.50, only.pos = TRUE, logfc.threshold = 0.25)
head(k.markers, n = 5)
jk.markers <- FindMarkers(alltf, ident.1 = c("GLIA1","GLIA2","GLIA3","HC","EC"), min.pct = 0.50, only.pos = TRUE, logfc.threshold = 0.25)
head(jk.markers, n = 5)
ijk.markers <- FindMarkers(alltf, ident.1 = c("GLIA1","GLIA2","GLIA3","HC","EC","TBA8","PREC"), min.pct = 0.50, only.pos = TRUE, logfc.threshold = 0.25)
head(ijk.markers, n = 5)

ab.markers <- FindMarkers(alltf, ident.1 = c("GLUT3","GLUT2","TBA4","OGL3-OA","IGL2-GLUT/DOP","PEP-Burs","OGL1"), min.pct = 0.50, only.pos = TRUE,logfc.threshold = 0.25)
head(ab.markers, n = 5)
abc.markers <- FindMarkers(alltf, ident.1 = c("GLUT3","GLUT2","TBA4","OGL3-OA","IGL2-GLUT/DOP","PEP-Burs","OGL1","IGL3","TBA2","TBA7"), min.pct = 0.50, only.pos = TRUE,logfc.threshold = 0.25)
head(abc.markers, n = 5)
de.markers <- FindMarkers(alltf, ident.1 = c("DOP1","OA","OGL2-DOP","DOP2","DOP3","IGL1-OA","ACH1","ACH3"), min.pct = 0.50, only.pos = TRUE, logfc.threshold = 0.25)
head(de.markers, n = 5)

show_col(hue_pal()(12))

#reorder clusters
levels(alltf)
my_levels <- c("PEP-APWG","IGL4-L11","TBA6","HC","EC", "GLIA3", "GLIA2", "GLIA1", "PREC", "TBA8", "GLUT4","VL", "TBA5", "CCAP", "TBA3", "TBA1", "PEP-Fmrfa3", "GABA", "GLUT1", "SUB", "SERT", "FBL", "PEP-Fmrfa1", "ACH2", "IGL1-OA","ACH1","ACH3","DOP3","DOP2","OGL2-DOP","OA","DOP1","TBA7","TBA2","IGL3", "OGL1","PEP-Burs","IGL2-GLUT/DOP","OGL3-OA","TBA4","GLUT2","GLUT3")
factor(Idents(alltf), levels= my_levels)
Idents(alltf) <- factor(Idents(alltf), levels= my_levels)
levels(alltf)
#dotplot important TF
tfasmarkers <- c("LOC115209324","LOC115219287","LOC115229657","LOC115213193","LOC115217234","LOC115224734","LOC115219411","LOC115213110","LOC115211031","LOC115220806","LOC115221485","LOC115210669","LOC115214204","LOC115210326","LOC115217836","LOC115217255", "LOC115210677","LOC115211946","LOC115224244", "LOC115218995","LOC115218419","LOC115222994")
x <- DotPlot(alltf, features = tfasmarkers) + RotatedAxis()
x + coord_flip() + theme(legend.position="top")
x
