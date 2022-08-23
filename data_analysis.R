#Subset data for stable clusters only
brain <- readRDS("brain_most_stable_clusters_w_annotation.rds")
stable_clusters <- subset(brain, idents = c("IGL2-GLUT/DOP","GLIA1","OGL2-DOP","DOP2","FBL","SUB","IGL2-GLUT/DOP","OGL1","GLUT1","DOP3","SERT","GLUT3","DOP1","IGL3","ACH2","PREC","PEP-Fmrfa3","PEP-Fmrfa1","VL","GABA","IGL1-OA","GLUT2","OA","ACH3","ACH1","EC","PEP-L11","GLIA2","GLUT4","TBA1","OGL3-OA","CCAP","PEP-APWG","TBA2","PRECL","HC","TBA3","TBA4","TBA5","PEP-Burs","TBA6","GLIA3","TBA7"))
saveRDS(stable_clusters, file = "brain_only_stable_clusters.rds")
brain <- readRDS("brain_only_stable_clusters.rds")

### Cells versus nuclei ###
library(Seurat)
library(tidyverse)
library(ggplot2)
library(viridis)

#read in data
stableclusters <- readRDS("stable_clusters.rds")
stableclusters[["annotation"]] <- Idents(object = stableclusters)

#Set orig.ident as active ident
stableclusters <- SetIdent(stableclusters, value = stableclusters@meta.data$orig.ident)
levels(stableclusters)

#DE gene expression based on orig.ident
cellsvsnuclei.markers <- FindAllMarkers(stableclusters, min.pct = 0.25)
write.csv(cellsvsnuclei.markers, "cellsvsnuclei_markers.csv")
head(cellsvsnuclei.markers, n =10)

DE_list_cells_nuclei <- c("LOC115209225","LOC115224220","LOC115232632","LOC115212946","LOC115231102","LOC115216616","LOC115218233","LOC115214029","LOC115219147","LOC115213862","LOC115228388","LOC118762678","LOC118767422","LOC115216394","LOC115216632","LOC115215203","LOC115209421","LOC118762564","LOC118763604","LOC115214900")

#dotplot top 10 DE genes between cells and nuclei
x <- DotPlot(stableclusters, features = DE_list_cells_nuclei) + RotatedAxis()
x  + coord_flip() + theme(legend.position="top")
stableclusters <- SetIdent(stableclusters, value = stableclusters@meta.data$annotation)

#doplot IEG
ieg <- c("LOC115222465","LOC115214910","LOC115212294","LOC115228662")
x <- DotPlot(stableclusters, features = ieg, split.by = 'orig.ident') + RotatedAxis()
x  + coord_flip() + theme(legend.position="top")
stableclusters <- SetIdent(stableclusters, value = stableclusters@meta.data$annotation)
stableclusters <- SetIdent(stableclusters, value = stableclusters@meta.data$orig.ident)
DoHeatmap(stableclusters, features = ieg) + scale_fill_viridis() + theme(legend.position="top")

#subset for nuclei and plot ieg
nuclei <- subset(stableclusters, idents = 'nuclei')
DimPlot(nuclei, reduction = "tsne", features = ieg, label = FALSE, pt.size = 0.5) + NoLegend()
FeaturePlot(nuclei, features = ieg, reduction ="tsne")
nuclei <- SetIdent(nuclei, value = nuclei@meta.data$annotation)
x <- DotPlot(nuclei, features = ieg) + RotatedAxis()
x  + coord_flip() + theme(legend.position="top")


### subcluster central constellation ###

allcells <- readRDS("seurat_object_with_annotation.rds")
#Store annotation in metadata
allcells[["annotation"]] <- Idents(object = allcells)
#set active identity to seurat clusters
allcells <- SetIdent(allcells, value = allcells@meta.data$seurat_clusters)
#subset for blob clusters
new.cluster.ids <- c("central","central","central","peripheral","peripheral","peripheral","central","central","central","peripheral","peripheral","peripheral","central","central","central","peripheral","central","central","peripheral","peripheral","central","peripheral","central","peripheral","central","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","central","peripheral","peripheral","central","peripheral","peripheral","central","peripheral","central","central","peripheral","peripheral","peripheral","central","peripheral","peripheral","peripheral","central","central","peripheral","peripheral","peripheral","peripheral","peripheral","central","peripheral","peripheral","central","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","peripheral","central","peripheral","peripheral","peripheral")
names(new.cluster.ids) <- levels(allcells)
allcells <- RenameIdents(allcells, new.cluster.ids)
#Store annotation in metadata
allcells[["centralperipheral"]] <- Idents(object = allcells)

central <- subset(allcells, idents = 'central')
peripheral <- subset(allcells, idents = 'peripheral')

#split and reintegrate
central[["percent.mt"]] <- PercentageFeatureSet(central), pattern = "^mt-")
central.list <- SplitObject(central, split.by = "orig.ident")
central.list <- lapply(central.list, SCTransform, vars.to.regress = "percent.mt")
DefaultAssay(central) <- "SCT"
features  <- SelectIntegrationFeatures(central.list, nfeatures = 3000)
central.list <- PrepSCTIntegration(central.list, anchor.features = features)
anchors <- FindIntegrationAnchors(central.list, normalization.method = "SCT", anchor.features = features)
central_cor <- IntegrateData(anchorset = anchors, normalization.method = "SCT",k.weight = 46)
central_cor <- RunPCA(central_cor,verbose = FALSE, npcs = 300)
ElbowPlot(central_cor, ndims = 150, reduction = "pca")
central_cor <- RunTSNE(central_cor, dims = 1:150, verbose = FALSE)
central_cor <- RunUMAP(central_cor, reduction = "pca", dims = 1:150)
DefaultAssay(central_cor) <- "integrated"
central_cor <- FindNeighbors(central_cor, dims = 1:150)
central_cor <- FindClusters(central_cor, verbose = FALSE, resolution = 5)
DimPlot(central_cor, reduction = "umap", label = TRUE)
label <- DimPlot(central_cor, reduction = "tsne", label = TRUE)
label
DefaultAssay(central_cor) <- "RNA"
central.markers <- FindAllMarkers(central_cor, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
central.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top1 <- central.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
write.csv(central.markers,file = "central.markers.csv")
central_cor <- FindVariableFeatures(central_cor, selection.method = "vst", nfeatures = 1000)

#Dotplot of the 20 most variable genes and transcription factors for the central constellation
top20 <- c("LOC115211916","LOC115222309","LOC115218027","LOC115220921","LOC115220458","LOC118764179","LOC115211114","LOC115216569","LOC115212296","LOC118761014","LOC115218445","LOC115211681","LOC115216744","LOC115214332","LOC118761714","LOC115209927","LOC115232466","LOC115210874","LOC115209753","LOC115223269")

top20TF <- c("LOC115229429","LOC115209430","LOC115222465","LOC115216550","LOC115213424","LOC115227132","LOC115219954","LOC115213959","LOC115225998","LOC115212850","LOC115226973","LOC115214910","LOC115214070","LOC115231286","LOC115210541","LOC115226597","LOC115222049","LOC115209324","LOC115213203","LOC115215081")

x <- DotPlot(central_cor, features = top20) + RotatedAxis()
x + coord_flip()
x
y <- DotPlot(central_cor, features = top20TF) + RotatedAxis()
x + coord_flip()+ theme(legend.position = "top") + y + coord_flip() + theme(legend.position = "top")

DefaultAssay(central_cor) <- "RNA"

# Was the integration successful?
DimPlot(central_cor, reduction = "umap",group.by = "orig.ident")
splitbysample <- DimPlot(central_cor, reduction = "tsne",group.by = "orig.ident")

# Quality metrics for the central constellation
qc <- FeaturePlot(object = central_cor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), reduction ="tsne", ncol =3)
qc
splitbysample + label

IEG <- c("LOC115222465","LOC115214910","LOC115212294","LOC115228662")

#plot ieg on nuclei only
brain[["annotation"]] <- Idents(object = brain)
#Set orig.ident as active ident
brain <- SetIdent(brain, value = brain@meta.data$orig.ident)
#specify IEG
#subset for nuclei
nuclei <- subset(brain, idents = 'nuclei')
FeaturePlot(brain, reduction = "tsne", features = ieg, label = FALSE, pt.size = 0.5) + NoLegend()
FeaturePlot(nuclei, reduction = "tsne", features = ieg, label = FALSE, pt.size = 0.5) + NoLegend()
plot <- FeaturePlot(nuclei, features = c("LOC115228662"), reduction ="tsne")


#Compare quality control metrics central versus periphery
allcells <- readRDS("seurat_object_with_annotation.rds")
splitbysample <- DimPlot(allcells, reduction = "tsne",group.by = "centralperipheral",  cols=c("#E52320", "#FBC381")) +theme(legend.position="top",legend.justification="left")
splitbysample 
qc <- FeaturePlot(object = allcells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), split.by = "centralperipheral", reduction ="tsne", ncol =3)
qc
VlnPlot(allcells, features =  c("nFeature_RNA", "nCount_RNA", "percent.mt"), n = 3, pt.size = 0, cols=c("#E52320", "#FBC381"), split.by ="orig.ident")
vln <- VlnPlot(allcells, features =  c("nFeature_RNA", "nCount_RNA", "percent.mt"), n = 2, pt.size = 0, cols=c("#E52320", "#FBC381"))
vln2 <- VlnPlot(nuclei, features =  c("LOC115228662"), pt.size = 0, cols=c("#E52320", "#FBC381")) + NoLegend()
vln2
violinplots <- vln + vln2
nuclei<- SetIdent(allcells, value = allcells@meta.data$orig.ident)
nucleionly <- SetIdent(allcells, value = allcells@meta.data$orig.ident)
nuclei <- subset(nucleionly, idents = 'nuclei')
nuclei<- SetIdent(nuclei, value = nuclei@meta.data$centralperipheral)
splitbysample + violinplots

#top ten differentially expressed genes between peripheral and central brain

centralperiphery.markers <- FindAllMarkers(allcells, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
centralperiphery.markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)
write.csv(centralperiphery.markers,file = "centralperiphery.markers.csv")
vlnTF <- VlnPlot(allcells, features = c("LOC115232381", "LOC115209324", "LOC115212946","LOC115217234"), n = 4, pt.size = 0, cols=c("#E52320", "#FBC381"))
vlnTF

### Neurotransmitters and neuropeptides ###
#Fig 2
#Subset for NT/NP
ntnp <- c("LOC115219340","LOC115215891","LOC115220799","LOC115213095","LOC115215674","LOC115230059","LOC115215438","LOC115220451","LOC115212373","LOC115229488","LOC115215366","LOC115216569","LOC115217507","LOC115218107","LOC115224553","LOC115209816","LOC115220458","LOC115218445","LOC115229259","LOC115215763","LOC115214332","LOC115212296","LOC115219332","LOC115231438")
NTNP <- subset(brain, features = ntnp)

#Order heatmap by Avg_logFC
x <- DotPlot(stable_clusters, features =NP, cluster.idents = TRUE)
x + coord_flip() + theme(axis.text.x = element_text(angle = 90))  


### Differential expression between glial subtypes ###
glia1.markers <- FindMarkers(stableclusters, ident.1 = "GLIA1", ident.2 = c("GLIA2","GLIA3"), min.pct = 0.50, only.pos = TRUE,logfc.threshold = 0.25)
head(glia1.markers, n = 5)
glia2.markers <- FindMarkers(stableclusters, ident.1 = "GLIA2", ident.2 = c("GLIA1","GLIA3"), min.pct = 0.50, only.pos = TRUE,logfc.threshold = 0.25)
head(glia2.markers, n = 5)
glia3.markers <- FindMarkers(stableclusters, ident.1 = "GLIA3", ident.2 = c("GLIA1","GLIA2"), min.pct = 0.50, only.pos = TRUE,logfc.threshold = 0.25)
head(glia3.markers, n = 5)

glia <- subset(stableclusters, idents = c("GLIA1","GLIA2","GLIA3"))

#dotplot for glial subtypes

glialmarkers <- c("LOC115214230","LOC115221317","LOC115215034","LOC115231765","LOC115210633","LOC115214418","LOC115214419","LOC115215814","LOC115217464","LOC115223080","LOC115212306","LOC115230501","LOC115229261","LOC115210961","LOC115218235")

x <- DotPlot(glia, features = glialmarkers) + RotatedAxis()
x + coord_flip() + theme(legend.position="top")
x
library(tidyverse)
DoHeatmap(glia, features = features, size = 3)


### protocadherins ###

#Produce dotplot with clustered features to plot all pcdh 
#subset for all pcdh
allpcdh <- c("LOC115219230","LOC118766122","LOC115219489","LOC115218755","LOC115218963","LOC115219524","LOC115219527","LOC115218752","LOC115218736","LOC115218968","LOC118766056","LOC115219539","LOC115218738","LOC115218739","LOC115218965","LOC115218966","LOC115219531","LOC115218947","LOC115219532","LOC115218756","LOC115219535","LOC115218740","LOC115218743","LOC115218750","LOC115218971","LOC115219480","LOC115219425","LOC115218744","LOC115218969","LOC115218949","LOC115218953","LOC115219540","LOC115218758","LOC115219434","LOC118765869","LOC115218977","LOC115219541","LOC115219402","LOC115219116","LOC115219391","LOC115219392","LOC115218999","LOC115218763","LOC115219395","LOC115218862","LOC115219408","LOC118766074","LOC115219398","LOC115219493","LOC115219494","LOC115219113","LOC118766075","LOC115219412","LOC115218959","LOC115219393","LOC115219401","LOC118766076","LOC115218957","LOC115219231","LOC115219399","LOC115219407","LOC118766087","LOC118766090","LOC118766121","LOC115219382","LOC115219429","LOC115218950","LOC115219114","LOC118766091","LOC115219112","LOC115218960","LOC115219368","LOC115219385","LOC115219388","LOC118766079","LOC118766080","LOC118766081","LOC118766082","LOC115218983","LOC115219232","LOC115219409","LOC115218954","LOC115219387","LOC115219233","LOC118766083","LOC115219117","LOC115219427","LOC118766125","LOC115219534","LOC115218988","LOC115219528","LOC115219370","LOC115218981","LOC115219236","LOC115218987","LOC115219237","LOC115219543","LOC115218982","LOC115219537","LOC115219238","LOC118766124","LOC115218955","LOC115218762","LOC115218741","LOC118766084","LOC115219216","LOC115218748","LOC118766123","LOC115219377","LOC115219428","LOC115218986","LOC115219542","LOC115218737","LOC115218780","LOC115219115","LOC115219371","LOC118766018","LOC115218747","LOC115218749","LOC115219365","LOC115219378","LOC115219529","LOC115219374","LOC115218781","LOC115218757","LOC115219405","LOC115219538","LOC115219183","LOC115219404","LOC115219544","LOC115218746","LOC115219375","LOC115219426","LOC118766077","LOC115219418","LOC115219474","LOC118760832","LOC115219410","LOC118766078","LOC115218958","LOC118766021","LOC115219390","LOC118766045","LOC115219394","LOC115219396","LOC115219413","LOC115219386","LOC115218978","LOC115219373","LOC115219120","LOC115218967","LOC115218975","LOC115219424","LOC118766093","LOC115218938","LOC115218956","LOC115219433","LOC118766019","LOC115219421","LOC115218974","LOC115219369","LOC118766058","LOC115219406","LOC115219545","LOC115221427","LOC115228361","LOC115228506","LOC115228963","LOC115229210","LOC115229310","LOC115229484","LOC115230295","LOC115230600","LOC115230633","LOC115230635","LOC115230634","LOC115230882","LOC118761750","LOC118761751","LOC118761753","LOC115231635","LOC115231636","LOC115232042","LOC115232041","LOC115231640","LOC115218437","LOC115218349")
PCDH_subset <- subset(stableclusters, features = allpcdh)
DefaultAssay(PCDH_subset) <- "RNA"
PCDH.cluster.averages <- AverageExpression(PCDH_subset)
PCDH_RNA <- PCDH.cluster.averages[["RNA"]]
write.csv(PCDH_RNA, "PCDH_RNA.csv")

#Order heatmap
all_pcdh_sort <- c("LOC115218968","LOC115219528","LOC118766058","LOC115218965","LOC115219541","LOC118766045","LOC115218750","LOC115218966","LOC118766074","LOC115218969","LOC115218744","LOC115219238","LOC115230635","LOC115219231","LOC115230633","LOC118761751","LOC115218739","LOC115219489","LOC115219406","LOC118765869","LOC115218977","LOC115219391","LOC115219392","LOC118766124","LOC118766122","LOC115219120","LOC115219405","LOC115219404","LOC115218978","LOC115219230","LOC115218757","LOC115219114","LOC115218756","LOC115219421","LOC115218762","LOC115218862","LOC115219409","LOC115218736","LOC115218741","LOC115218740","LOC115219116","LOC115219236","LOC115218971","LOC115218737","LOC115218743","LOC115219369","LOC115219115","LOC115219544","LOC115218955","LOC115219535","LOC118766125","LOC115218956","LOC115219494","LOC115219387","LOC118766078","LOC115219407","LOC115218738","LOC118761753","LOC115219429","LOC115219413","LOC115218781","LOC118766075","LOC115218960","LOC115218947","LOC115219537","LOC115219393","LOC115219396","LOC115219368","LOC115218950","LOC115218938","LOC118760832","LOC115219538","LOC115219534","LOC115218949","LOC118766082","LOC115219388","LOC115219112","LOC118766121","LOC115219424","LOC115219527","LOC115219539","LOC118766077","LOC115219386","LOC115219237","LOC115219474","LOC115219398","LOC115232041","LOC115219532","LOC115219117","LOC115219433","LOC115218974","LOC115219410","LOC118766080","LOC115219412","LOC118766081","LOC118766083","LOC115219385","LOC115219480","LOC115219399","LOC115218747","LOC115219390","LOC118766076","LOC115219426","LOC115230882","LOC118766087","LOC115219402","LOC115219427","LOC115219543","LOC115229210","LOC118766123","LOC115219418","LOC115218958","LOC115218986","LOC115219232","LOC115219216","LOC115219434","LOC115218983","LOC115219183","LOC115219395","LOC115219394","LOC115218746","LOC115219233","LOC118766084","LOC115218749","LOC115218748","LOC115230600","LOC115218763","LOC115219113","LOC115219408","LOC115230634","LOC115218975","LOC115219542","LOC115219529","LOC115218959","LOC115219401","LOC115219428","LOC118766091","LOC115218981","LOC118766093","LOC115218349","LOC115218982","LOC118766079","LOC115219382","LOC115219375","LOC115218957","LOC118766090","LOC115219374","LOC115218758","LOC115218988","LOC115218437","LOC115221427","LOC115219365","LOC118766018","LOC115218954","LOC115219545","LOC115218953","LOC115219540","LOC118766019","LOC115219531","LOC115219377","LOC115219378","LOC115218987","LOC115218755","LOC115219373","LOC115218780","LOC115219370","LOC115219371","LOC115219425","LOC115232042","LOC115218999","LOC115219493","LOC118766021")

#produce dotplot plot in different colors
octopus_colors <- c('black','#29B5BC','#ED756F')

pdf(file = "PCDH_dotplot_clustered_expression_colors.pdf",
    width = 25, # The width of the plot in inches
    height = 10)


for (clrs in 1:length(octopus_colors)){
  
  dots <- DotPlot(stableclusters, features = all_pcdh_sort, cluster.idents = TRUE, cols = c("white", octopus_colors[clrs]))  + RotatedAxis() + theme(axis.text.x = element_text(size = 12)) + theme(axis.text.y = element_text(size = 12)) + theme(legend.title = element_text(size = 10)) + theme(legend.text = element_text(size = 7)) +theme(legend.position = "top")
   print(dots)
}
print(dots)
dev.off()

#Protocadherin expression per cell type and per cell
#read in data
stableclusters <- readRDS("stable_clusters.rds")
#install packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(dplyr)

#Set output directory
setwd("~/pcdh dotplot")
#read in data to make stacked barplot
pcdh_per_celltype <- read.csv("~/pcdh dotplot/pcdh_per_celltype.csv", sep=";")
#basic plot
p <- ggplot(df, aes(x = Cell.type, y = X.PCDH))+
  geom_col(aes(fill = type), width = 0.7)
p

head(pcdh_per_celltype)

df <- pcdh_per_celltype %>%
  group_by(type) %>%
  arrange(type,"clustered","not clustered")


p2 <- ggplot(pcdh_per_celltype, aes(fill=type, y=X.PCDH, x=Cell.type)) + geom_bar(position="stack", stat="identity") + RotatedAxis()
p3 <- p2 + scale_y_continuous(expand = expansion(mult = c(0, 0))) + labs(x= "Cell type", y = "Number of PCDHs") + theme(legend.position = "top")
p3 + scale_fill_manual(values = c("#29B5BC", "#ED756F"))

#save plot as a pdf landscape width 25, height 10

#produce the same plot but with the sum of the average pcdh expression per cell type

sum_average_expression_pcdh_per_celltype <- read.csv2("~/pcdh dotplot/sum_average_expression_pcdh_per_celltype.csv")

p2 <- ggplot(sum_average_expression_pcdh_per_celltype, aes(fill=type, y=X.PCDH, x=Cell.type)) + geom_bar(position="stack", stat="identity") + RotatedAxis()
p3 <- p2 + scale_y_continuous(expand = expansion(mult = c(0, 0))) + labs(x= "Cell type", y = "Sum of average PCDH expression") + theme(legend.position = "top")
p3 + scale_fill_manual(values = c("#29B5BC", "#ED756F"))
#save plot as a pdf landscape width 25, height 10

#read in number of PCDHs/cell 
number.of.pcdh.per.cell <- read.table("~/pcdh dotplot/number of pcdh per cell.csv", quote="\"", comment.char="")

plot_allstablecells <- ggplot(number.of.pcdh.per.cell, aes(V1)) + geom_histogram(binwidth=2, fill="#575959", col="grey")  + scale_y_continuous(expand = expansion(mult = c(0, 0)))  + scale_x_continuous(expand = expansion(mult = c(0, 0)))  + labs(x= "Number of PCDHs", y = "Number of Cells") 
p1 <-plot_allstablecells %+% geom_vline(aes(xintercept = mean(V1)),col='#29B5BC',size=1)
p1

#subset stable clusters for neuronal versus non neuronal cell types
allpcdh <- c("LOC115219230","LOC118766122","LOC115219489","LOC115218755","LOC115218963","LOC115219524","LOC115219527","LOC115218752","LOC115218736","LOC115218968","LOC118766056","LOC115219539","LOC115218738","LOC115218739","LOC115218965","LOC115218966","LOC115219531","LOC115218947","LOC115219532","LOC115218756","LOC115219535","LOC115218740","LOC115218743","LOC115218750","LOC115218971","LOC115219480","LOC115219425","LOC115218744","LOC115218969","LOC115218949","LOC115218953","LOC115219540","LOC115218758","LOC115219434","LOC118765869","LOC115218977","LOC115219541","LOC115219402","LOC115219116","LOC115219391","LOC115219392","LOC115218999","LOC115218763","LOC115219395","LOC115218862","LOC115219408","LOC118766074","LOC115219398","LOC115219493","LOC115219494","LOC115219113","LOC118766075","LOC115219412","LOC115218959","LOC115219393","LOC115219401","LOC118766076","LOC115218957","LOC115219231","LOC115219399","LOC115219407","LOC118766087","LOC118766090","LOC118766121","LOC115219382","LOC115219429","LOC115218950","LOC115219114","LOC118766091","LOC115219112","LOC115218960","LOC115219368","LOC115219385","LOC115219388","LOC118766079","LOC118766080","LOC118766081","LOC118766082","LOC115218983","LOC115219232","LOC115219409","LOC115218954","LOC115219387","LOC115219233","LOC118766083","LOC115219117","LOC115219427","LOC118766125","LOC115219534","LOC115218988","LOC115219528","LOC115219370","LOC115218981","LOC115219236","LOC115218987","LOC115219237","LOC115219543","LOC115218982","LOC115219537","LOC115219238","LOC118766124","LOC115218955","LOC115218762","LOC115218741","LOC118766084","LOC115219216","LOC115218748","LOC118766123","LOC115219377","LOC115219428","LOC115218986","LOC115219542","LOC115218737","LOC115218780","LOC115219115","LOC115219371","LOC118766018","LOC115218747","LOC115218749","LOC115219365","LOC115219378","LOC115219529","LOC115219374","LOC115218781","LOC115218757","LOC115219405","LOC115219538","LOC115219183","LOC115219404","LOC115219544","LOC115218746","LOC115219375","LOC115219426","LOC118766077","LOC115219418","LOC115219474","LOC118760832","LOC115219410","LOC118766078","LOC115218958","LOC118766021","LOC115219390","LOC118766045","LOC115219394","LOC115219396","LOC115219413","LOC115219386","LOC115218978","LOC115219373","LOC115219120","LOC115218967","LOC115218975","LOC115219424","LOC118766093","LOC115218938","LOC115218956","LOC115219433","LOC118766019","LOC115219421","LOC115218974","LOC115219369","LOC118766058","LOC115219406","LOC115219545","LOC115221427","LOC115228361","LOC115228506","LOC115228963","LOC115229210","LOC115229310","LOC115229484","LOC115230295","LOC115230600","LOC115230633","LOC115230635","LOC115230634","LOC115230882","LOC118761750","LOC118761751","LOC118761753","LOC115231635","LOC115231636","LOC115232042","LOC115232041","LOC115231640","LOC115218437","LOC115218349")

#subset for non neuronal
nonneuronal <- subset(x = stableclusters, idents = c("GLIA1", "GLIA2","GLIA3", "TBA8", "HC", "EC", "FBL", "PREC"))
#subset for pcdh
PCDH_noneuronal <- subset(nonneuronal, features = allpcdh)
#subset for neuronal clusters
neuronal <- subset(x = stableclusters, idents = c("GLIA1", "GLIA2","GLIA3", "TBA8", "HC", "EC", "FBL", "PREC"), invert = TRUE)
#subset for pcdh
PCDH_neuronal <- subset(neuronal, features = allpcdh)
#save raw counts as csv and open in excel
write.csv(PCDH_neuronal@assays[["RNA"]]@counts, "PCDH_neuronal_RNA.csv")
write.csv(PCDH_noneuronal@assays[["RNA"]]@counts, "PCDH_nonneuronal_RNA.csv")

#read in number of PCDHs/cell 
number.of.pcdh.per.neuronal.cell <- read.table("~/temporary files R studio/pcdh dotplot/number of pcdh per neuronal cell.csv", quote="\"", comment.char="")

plot_allneurons <- ggplot(number.of.pcdh.per.neuronal.cell, aes(V1)) + geom_histogram(binwidth=2, fill="#575959", col="grey")  + scale_y_continuous(expand = expansion(mult = c(0, 0)))  + scale_x_continuous(expand = expansion(mult = c(0, 0)))  + labs(x= "Number of PCDHs", y = "Number of Neurons")
p2 <- plot_allneurons %+% geom_vline(aes(xintercept = mean(V1)),col='#29B5BC',size=1)
p2

number.of.pcdh.per.non.neuronal.cell <- read.table("~/temporary files R studio/pcdh dotplot/number of pcdh per non neuronal cell.csv", quote="\"", comment.char="")
plot_allnonneurons <- ggplot(number.of.pcdh.per.non.neuronal.cell, aes(V1)) + geom_histogram(binwidth=2, fill="#575959", col="grey")  + scale_y_continuous(expand = expansion(mult = c(0, 0)))  + scale_x_continuous(expand = expansion(mult = c(0, 0)))  + labs(x= "Number of PCDHs", y = "Number of non-neuronal cells")
p3 <- plot_allnonneurons %+% geom_vline(aes(xintercept = mean(V1)),col='#29B5BC',size=1)
p3

p1+p2+p3

#save this figure in the same dimensions; a pdf landscape width 25, height 10


### Script to calculate the % expression of a gene ###
#https://github.com/satijalab/seurat/issues/371

PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}
#Elav, Onecut
PrctCellExpringGene(allcells ,genes =c("LOC115215662","LOC115220621"))
#Gs2, Apolpp
PrctCellExpringGene(stable_clusters ,genes =c("LOC115222946","LOC115230720"))
#Vglut, Vacht, Th, DAT
PrctCellExpringGene(brain ,genes =c("LOC115220451","LOC115212373", "LOC115215438", "LOC115230059"))

#Heatmaps for PCDHs, GPCRs, ZnF
#PCDHs
PCDHs <- c(Pcdhs$V1)

stable_clusters <- brain
HVGPCDH <- c("LOC115218437","LOC115218736","LOC115218737","LOC115218738","LOC115218739","LOC115218740","LOC115218741","LOC115218743","LOC115218746","LOC115218747","LOC115218748","LOC115218749","LOC115218750","LOC115218755","LOC115218757","LOC115218758","LOC115218762","LOC115218763","LOC115218862","LOC115218947","LOC115218949","LOC115218950","LOC115218955","LOC115218969","LOC115218971","LOC115218974","LOC115218975","LOC115218977","LOC115218978","LOC115218982","LOC115218983","LOC115218986","LOC115218987","LOC115218988","LOC115219112","LOC115219115","LOC115219116","LOC115219183","LOC115219216","LOC115219230","LOC115219231","LOC115219232","LOC115219233","LOC115219236","LOC115219365","LOC115219368","LOC115219370","LOC115219371","LOC115219377","LOC115219404","LOC115219405","LOC115219474","LOC115219480","LOC115219489","LOC115219527","LOC115219528","LOC115219529","LOC115219531","LOC115219532","LOC115219534","LOC115219535","LOC115219537","LOC115219538","LOC115219539","LOC115219540","LOC115219541","LOC115219542","LOC115219543","LOC115219544","LOC115219545","LOC115229210","LOC118760832","LOC118765869","LOC118766075","LOC118766080")
DefaultAssay(stable_clusters) <- "SCT"
PCDHs_S <- subset(stable_clusters, features = PCDHs)
DefaultAssay(PCDHs_S) <- "RNA"
PCDH.cluster.averages <- AverageExpression(PCDHs)
PCDH_RNA <- PCDH.cluster.averages[["integrated"]]
write.csv(PCDH_RNA, "PCDH_RNA.csv")

write.csv(PCDH_subset@assays[["RNA"]]@counts, "PCDH_RNA_counts_average_exp.csv")


library(ComplexHeatmap)
my_matrix <- as.matrix(PCDH_RNA)
scaled_mat = t(scale(t(my_matrix)))
library(viridis)
my_pallete <- viridis(256, option = "D")
fontsize <- 1

PCDH.heatmap <- Heatmap(scaled_mat, cluster_rows = TRUE, 
        cluster_columns= TRUE,
        column_names_side = "top", 
        row_names_side = "left", 
        column_names_gp=gpar(cex=1,1), 
        row_names_gp=gpar(cex=0.7),
        show_row_names = TRUE, show_column_names = TRUE, name ="Avg Exp",col = my_pallete,  column_title = "Cell types", 
        row_title = "PCDHs", show_column_dend = TRUE, show_row_dend = FALSE)

PCDH.heatmap

#Znf
HVGZnF <- c("LOC115208765","LOC115209070","LOC115209136","LOC115209253","LOC115209324","LOC115209537","LOC115209715","LOC115209797","LOC115210091","LOC115210212","LOC115210281","LOC115210289","LOC115210709","LOC115211753","LOC115212031","LOC115212047","LOC115212072","LOC115212762","LOC115212984","LOC115213034","LOC115213193","LOC115214010","LOC115214070","LOC115214133","LOC115214159","LOC115214230","LOC115214501","LOC115214759","LOC115214896","LOC115214953","LOC115214961","LOC115215412","LOC115215425","LOC115215499","LOC115215516","LOC115215583","LOC115215677","LOC115215688","LOC115216354","LOC115216394","LOC115216631","LOC115216706","LOC115217000","LOC115217040","LOC115217077","LOC115217234","LOC115217255","LOC115217322","LOC115217707","LOC115218143","LOC115218350","LOC115218768","LOC115218805","LOC115218818","LOC115218821","LOC115218822","LOC115218828","LOC115218829","LOC115218832","LOC115218838","LOC115218843","LOC115218852","LOC115218854","LOC115218855","LOC115218874","LOC115218876","LOC115218900","LOC115218903","LOC115218904","LOC115218910","LOC115218921","LOC115218933","LOC115219001","LOC115219031","LOC115219047","LOC115219052","LOC115219107","LOC115219135","LOC115219187","LOC115219189","LOC115219273","LOC115219325","LOC115219331","LOC115219332","LOC115219345","LOC115219354","LOC115219445","LOC115219487","LOC115219492","LOC115219669","LOC115219871","LOC115220292","LOC115220746","LOC115220963","LOC115221075","LOC115221120","LOC115221329","LOC115221342","LOC115221376","LOC115221420","LOC115221476","LOC115221666","LOC115221708","LOC115221776","LOC115221811","LOC115222524","LOC115222603","LOC115222733","LOC115223641","LOC115224245","LOC115224323","LOC115224503","LOC115224742","LOC115224744","LOC115224869","LOC115225497","LOC115225746","LOC115225843","LOC115225892","LOC115226337","LOC115226572","LOC115226582","LOC115227205","LOC115227486","LOC115228118","LOC115228662","LOC115229099","LOC115229234","LOC115232434","LOC115232482","LOC115232562","LOC115232597","LOC115232598","LOC118762166","LOC118762983","LOC118765768","LOC118765937","LOC118765939","LOC118767207","LOC118767500","LOC118768437")  
DefaultAssay(stable_clusters) <- "SCT"
ZnFs <- subset(stable_clusters, features = Znfs)
DefaultAssay(ZnFs) <- "RNA"
ZnF.cluster.averages <- AverageExpression(ZnFs)
ZnF_SCT <- ZnF.cluster.averages[["integrated"]]
write.csv(ZnF_SCT, "ZnF_SCT.csv")


library(ComplexHeatmap)
my_matrix <- as.matrix(ZnF_SCT)
scaled_mat = t(scale(t(my_matrix)))
library(viridis)
my_pallete <- viridis(256, option = "D")
fontsize <- 1

ZnF.heatmap <- Heatmap(scaled_mat, cluster_rows = TRUE, 
                        cluster_columns= TRUE,
                        column_names_side = "top", 
                        row_names_side = "left", 
                        column_names_gp=gpar(cex=1,1), 
                        row_names_gp=gpar(cex=0.7),
                        show_row_names = TRUE, show_column_names = TRUE, name ="Avg Exp",col = my_pallete,  column_title = "Cell types", 
                        row_title = "ZnFs", show_column_dend = TRUE, show_row_dend = FALSE)
ZnF.heatmap

#GPCRs
GPCRs <- c(GPCRs$V1)
HVGGPCR <- c("LOC115209262","LOC115209264","LOC115209347","LOC115209555","LOC115209684","LOC115209755","LOC115209874","LOC115210334","LOC115210453","LOC115210524","LOC115210538","LOC115210745","LOC115210951","LOC115211759","LOC115212036","LOC115212179","LOC115212306","LOC115213029","LOC115213121","LOC115213189","LOC115213266","LOC115213353","LOC115213422","LOC115213470","LOC115213477","LOC115213895","LOC115214184","LOC115214417","LOC115214436","LOC115214450","LOC115214515","LOC115214785","LOC115214856","LOC115214916","LOC115214943","LOC115214950","LOC115215078","LOC115215151","LOC115215156","LOC115215171","LOC115215592","LOC115215745","LOC115215842","LOC115215845","LOC115215927","LOC115216181","LOC115216203","LOC115216255","LOC115216411","LOC115216420","LOC115216428","LOC115216455","LOC115216514","LOC115216603","LOC115216703","LOC115216729","LOC115216781","LOC115216961","LOC115216988","LOC115217118","LOC115217159","LOC115217305","LOC115217464","LOC115217869","LOC115217894","LOC115217931","LOC115217992","LOC115218008","LOC115218025","LOC115218041","LOC115218067","LOC115218068","LOC115218117","LOC115218246","LOC115218300","LOC115218942","LOC115219158","LOC115219180","LOC115219201","LOC115219274","LOC115219515","LOC115219701","LOC115219741","LOC115219757","LOC115219772","LOC115219774","LOC115219911","LOC115220372","LOC115220511","LOC115220541","LOC115220560","LOC115220758","LOC115220891","LOC115220974","LOC115221220","LOC115222127","LOC115222207","LOC115222286","LOC115223004","LOC115223210","LOC115223272","LOC115223637","LOC115224446","LOC115224729","LOC115224846","LOC115224862","LOC115225199","LOC115225647","LOC115225832","LOC115225961","LOC115226204","LOC115226896","LOC115227278","LOC115227327","LOC115227336","LOC115227615","LOC115230791","LOC115230855","LOC115231007","LOC115231120","LOC115231288","LOC115231298","LOC115231777","LOC115232411","LOC115232543","LOC115232565","LOC115232578","LOC115232618","LOC115232653","LOC118767411")
DefaultAssay(stable_clusters) <- "SCT"
GPCRs <- subset(stable_clusters, features = GPCRs)
DefaultAssay(GPCRs) <- "RNA"
GPCRs.cluster.averages <- AverageExpression(GPCRs)
GPCRs_SCT <- GPCRs.cluster.averages[["integrated"]]
write.csv(GPCRs_SCT, "GPCR_SCT.csv")
GPCRs_n <- SCTransform(GPCRs, verbose = FALSE)

library(ComplexHeatmap)
my_matrix <- as.matrix(GPCRs_SCT)
scaled_mat = t(scale(t(my_matrix)))
library(viridis)
my_pallete <- viridis(256, option = "D")
fontsize <- 1

GPCRs.heatmap <- Heatmap(scaled_mat, cluster_rows = TRUE, 
                       cluster_columns= TRUE,
                       column_names_side = "top", 
                       row_names_side = "left", 
                       column_names_gp=gpar(cex=1,1), 
                       row_names_gp=gpar(cex=0.7),
                       show_row_names = TRUE, show_column_names = TRUE, name ="Avg Exp",col = my_pallete,  column_title = "Cell types", 
                       row_title = "GPCRs", show_column_dend = TRUE, show_row_dend = FALSE)
GPCRs.heatmap

#Prepare figure on cell cycle
DefaultAssay(brain) <- "RNA"
brain <- RenameIdents(object = brain, 'PRECL'="TBA8")
library(tidyverse)
library(ggplot2)
sphase <- c("LOC115222129","LOC115230618","LOC115209366","LOC115226207","LOC115223095","LOC115210989","LOC115219310","LOC115223264","LOC115214153","LOC115210734","LOC115210722","LOC115232599","LOC115217460","LOC115225861","LOC115217680")
Sg2m <- c("LOC115223887","LOC115213349","LOC115212935","LOC115223276","LOC115225396","LOC115228146","LOC115221018","LOC115221466","LOC115216795","LOC115219111","LOC115223993","LOC115224806","LOC115218498","LOC115210279","LOC115216844","LOC115221619","LOC115230618","LOC115209366","LOC115226207","LOC115223095","LOC115210989","LOC115219310","LOC115223264","LOC115214153","LOC115210734","LOC115210722","LOC115232599","LOC115217460","LOC115225861","LOC115217680","LOC115222129")
x <- DotPlot(brain, features = sphase, cluster.idents =TRUE)
x + coord_flip() + theme(axis.text.x = element_text(angle = 90))  
y <- DotPlot(brain, features = Sg2m, cluster.idents =TRUE)
y + coord_flip() + theme(axis.text.x = element_text(angle = 90))  

#all homeobox tf dotplot 
HD <- c("LOC118763461","LOC115210569","LOC115219954","LOC115222049","LOC115223154","LOC115210984","LOC115219512","LOC115216272","LOC115219178","LOC115209922","LOC115216303","LOC115214412","LOC115214327","LOC118765224","LOC115213203","LOC115212850","LOC115219234","LOC115213959","LOC115227132","LOC115231385","LOC115221231","LOC115230438","LOC115218560","LOC115230509","LOC115214240","LOC115230771","LOC115227693","LOC115226966","LOC115226973","LOC115209827","LOC115226048","LOC115210541","LOC118767514","LOC118766466","LOC115215778","LOC115227652","LOC115231183","LOC115226895","LOC115227008","LOC115209485","LOC115225998","LOC115221723","LOC115209521","LOC118767065","LOC115215153","LOC115231384","LOC115227071","LOC115215102","LOC115231388","LOC115220811","LOC115210669","LOC115211005","LOC115218558")
DefaultAssay(brain) <- "RNA"
a <- DotPlot(brain, features = HD, cluster.idents =TRUE)
a + coord_flip() + theme(axis.text.x = element_text(angle = 90))  

#Dotplot ISH genes
ish <- c("LOC115215438","LOC115212373","LOC115220451","LOC115215891","LOC115220799","LOC115231765","LOC115215674","LOC115218755","LOC115219340","LOC115218977","LOC115215570","LOC115209882","LOC115213270","LOC115215366","LOC118766118","LOC115227546","LOC115220878","LOC115224224","LOC115220621","LOC115218107","LOC115222049","LOC115212290","LOC115209866","LOC115230720","LOC115222946","LOC115211706","LOC115211717","LOC115218638","LOC115212850","LOC118765224")
x <- DotPlot(stable_clusters, features = ish, cluster.idents =TRUE)
x  + coord_flip() + theme(axis.text.x = element_text(angle = 90)) + theme(legend.position="top") + scale_colour_gradient2(low = "grey", mid = "grey", high = "#00BFC4")


