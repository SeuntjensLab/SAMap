#### Code used for cross-species analysis ####

### Plot SAMap annotation on tSNE
brain_most_stable_clusters <- readRDS("brain_most_stable_clusters.rds")

#Plot fly annotation
fly <- c("Not mapped","Not mapped","Not mapped","T1 neurons (0,31)","ensheathing glia (0,54)","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","T1 neurons (0,31)","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Dopaminergic PAM neurons (0,71)","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Gamma-Kenyon cells (0,27)","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","C3 neurons (0,27)","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Alpha/beta-Kenyon cells (0,25)","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Plasmocytes (0,69)","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped","Not mapped")
names(fly) <- levels(brain_most_stable_clusters)
brain_most_stable_clusters_fly <- RenameIdents(brain_most_stable_clusters, fly)
saveRDS(brain_most_stable_clusters_fly, file = "brain_most_stable_clusters_fly.rds")

#mouse
mouse <- c("not mapped","not mapped","not mapped","not mapped","Non-telencephalon astrocytes    Hypothalamus, Thalamus (0,68)","not mapped","not mapped","not mapped","not mapped","not mapped","Endothelial cells (0,27)","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","Dopaminergic interneurons (0,27)","not mapped","not mapped","not mapped","not mapped","not mapped","Oligodendrocytes (0,80)","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","Microglia (0,40)","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","Telencephalon astrocytes    Cortex, Hippocampus, Striatum, Amygdala (0,32)","not mapped","not mapped","not mapped","not mapped")
names(mouse) <- levels(brain_most_stable_clusters)
brain_most_stable_clusters_mouse <- RenameIdents(brain_most_stable_clusters, mouse)
saveRDS(brain_most_stable_clusters_mouse, file = "brain_most_stable_clusters_mouse.rds")

#read in same colours as for the original tsne
brain_most_stable_clusters <- readRDS("brain_most_stable_clusters.rds")
annotation <- c("not separated","unstable","unstable","IGL2- GLUT/DOP","GLIA1","OGL2-DOP","unstable","unstable","not separated","DOP2","FBL","SUB","not separated","unstable","unstable","IGL2- GLUT/DOP","unstable","not separated","unstable","OGL1","unstable","GLUT1","unstable","unstable","unstable","DOP3","SERT","GLUT3","DOP1","IGL3","TBA1","unstable","PROG","TBA2","unstable","TBA3","unstable","unstable","VL","unstable","unstable","unstable","unstable","GABA","unstable","unstable","IGL1-OA","unstable","not separated","unstable","GLUT2","OA","TBA4","ACH","EC","unstable","PEP2","GLIA2","not separated","GLUT4","TBA5","OGL3-OA","CCAP","unstable","TBA6","unstable","TBA7","TBA8","unstable","unstable","HC","unstable","TBA9","unstable","TBA10","TBA11","PEP1","TBA12","unstable","unstable","unstable","unstable","TBA13","unstable","unstable","TBA14","unstable")
names(annotation) <- levels(brain_most_stable_clusters)
brain_most_stable_clusters_annotated <- RenameIdents(brain_most_stable_clusters, annotation)
saveRDS(brain_most_stable_clusters_annotated, file = "brain_most_stable_clusters_annotated.rds")
brain_most_stable_clusters<- readRDS("brain_most_stable_clusters_annotated.rds")
memory.limit(size = 15000)
p <- DimPlot(brain_most_stable_clusters_annotated, reduction = 'tsne', pt.size = 1, label=FALSE)
p
library(ggplot2)
pbuild <- ggplot_build(p)
pdata <- pbuild$data[[1]] 
pdata <-  pdata[order(pdata$group), ] # Order the plot data by group
ucols <- unique(pdata$colour) # Get a vector of unique colors
names(ucols) <- unique(pdata$group) # Add the groups to the vector of colors as names
ucols
brain_most_stable_clusters_annotated <- readRDS("brain_most_stable_clusters_annotated.rds")
octopus <- DimPlot(brain_most_stable_clusters_annotated, reduction = 'tsne', pt.size = 1, label=FALSE, cols = c('not separated' = '#aeadb3','unstable' = '#aeadb3','IGL2- GLUT/DOP' = '#EC8239','GLIA1' = '#E48800','OGL2-DOP' = '#DB8E00','DOP2' = '#D29300','FBL' = '#C79800','SUB' = '#BB9D00','OGL1' = '#AEA200','GLUT1' = '#A0A600','DOP3' = '#8FAA00','SERT' = '#7CAE00','GLUT3' = '#64B200','DOP1' = '#43B500','IGL3' = '#00B81B','TBA1' = '#00BA43','PROG' = '#00BD5C','TBA2' = '#00BE72','TBA3' = '#00C085','VL' = '#00C096','GABA' = '#00C1A7','IGL1-OA' = '#00C0B6','GLUT2' = '#00BFC4','OA' = '#00BDD2','TBA4' = '#00BADE','ACH' = '#00B7E9','EC' = '#00B2F3','PEP2' = '#00ACFB','GLIA2' = '#00A6FF','GLUT4' = '#509FFF','TBA5' = '#7C96FF','OGL3-OA' = '#9A8EFF','CCAP' = '#B385FF','TBA6' = '#C77CFF','TBA7' = '#D874FD','TBA8' = '#E56DF5','HC' = '#EF67EB','TBA9' = '#F763E0','TBA10' = '#FD61D3','TBA11' = '#FF61C5','PEP1' = '#FF63B6','TBA12' = '#FF67A5','TBA13' = '#FF6B94','TBA14' = '#FF6C91'))
brain_most_stable_clusters_fly<- readRDS("brain_most_stable_clusters_fly.rds")
fly.plot <- DimPlot(brain_most_stable_clusters_fly, reduction = 'tsne', pt.size = 1, label=TRUE, cols = c('Not mapped'='#aeadb3', 'T1 neurons (0,31)' ='#EC8239', 'ensheathing glia (0,54)'= '#E48800','Dopaminergic PAM neurons (0,71)'='#7CAE00', 'Gamma???Kenyon cells (0,27)'='#00C096', 'C3 neurons (0,27)'='#00B7E9', 'Alpha/beta???Kenyon cells (0,25)'='#7C96FF'))
fly.plot
brain_most_stable_clusters_mouse <- readRDS("brain_most_stable_clusters_mouse.rds")
mouse.plot <- DimPlot(brain_most_stable_clusters_mouse, reduction = 'tsne', pt.size = 1, label=FALSE, cols = c('not mapped'='#aeadb3', 'Endothelial cells (0,27)'= '#C79800', 'Non-telencephalon astrocytes    Hypothalamus, Thalamus (0,68)'='#E48800', 'Dopaminergic interneurons (0,27)'='#7CAE00', 'Oligodendrocytes (0,80)'='#00BD5C', 'Microglia (0,40)'='#EF67EB', 'Telencephalon astrocytes    Cortex, Hippocampus, Striatum, Amygdala (0,32)'='#FF6B94'))
mouse.plot
samapfigure <- mouse.plot + fly.plot 
samapfigure
library(Seurat)
#octopus tsne cell types that are mapped
octopus <- DimPlot(brain_most_stable_clusters_annotated, reduction = 'tsne', pt.size = 1, label=FALSE, cols = c('not separated' = '#aeadb3','unstable' = '#aeadb3','IGL2- GLUT/DOP' = '#EC8239','GLIA1' = '#E48800','OGL2-DOP' = '#aeadb3','DOP2' = '#aeadb3','FBL' = '#C79800','SUB' = '#aeadb3','OGL1' = '#aeadb3','GLUT1' = '#aeadb3','DOP3' = '#aeadb3','SERT' = '#7CAE00','GLUT3' = '#aeadb3','DOP1' = '#aeadb3','IGL3' = '#aeadb3','TBA1' = '#aeadb3','PROG' = '#00BD5C','TBA2' = '#aeadb3','TBA3' = '#aeadb3','VL' = '#00C096','GABA' = '#aeadb3','IGL1-OA' = '#aeadb3','GLUT2' = '#aeadb3','OA' = '#aeadb3','TBA4' = '#aeadb3','ACH' = '#aeadb3','EC' = '#aeadb3','PEP2' = '#aeadb3','GLIA2' = '#aeadb3','GLUT4' = '#aeadb3','TBA5' = '#7C96FF','OGL3-OA' = '#aeadb3','CCAP' = '#B385FF','TBA6' = '#aeadb3','TBA7' = '#aeadb3','TBA8' = '#aeadb3','HC' = '#EF67EB','TBA9' = '#aeadb3','TBA10' = '#aeadb3','TBA11' = '#aeadb3','PEP1' = '#aeadb3','TBA12' = '#aeadb3','TBA13' = '#aeadb3','TBA14' = '#aeadb3'))
octopus        

#FLY
brain_most_stable_clusters_fly<- readRDS("brain_most_stable_clusters_fly.rds")
fly.plot <- DimPlot(brain_most_stable_clusters_fly, reduction = 'tsne', pt.size = 0.5, label=FALSE, cols = c('Not mapped'='#CCCCCC', 'T1 neurons (0,31)' ='#3C5488FF', 'ensheathing glia (0,54)'= '#00A087FF','Dopaminergic PAM neurons (0,71)'='#4DBBD5FF', 'Gamma-Kenyon cells (0,27)'='#8491B4FF', 'c3 neurons (0,27)'='#7E6148FF', 'Alpha/beta-Kenyon cells (0,25)'='#DC0000FF', 'Plasmocytes (0,69)' = '#CCCCCC'))
fly <- fly.plot + theme(legend.position = "bottom")
fly

#MOUSE
brain_most_stable_clusters_mouse <- readRDS("brain_most_stable_clusters_mouse.rds")
mouse.plot <- DimPlot(brain_most_stable_clusters_mouse, reduction = 'tsne', pt.size = 0.5, label=FALSE, cols = c('not mapped'='#CCCCCC', 'Endothelial cells (0,27)'= '#91D1C2FF', 'Non-telencephalon astrocytes    Hypothalamus, Thalamus (0,68)'='#00A087FF', 'Dopaminergic interneurons (0,27)'='#4DBBD5FF', 'Oligodendrocytes (0,80)'='#CC0099', 'Microglia (0,40)'='#F39B7FFF', 'Telencephalon astrocytes    Cortex, Hippocampus, Striatum, Amygdala (0,32)'='#FFCC33'))
mouse <- mouse.plot + theme(legend.position = "bottom")
mouse
samapfigure <- mouse + fly
samapfigure

#OCTOPUS
#read in same colours as for the original tsne
annotation <- c("not mapped","not mapped","not mapped","IGL2-GLUT/DOP","GLIA1","not mapped","not mapped","not mapped","not mapped","not mapped","FBL","not mapped","not mapped","not mapped","not mapped","IGL2-GLUT/DOP","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","SERT","not mapped","not mapped","not mapped","not mapped","not mapped","PREC","not mapped","not mapped","not mapped","not mapped","not mapped","VL","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","ACH1","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","TBA1","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","HC","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","GLIA3","not mapped","not mapped","not mapped","not mapped")
annotation <- c("not mapped","not mapped","not mapped","IGL2-GLUT/DOP","GLIA1","not mapped","not mapped","not mapped","not mapped","not mapped","FBL","not mapped","not mapped","not mapped","not mapped","IGL2-GLUT/DOP","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","SERT","not mapped","not mapped","not mapped","not mapped","not mapped","PREC","not mapped","not mapped","not mapped","not mapped","not mapped","VL","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","ACH1","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","TBA1","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","HC","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","not mapped","GLIA3","not mapped","not mapped","not mapped","not mapped")
names(annotation) <- levels(brain_most_stable_clusters_annotated)
brain_most_stable_clusters_octopus <- RenameIdents(brain_most_stable_clusters_annotated, annotation)


#change color scheme to more distinctive colors
install.packages("scales")
library(scales)
install.packages("ggsci")
library(ggsci)
show_col(pal_npg("nrc")(10))
list <- pal_npg("nrc")9)
list
octopus.plot <- DimPlot(brain_most_stable_clusters, reduction = 'tsne', pt.size = 0.5, label=FALSE, cols = c('not mapped'='#CCCCCC', 'unstable' = '#CCCCCC', 'VL' = '#DC0000FF', 'SUB' = '#3C5F88FF', 'IGL3' = '#FFCC33', 'IGL1-OA'='#F39B7FFF', 'OGL2-DOP'='#CC0099', 'OGL1'='#4DBBD5FF', 'OGL3-OA' ='#00A087FF'))

### Prepare heatmaps of glial signature in flies and mouse ###
library(Seurat)

#read in fly seurat object
flyso <- readRDS("Davieetal2018flybraindataset57kwannotation.rds")
setwd("C:/Users/SZM/OneDrive - KU Leuven/Single-cell paper 2021/first draft/seurat/samap/create heatmaps mouse and fly")
#Import metadata
GSE107451_DGRP.551_w1118_WholeBrain_57k_Metadata <- read.delim("~/temporary files R studio/samap/create heatmaps mouse and fly/GSE107451_DGRP-551_w1118_WholeBrain_57k_Metadata.tsv", row.names=1)
#Import matrix
expression_matrix <- ReadMtx(mtx = "matrix.mtx",cells = "barcodes.tsv", features = "genes.tsv", feature.column = 2)

flyso <- CreateSeuratObject(counts = expression_matrix, assay = "RNA", meta.data = GSE107451_DGRP.551_w1118_WholeBrain_57k_Metadata)
flyso <- NormalizeData(flyso)
flyso <- ScaleData(flyso)
DefaultAssay(flyso) <- "RNA"
Idents(flyso) <- "annotation"
list(flyso@active.ident)
list(flyso@meta.data$annotation)
celltypes_list <- table(Idents(flyso))
write.csv(celltypes_list,'celltypes_list_fly.csv')

#save rds 
saveRDS(flyso, file = "Davieetal2018flybraindataset57kwannotation.rds")

#Prepare plots
x <- DoHeatmap(subset(flyso, downsample = 100), features = fly_new, size = 3)
x + guides(col = "none") + theme(legend.position="top")
x 
levels(flyso)

x <- DoHeatmap(flyso, features = fly_new, size = 3)
x + guides(col = "none") + theme(legend.position="top")


fly_new <- c("pnt", "Tis11","CAP","N","whd","Irk2","Gs2","Fit1","Pdk","if","LRP1","RhoGEF2","Sap-r","Eaat1","ced-6","CG30069","CG16974","CG8745","kibra","CG7766")

#ensheating glia versus all
new.cluster.ids <- c("Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Ensheating glia","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other")
names(new.cluster.ids) <- levels(flyso)
flyso <- RenameIdents(flyso, new.cluster.ids)
DoHeatmap(subset(flyso, downsample = 100), features = fly_new, size = 3)

#read in mouse seurat object
mouseso <- readRDS("Kleshchevnikov_snRNA_brain_slice_8_All.rds.gz")
DefaultAssay(mouseso) <- "RNA"
mouse_new <- c("Etv4","Zfp36l1", "Sorbs1","Notch1", "Cpt1a", "Kcnj10", "Glul","Fermt2","Pdk2","Itga6","Lrp1","Arhgef12","Psap","Slc1a2","Gulp1","Epb41l2","Lrig1", "Etnppl","Wwc1","Phka1")
mouseso <- SetIdent(mouseso, value = mouseso@meta.data$annotation_carmen)
levels(mouseso)
celltypes_list_mouse <- table(Idents(mouseso))
write.csv(celltypes_list_mouse,'celltypes_list_mouse.csv')
DoHeatmap(subset(mouseso, downsample = 100), features = mouse_new, size = 3)

#group all glial and other cells to plot expression
#non-telencephalic astrocytes versus all other tissues
new.cluster.ids <- c("Other","Other","Other","ACNT","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other")
#all glia : ACNT, ACTE,MGL,IOL,OGC, OPC +ENDO
new.cluster.ids <- c("Other","Glia","Glia","Glia","Other","Other","Other","Other","Glia","Other","Other","Glia","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Glia","Other","Other","Other","Other","Other","Other","Glia","Other","Other","Other","Other","Other","Other","Other","Glia","Other","Other","Other","Other","Other","Other","Glia","Other","Other","Other")
names(new.cluster.ids) <- levels(mouseso)
mouseso <- RenameIdents(mouseso, new.cluster.ids)
DoHeatmap(subset(mouseso, downsample = 100), features = mouse_new, size = 3)
saveRDS(mouseso, file = "mousesoACNTvsall.rds")

#read in octopus seurat object
octopus_new <- c("LOC115211946","LOC115229692","LOC115212426","LOC115222834","LOC115230435","LOC115222850","LOC115222946","LOC115210280","LOC115220766","LOC115215786","LOC115223367","LOC115211684","LOC115209891","LOC115216632","LOC115220978","LOC115218249","LOC115217235","LOC115213940","LOC115222587","LOC115213880")
octoso <- readRDS("seurat_object_with_annotation.rds")
celltypes_list_octo <- table(Idents(octoso))
write.csv(celltypes_list_octo,'celltypes_list_octo.csv')

new.cluster.ids <- c("Other","Other","Other","Other","Other","GLIA1","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other","Other")
names(new.cluster.ids) <- levels(octoso)
octoso <- RenameIdents(octoso, new.cluster.ids)
DoHeatmap(subset(octoso, downsample = 100), features = octopus_new, size = 3)

### Creating plots for genes underlying mapping ###
#Creating fly/octopus plots

octopus_fly_colors <- c('#3C5588', '#4DBAD4','#7F6249', '#DB0C15','#838FB2','black')

#Read in fly data
flyso <- readRDS("Davieetal2018flybraindataset57kwannotation.rds")
allflygenes_r <- c("CG42492","Frq1","CG42450","RapGAP1","crol","CG9650","emc","Mnt","SoxN","cic","Vmat","JhI-21","Ddc","Pu","5-HT1A","nAChRalpha1","pdm3","Rbp6","Argk","CG44422","CG2082","sr","scrt","CG34357","CG30158","mAChR-A","CG42322","Hrb27C","Pka-R2","rut","dnc","Deaf1","Eip93F","Mef2","zfh2")

#Create the plot
pdf(file = "heatmap_fly_colors.pdf",
    width = 23, # The width of the plot in inches
    height = 10)

library(Seurat)
library(tidyverse)
for (clrs in 1:length(octopus_fly_colors)){
  
  dots <- DotPlot(flyso, features = allflygenes_r, cols = c("white", octopus_fly_colors[clrs]))  + RotatedAxis() + theme(axis.text.x = element_text(size = 12)) + theme(axis.text.y = element_text(size = 12)) + theme(legend.title = element_text(size = 10)) + theme(legend.text = element_text(size = 7)) +theme(legend.position = "top")
  dots <- dots + coord_flip()
  print(dots)
}
print(dots)
dev.off()

alloctoflygenes_r <- c("LOC115218580","LOC115215464","LOC115231120","LOC115224716","LOC115218829","LOC115217234","LOC115228298","LOC115215103","LOC115229657","LOC115221753","LOC115231184","LOC115227414","LOC115209145","LOC115209385","LOC115216428","LOC115217980","LOC115219688","LOC115212906","LOC115209225","LOC115219879","LOC115219289","LOC115228662","LOC115212072","LOC115209871","LOC115232520","LOC115209874","LOC115215688","LOC115209444","LOC115212290","LOC115219753","LOC115219746","LOC115226291","LOC115219287","LOC115211031","LOC115220811")

#Create the plot
pdf(file = "heatmap_octopus_fly_colors.pdf",
    width = 23, # The width of the plot in inches
    height = 10)

library(Seurat)
library(tidyverse)
for (clrs in 1:length(octopus_fly_colors)){
  
  dots <- DotPlot(stableclusters, features = alloctoflygenes_r, cols = c("white", octopus_fly_colors[clrs]))  + RotatedAxis() + theme(axis.text.x = element_text(size = 12)) + theme(axis.text.y = element_text(size = 12)) + theme(legend.title = element_text(size = 10)) + theme(legend.text = element_text(size = 7)) +theme(legend.position = "top")
  dots <- dots + coord_flip()
  print(dots)
}
print(dots)
dev.off()

#Creating mouse/octopus plots

octopus_mouse_colors <- c('black','#B83289','#4DBAD4','#019F86','#F39A7F','#FDCA38','#92CDC0' )

#Read in mouse data
setwd("C:/Users/SZM/OneDrive - KU Leuven/Single-cell paper 2021/first draft/seurat/samap/create heatmaps mouse and fly")
mouseso <- readRDS("Kleshchevnikov_snRNA_brain_slice_8_All.rds.gz")
DefaultAssay(mouseso) <- "RNA"
mouseso <- SetIdent(mouseso, value = mouseso@meta.data$annotation_carmen)

allmousegenes_r <- c("Dhrs3","Pdk2","Kcnj16","Cdh20","Epb41l2","Rorb","Rora","Pax5","Slc6a3","Chrna6","Ddc","Gch1","Fryl","Pbx3","Satb1","Adgre1","Col15a1","Stab1","Msn","Mef2a","Mef2c","Ptpn13","Col1a2","Col4a5","Igfbp7","Svil","Ebf1","Foxp1","Actb","Piezo2","Phldb1","Neo1","Arnt2","Baz2b","Tet1","Mta3","Rere","Smad7","Tcf12","Zeb2")
#Create the plot
pdf(file = "octopus_heatmap_mouse_colors.pdf",
    width = 20, # The width of the plot in inches
    height = 10)


for (clrs in 1:length(octopus_mouse_colors)){

dots <- DotPlot(mouseso, features = allmousegenes_r, cols = c("white", octopus_mouse_colors[clrs]))  + RotatedAxis() + theme(axis.text.x = element_text(size = 12)) + theme(axis.text.y = element_text(size = 12)) + theme(legend.title = element_text(size = 10)) + theme(legend.text = element_text(size = 7)) +theme(legend.position = "top")
dots <- dots + coord_flip()
print(dots)
}
print(dots)
dev.off()

#Read in octopus data
setwd("~/temporary files R studio/preparing files for GEO")
stableclusters <- readRDS("stable_clusters.rds")

alloctomousegenes_r <- c("LOC115210961","LOC115220766","LOC115222850","LOC115216961","LOC115218249","LOC115209615","LOC115232352","LOC115216435","LOC115213095","LOC115217980","LOC115209145","LOC115209385","LOC115221035","LOC115210669","LOC115215778","LOC115214950","LOC115223963","LOC115219218","LOC115225885","LOC115211138","LOC115211031","LOC115232620","LOC115213426","LOC115221996","LOC115219103","LOC115211517","LOC115223202","LOC115210326","LOC115214120","LOC115221801","LOC115211639","LOC115216496","LOC115220540","LOC115220806","LOC115214204","LOC115210243","LOC115224734","LOC115218995","LOC115219652","LOC115221075")


#Create the plot
pdf(file = "octopus_heatmap_m_colors.pdf",
    width = 20, # The width of the plot in inches
    height = 10)


for (clrs in 1:length(octopus_mouse_colors)){
  
  dots <- DotPlot(stableclusters, features = alloctomousegenes_r, cols = c("white", octopus_mouse_colors[clrs]))  + RotatedAxis() + theme(axis.text.x = element_text(size = 12)) + theme(axis.text.y = element_text(size = 12)) + theme(legend.title = element_text(size = 10)) + theme(legend.text = element_text(size = 7)) +theme(legend.position = "top")
  dots <- dots + coord_flip()
  print(dots)
}
print(dots)
dev.off()


