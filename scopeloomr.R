#module load R/4.0.0-foss-2018a-X11-20180604
#module load GCC/6.4.0-2.28
#module load OpenBLAS/0.2.20-GCC-6.4.0-2.28
#module load binutils/2.32-GCCcore-6.4.0
#module load HDF5/1.10.5-foss-2018a
#module load Pandoc/2.5
#module load OpenSSL/1.1.1g-GCCcore-6.4.0
#module load cURL/7.65.3-GCCcore-6.4.0
options(width = 200)
library(data.table)
#############################################################################################
library(SCopeLoomR)	#v0.13.0
library(Seurat)

wd <- file.path(staging, "/OEB/OEBlarvaeClustRuth")
setwd(wd)

brain <- readRDS("brain_most_stable_clusters_w_annotation.rds")
brain <- brain_most_stable_clusters_annotated

as.data.table(brain@active.ident)[, .N, by = V1][order(N, decreasing = T)]
#               V1    N
# 1:      unstable 7119
# 2: not separated 2236
# 3: IGL2-GLUT/DOP  917
# 4:         GLIA1  537
# 5:      OGL2-DOP  520
# 6:          DOP2  450
# 7:           FBL  430
# 8:           SUB  420
#...

markers <- readRDS("brain_seurat_markers.rds.gz")
as.data.table(markers)[, .N, by = cluster][order(N, decreasing = T)]
#          cluster    N
# 1:          PREC 1954
# 2:         GLIA2 1853
# 3:         GLIA1 1806
# 4:         PRECL 1716
# 5:            EC 1617
# 6:          IGL3 1400
#...

build_loom(file.name = "paralarval_brain_annotated.loom",
           dgem = brain@assays$RNA@counts,
           title = "paralarval_brain_annotated",
           default.embedding = brain@reductions$tsne@cell.embeddings,
           default.embedding.name = "tsne")

loom <- open_loom("paralarval_brain_annotated.loom", mode = "r+")


add_embedding(loom = loom, 
              embedding = brain@reductions$umap@cell.embeddings,
              name = "umap")
add_embedding(loom = loom, 
              embedding = brain@reductions$pca@cell.embeddings,
              name = "pca")
add_col_attr(loom = loom, key = "percent.mito",
             value = brain@meta.data$percent.mt, as.metric = TRUE)
add_col_attr(loom = loom, key = "nCount",
             value = brain@meta.data$nCount_RNA, as.metric = TRUE)
add_col_attr(loom = loom, key = "nFeature",
             value = brain@meta.data$nFeature_RNA, as.metric = TRUE)
add_col_attr(loom = loom, key = "batch",
             value = brain@meta.data$orig.ident, as.annotation = TRUE)


annotation <- data.table(celltype = levels(brain@active.ident),
	ct_id = 0:(length(levels(brain@active.ident))-1))
annot_list <- setNames(annotation[, ct_id], annotation[, celltype])

ct_all <- setNames(as.vector(brain@active.ident), rownames(brain@meta.data))
cl_all <- setNames(annot_list[as.vector(ct_all)], rownames(brain@meta.data))
brain@meta.data$SCT_snn_res.2 <- as.vector(cl_all)

add_annotated_clustering(loom = loom,
	group = "Seurat",
	name = "SCT_snn_res.2",
	clusters = cl_all,
	annotation = ct_all,
	is.default = T,
	overwrite.default = T)

#rename markers to cluster ids
markers$cluster <- as.vector(annot_list[markers$cluster])
add_clustering_markers(loom = loom,
	clustering.id = 0,
	clustering.markers = split(markers, markers$cluster),
	marker.metric.accessors = c("avg_log2FC", "p_val_adj"),
	marker.metric.names = c("Avg. logFC", "adjusted P-value"),
	marker.metric.description = c("Average log fold change", "Adjusted p-value (BF)"))


add_hierarchy(loom = loom, hierarchy = create_hierarchy(level.1.name = "Octopus_vulgaris", level.2.name = "brain"))
close_loom(loom)


###################################################Create second loom file with seurat clusterings
brain@active.ident <- brain@meta.data$"seurat_clusters"


brain <- readRDS("brain_most_stable_clusters_w_annotation.rds")
brain <- brain_most_stable_clusters_annotated

as.data.table(brain@active.ident)[, .N, by = V1][order(N, decreasing = T)]
#               V1    N
# 1:      unstable 7119
# 2: not separated 2236
# 3: IGL2-GLUT/DOP  917
# 4:         GLIA1  537
# 5:      OGL2-DOP  520
# 6:          DOP2  450
# 7:           FBL  430
# 8:           SUB  420
#...

markers <- readRDS("brain_seurat_markers_or.rds.gz")
as.data.table(markers)[, .N, by = cluster][order(N, decreasing = T)]
#          cluster    N
# 1:          PREC 1954
# 2:         GLIA2 1853
# 3:         GLIA1 1806
# 4:         PRECL 1716
# 5:            EC 1617
# 6:          IGL3 1400
#...

build_loom(file.name = "paralarval_brain_seurat.loom",
           dgem = brain@assays$RNA@counts,
           title = "paralarval_brain_seurat",
           default.embedding = brain@reductions$tsne@cell.embeddings,
           default.embedding.name = "tsne")

loom <- open_loom("paralarval_brain_seurat.loom", mode = "r+")


add_embedding(loom = loom, 
              embedding = brain@reductions$umap@cell.embeddings,
              name = "umap")
add_embedding(loom = loom, 
              embedding = brain@reductions$pca@cell.embeddings,
              name = "pca")
add_col_attr(loom = loom, key = "percent.mito",
             value = brain@meta.data$percent.mt, as.metric = TRUE)
add_col_attr(loom = loom, key = "nCount",
             value = brain@meta.data$nCount_RNA, as.metric = TRUE)
add_col_attr(loom = loom, key = "nFeature",
             value = brain@meta.data$nFeature_RNA, as.metric = TRUE)
add_col_attr(loom = loom, key = "batch",
             value = brain@meta.data$orig.ident, as.annotation = TRUE)


annotation <- data.table(celltype = levels(brain@active.ident),
                         ct_id = 0:(length(levels(brain@active.ident))-1))
annot_list <- setNames(annotation[, ct_id], annotation[, celltype])

ct_all <- setNames(as.vector(brain@active.ident), rownames(brain@meta.data))
cl_all <- setNames(annot_list[as.vector(ct_all)], rownames(brain@meta.data))
brain@meta.data$SCT_snn_res.2 <- as.vector(cl_all)

add_annotated_clustering(loom = loom,
                         group = "Seurat",
                         name = "SCT_snn_res.2",
                         clusters = cl_all,
                         annotation = ct_all,
                         is.default = T,
                         overwrite.default = T)

#rename markers to cluster ids
markers$cluster <- as.vector(annot_list[markers$cluster])
add_clustering_markers(loom = loom,
                       clustering.id = 0,
                       clustering.markers = split(markers, markers$cluster),
                       marker.metric.accessors = c("avg_log2FC", "p_val_adj"),
                       marker.metric.names = c("Avg. logFC", "adjusted P-value"),
                       marker.metric.description = c("Average log fold change", "Adjusted p-value (BF)"))


add_hierarchy(loom = loom, hierarchy = create_hierarchy(level.1.name = "Octopus_vulgaris", level.2.name = "brain"))
close_loom(loom)

sessionInfo() #{
#R version 4.0.0 (2020-04-24)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: CentOS Linux 7 (Core)
#
#Matrix products: default
#BLAS/LAPACK: /vsc-hard-mounts/leuven-apps/skylake/2018a/software/OpenBLAS/0.2.20-GCC-6.4.0-2.28/lib/libopenblas_haswellp-r0.2.20.so
#
#locale:
# [1] LC_CTYPE=en_US.utf8       LC_NUMERIC=C              LC_TIME=en_US.utf8       
# [4] LC_COLLATE=C              LC_MONETARY=en_US.utf8    LC_MESSAGES=en_US.utf8   
# [7] LC_PAPER=en_US.utf8       LC_NAME=C                 LC_ADDRESS=C             
#[10] LC_TELEPHONE=C            LC_MEASUREMENT=en_US.utf8 LC_IDENTIFICATION=C      
#
#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     
#
#other attached packages:
#[1] SeuratObject_4.0.2 Seurat_4.0.4       SCopeLoomR_0.13.0  data.table_1.14.2 
#
#loaded via a namespace (and not attached):
#  [1] nlme_3.1-147          matrixStats_0.61.0    spatstat.sparse_2.0-0 bit64_4.0.5          
#  [5] RcppAnnoy_0.0.19      RColorBrewer_1.1-2    httr_1.4.2            sctransform_0.3.2    
#  [9] tools_4.0.0           utf8_1.2.2            R6_2.5.1              irlba_2.3.3          
# [13] rpart_4.1-15          KernSmooth_2.23-17    uwot_0.1.10.9000      mgcv_1.8-31          
# [17] DBI_1.1.1             lazyeval_0.2.2        colorspace_2.0-2      tidyselect_1.1.1     
# [21] gridExtra_2.3         bit_4.0.4             compiler_4.0.0        hdf5r_1.3.4          
# [25] plotly_4.9.4.1        scales_1.1.1          lmtest_0.9-38         spatstat.data_2.1-0  
# [29] ggridges_0.5.3        pbapply_1.5-0         goftest_1.2-2         stringr_1.4.0        
# [33] digest_0.6.28         spatstat.utils_2.2-0  pkgconfig_2.0.3       htmltools_0.5.2      
# [37] parallelly_1.28.1     fastmap_1.1.0         htmlwidgets_1.5.4     rlang_0.4.11         
# [41] shiny_1.7.1           generics_0.1.0        zoo_1.8-9             jsonlite_1.7.2       
# [45] ica_1.0-2             dplyr_1.0.7           magrittr_2.0.1        patchwork_1.1.1      
# [49] Matrix_1.3-4          Rcpp_1.0.7            munsell_0.5.0         fansi_0.5.0          
# [53] abind_1.4-5           reticulate_1.22       lifecycle_1.0.1       stringi_1.7.5        
# [57] MASS_7.3-51.6         Rtsne_0.15            plyr_1.8.6            grid_4.0.0           
# [61] parallel_4.0.0        listenv_0.8.0         promises_1.2.0.1      ggrepel_0.9.1        
# [65] crayon_1.4.1          deldir_1.0-2          miniUI_0.1.1.1        lattice_0.20-41      
# [69] cowplot_1.1.1         splines_4.0.0         tensor_1.5            pillar_1.6.3         
# [73] igraph_1.2.8          rjson_0.2.20          spatstat.geom_2.2-2   future.apply_1.8.1   
# [77] reshape2_1.4.4        codetools_0.2-16      leiden_0.3.9          glue_1.4.2           
# [81] png_0.1-7             vctrs_0.3.8           httpuv_1.6.3          gtable_0.3.0         
# [85] RANN_2.6.1            purrr_0.3.4           spatstat.core_2.3-0   polyclip_1.10-0      
# [89] tidyr_1.1.4           scattermore_0.7       future_1.22.1         assertthat_0.2.1     
# [93] ggplot2_3.3.5         mime_0.12             xtable_1.8-4          later_1.3.0          
# [97] survival_3.1-12       viridisLite_0.4.0     tibble_3.1.5          cluster_2.1.0        
#[101] globals_0.14.0        fitdistrplus_1.1-6    ellipsis_0.3.2        ROCR_1.0-11        
#}
