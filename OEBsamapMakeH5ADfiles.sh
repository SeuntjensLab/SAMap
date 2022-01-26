#<[[OEB]]
#<[[OEBsamap]]
#<[[sceasy]]

# task: convert seurat files to h5ad for SAMap with [[sceasy]]
conda_initialize /staging/leuven/stg_00002/lcb/kspan/software/anaconda3
conda activate sceasy

#module load R/4.0.0-foss-2018a-X11-20180604
#module load GCC/6.4.0-2.28
#module load OpenBLAS/0.2.20-GCC-6.4.0-2.28
#module load binutils/2.32-GCCcore-6.4.0
#module load HDF5/1.10.5-foss-2018a
#module load Pandoc/2.5
.libPaths("/staging/leuven/stg_00002/lcb/kspan/R/x86_64-pc-linux-gnu-library/4.0")
staging <- "/staging/leuven/stg_00002/lcb"
library(data.table)
#############################################################################################
library(sceasy)
library(reticulate)
use_condaenv('sceasy')
loompy <- reticulate::import('loompy')
library(Seurat)

wd <- file.path(staging, "kspan/OEB/OEBlarvaeClustRuth")
setwd(wd)

## octopus data
## 17961 genes x 17081 cells
#############################
so <- readRDS("brain_most_stable_clusters.rds")
#'Idents' are ruths prefered clustering -> not in meta.data -> add
so <- AddMetaData(so, Idents(so), col.name = "stable_clusters")
sceasy::convertFormat(so, from="seurat", to="anndata", outFile='brain_most_stable_clusters.h5ad')

## fly data
###########
#[[OEBsamapClustRuthFlyI]]
fso <- readRDS("v2_fca_biohub_head_10x.rds.gz")
sceasy::convertFormat(fso, from="seurat", to="anndata", outFile= "v2_fca_biohub_head_10x.h5ad")

#does not work:
loomfile <- file.path(staging, "lcb_projects/fca/atlas/v2/loom/v2_fca_biohub_head_10x.loom")
sceasy::convertFormat(loomfile, from="loom", to="anndata",
                       outFile='v2_fca_biohub_head_10x.loom.h5ad')
#Error in py_call_impl(callable, dots$args, dots$keywords) : 
#  ValueError: Row attribute 'ClusterMarkers_-1' dtype [('0', '<i8'), ('...

#[[OEBsamapClustRuthFlyII]]
fso <- readRDS("Davie_Janssens_Koldere_et_al_2018_AdultBrain.rds.gz")
sceasy::convertFormat(fso, from="seurat", to="anndata", outFile= "Davie_Janssens_Koldere_et_al_2018_AdultBrain.h5ad")

#[[OEBsamapClustRuthFlyIII]]
fso <- readRDS("Brain_VNC_Ravenscroft_et_al_2019.rds.gz")
sceasy::convertFormat(fso, from="seurat", to="anndata", outFile= "Brain_VNC_Ravenscroft_et_al_2019.h5ad")

sessionInfo() #{
R version 4.0.0 (2020-04-24)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /vsc-hard-mounts/leuven-apps/skylake/2018a/software/OpenBLAS/0.2.20-GCC-6.4.0-2.28/lib/libopenblas_haswellp-r0.2.20.so

locale:
 [1] LC_CTYPE=en_US.utf8       LC_NUMERIC=C
 [3] LC_TIME=en_US.utf8        LC_COLLATE=C
 [5] LC_MONETARY=en_US.utf8    LC_MESSAGES=en_US.utf8
 [7] LC_PAPER=en_US.utf8       LC_NAME=C
 [9] LC_ADDRESS=C              LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.utf8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] SeuratObject_4.0.1 Seurat_4.0.1       data.table_1.12.8  sceasy_0.0.6
[5] reticulate_1.20

loaded via a namespace (and not attached):
  [1] nlme_3.1-147          spatstat.sparse_2.0-0 bitops_1.0-6
  [4] matrixStats_0.58.0    bit64_0.9-7           RcppAnnoy_0.0.18
  [7] RColorBrewer_1.1-2    httr_1.4.1            sctransform_0.3.2
 [10] tools_4.0.0           R6_2.4.1              irlba_2.3.3
 [13] rpart_4.1-15          KernSmooth_2.23-17    uwot_0.1.10
 [16] mgcv_1.8-31           lazyeval_0.2.2        colorspace_1.4-1
 [19] npsurv_0.4-0          tidyselect_1.0.0      gridExtra_2.3
 [22] bit_1.1-15.2          compiler_4.0.0        cli_2.0.2
 [25] plotly_4.9.2.1        caTools_1.18.0        scales_1.1.0
 [28] spatstat.data_2.1-0   lmtest_0.9-37         ggridges_0.5.2
 [31] pbapply_1.4-2         rappdirs_0.3.1        goftest_1.2-2
 [34] stringr_1.4.0         digest_0.6.25         spatstat.utils_2.1-0
 [37] pkgconfig_2.0.3       htmltools_0.5.1.1     fastmap_1.0.1
 [40] htmlwidgets_1.5.1     rlang_0.4.11.9000     shiny_1.4.0.2
 [43] zoo_1.8-7             jsonlite_1.6.1        ica_1.0-2
 [46] gtools_3.8.2          dplyr_0.8.5           magrittr_1.5
 [49] patchwork_1.1.1       Matrix_1.3-3          Rcpp_1.0.6
 [52] munsell_0.5.0         fansi_0.4.1           abind_1.4-5
 [55] lifecycle_0.2.0       stringi_1.4.6         MASS_7.3-51.6
 [58] gplots_3.0.3          Rtsne_0.15            plyr_1.8.6
 [61] grid_4.0.0            parallel_4.0.0        gdata_2.18.0
 [64] listenv_0.8.0         promises_1.1.0        ggrepel_0.8.2
 [67] crayon_1.3.4          deldir_0.1-25         miniUI_0.1.1.1
 [70] lattice_0.20-41       cowplot_1.0.0         splines_4.0.0
 [73] tensor_1.5            pillar_1.4.3          igraph_1.2.5
 [76] spatstat.geom_2.1-0   future.apply_1.7.0    reshape2_1.4.4
 [79] codetools_0.2-16      leiden_0.3.7          glue_1.4.0
 [82] lsei_1.2-0            png_0.1-7             vctrs_0.2.4
 [85] httpuv_1.5.2          polyclip_1.10-0       gtable_0.3.0
 [88] RANN_2.6.1            purrr_0.3.4           spatstat.core_2.1-2
 [91] tidyr_1.0.2           scattermore_0.7       future_1.17.0
 [94] assertthat_0.2.1      ggplot2_3.3.0         mime_0.9
 [97] xtable_1.8-4          later_1.0.0           survival_3.1-12
[100] viridisLite_0.3.0     tibble_3.0.1          cluster_2.1.0
[103] globals_0.12.5        fitdistrplus_1.0-14   ellipsis_0.3.0
[106] ROCR_1.0-7
#}
