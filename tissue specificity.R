#Cell type specificity

#Average gene expression per cell type
setwd("~/first draft/seurat/integration/most_stable")
brain_most_stable_clusters<- readRDS("brain_most_stable_clusters.rds")
DefaultAssay(brain_most_stable_clusters) <- "SCT"
cluster.averages <- AverageExpression(brain_most_stable_clusters)
cluster.averages.SCT <- cluster.averages$SCT
setwd("~/first draft/tau")
write.csv(cluster.averages.SCT, "clusteraveragessct.csv")
df_avEx_celltype <-  as.data.frame(cluster.averages.SCT)
df_avEx_celltype$gene_id <- row.names(cluster.averages.SCT)
df_avEx_celltype <- df_avEx_celltype[,c(ncol(df_avEx_celltype), 1:ncol(df_avEx_celltype)-1)]
# remove low quality entry
df_avEx_celltype$LowQual <- NULL
df_avEx_celltype
write.csv(df_avEx_celltype, "df_avEx_celltype.csv")
dgem <- brain_most_stable_clusters$SCT@counts
#get cells per gene
ncells <- apply(dgem, 1, function(x){ length(x[x > 0]) } )
df_ncells <- data.frame(gene_id = names(ncells), ncells = ncells)
write.csv(df_ncells, "df_ncellse.csv")

path_out_mtx_celltype <- "~/matrix_avExp_cellType_snscRNA_brain.tsv"
write.table(df_avEx_celltype, path_out_mtx_celltype, row.names=FALSE, col.names=TRUE, quote=F, sep="\t")

#get cells per gene
dgem <- brain_most_stable_clusters$SCT@counts
ncells <- apply(dgem, 1, function(x){ length(x[x > 0]) } )

df_ncells <- data.frame(gene_id = names(ncells), ncells = ncells)

write.table(df_ncells, path_out_ncells, row.names=FALSE, col.names=TRUE, quote=F, sep="\t")
