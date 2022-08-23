#Install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)
install.packages("gplots")
library(gplots)

data = read.table("counts.txt", header=TRUE, row.names=1)
data <- data[,-1]
dim(as.matrix(new_df))

#samples
#IZ10.bam SEM
#IZ11.bam SUB
#IZ12.bam OL
#IZ13.bam ARM
#IZ15.bam SEM
#IZ16.bam SUB
#IZ17.bam OL
#IZ18.bam ARM
#IZ20.bam SEM
#IZ21.bam SUB
#IZ22.bam OL
#IZ23.bam ARM

#remove samples that are not of interest
new_df = subset(data, select = -c(13:25)) 
new_df <- subset(new_df, select = -c(IZ13.bam, IZ18.bam,IZ23.bam))
new_df <- subset(new_df, select = -c(GDC1.bam, GDC2.bam,GDC3.bam,GDC7.bam,GDC8.bam,GDC9.bam,GDC10.bam,GDC11.bam,GDC12.bam))

#prefilter data (detected in 2 samples with more than one read)
filter <- apply(new_df, 1, function(x) length(x[x>1])>=2)
df <- new_df[filter,]

samples <- read.csv("~/conditions.txt", sep="")

dds <- DESeqDataSetFromMatrix(countData=df, colData=samples, design= ~ lobe) 
dds <- DESeqDataSetFromMatrix(countData=df, colData=samples, design= ~ condition) 

dds$lobe <- relevel(dds$lobe,"MU")
dds$condition <- relevel(dds$condition,"central")
dds <- DESeq(dds)

sizeFactors(dds) #which number is used for normalisation
colData(dds) #each sample with its size factor per condition
normalized_counts <- counts(dds, normalized=TRUE) #you look for the counts that are normalized

write.table(normalized_counts, file="normalized_counts.txt")
write.csv(normalized_counts, "normalized_counts.csv")
write.csv(df, "raw_counts.csv")

###check for variation sources/explore data with PCA and heatmap
rld <- rlog(dds, blind=TRUE) #data transformation for visualization: log transformation
plotPCA(rld, intgroup="condition") #results are moved to the vector rld
pcaData <- plotPCA(rld, intgroup=c("condition", "lobe"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=lobe)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
library(ggplot2)
rld_mat <- assay(rld)

##making a heatmap to visualise within sample variation: lighter colors: higher correlation. Red : lower correlation 
rld_cor <- cor(rld_mat) #computes pairwise correlation values between the samples: vulgaris-vulgaris, vulgaris-bimaculoides,...
heatmap.2(rld_cor, key=TRUE, trace="none")

res <- results(dds)
res
write.table(res, file="res.txt")

#extract upregulated genes
Res$Expression
data_frame_mod <- filter(Res,Expression=="Up-regulated")
DE_optic_vs_central <- rownames(data_frame_mod)
write.table(DE_optic_vs_central, file="Up_optic_vs_central.txt")

#filter for genes detected in single cell
optic_highly_expressed <- read.table("~/optic_highly_expressed.txt", quote="\"", comment.char="")
optic <- optic_highly_expressed$V1

#extract downregulated genes
Res$Expression
data_frame_mod <- filter(Res,Expression=="Down-regulated")
Down_optic_vs_central <- rownames(data_frame_mod)
write.table(Down_optic_vs_central, file="Down_optic_vs_central.txt")
#filter for genes detected in single cell
central_highly_expressed <- read.table("~/central_highly_expressed.txt", quote="\"", comment.char="")
central <- central_highly_expressed$V1
sessionInfo("DESeq2")

#create module scoring for genes upregulated in optic versus central
setwd("~/temporary files R studio/preparing files for GEO")
brain <- readRDS("seurat_object_with_annotation.rds")

##OPTIC
#100 genes,smallest p adj + most neg fold change
OL <- c("LOC115222757","LOC115215914","LOC115212921","LOC115215570","LOC115222878","LOC115213528","LOC118762162","LOC115218923","LOC115224902","LOC115210951","LOC115216514","LOC115216420","LOC118764039","LOC115209669","LOC118762160","LOC115216765","LOC115224021","LOC115220952","LOC115232533","LOC115223946","LOC115211073","LOC115219879","LOC118766480","LOC115209497","LOC115230791","LOC115223700","LOC115218977","LOC115223080","LOC115219720","LOC115216876","LOC115222357","LOC118761763","LOC115209184","LOC115211044","LOC115211675","LOC115209653","LOC115232563","LOC118767514","LOC118765869","LOC118765039","LOC115219369","LOC115224741","LOC118762414","LOC115231020","LOC115218238","LOC115220538","LOC115221267","LOC115213418","LOC115215738","LOC115217871","LOC115209714","LOC115219813","LOC115232578","LOC115232597","LOC115211998","LOC115232618","LOC115215208","LOC118765160","LOC115211969","LOC115216584","LOC115211100","LOC115215596","LOC118760920","LOC115215679","LOC118761015","LOC115216729","LOC115224495","LOC115218927","LOC115209434","LOC118763813","LOC115213213","LOC115219368","LOC115216690","LOC115209420","LOC115218974","LOC115227516","LOC115209806","LOC115213422","LOC115214993","LOC115230030","LOC115218955","LOC115209262","LOC115218997","LOC115209122","LOC115209324","LOC115216384","LOC115217069","LOC118764325","LOC118764777","LOC115232562","LOC115218736","LOC115225708","LOC118767594","LOC115223212","LOC115227587","LOC115218045","LOC115221800","LOC115218983","LOC115214842")
library(Seurat)
brain <- AddModuleScore(object = brain, features = optic, name ="Optic")
plot <- FeaturePlot(object = brain, features ="Optic1", reduction = "tsne",min.cutoff = 0, max.cutoff =2.5)
plot + ggtitle("Optic lobe")

##CENTRAL
#100 genes, smallest p adj + most pos fold change
CB <- c("LOC115227037","LOC115227860","LOC115214113","LOC115220836","LOC115215075","LOC115213033","LOC115221654","LOC115213194","LOC115214132","LOC115222993","LOC115218782","LOC115212087","LOC115213084","LOC115209753","LOC115211758","LOC115213836","LOC115209311","LOC115209331","LOC115215588","LOC115212402","LOC115214226","LOC118763170","LOC115222738","LOC118763027","LOC115232232","LOC115220974","LOC118763091","LOC115210485","LOC115226112","LOC115215674","LOC118762907","LOC115222573","LOC115214172","LOC115213157","LOC115213333","LOC118765708","LOC115230630","LOC118763633","LOC115218422","LOC115222298","LOC115227414","LOC115224826","LOC115211890","LOC115217300","LOC118761157","LOC115210249","LOC115211767","LOC115210843","LOC115220808","LOC115215372","LOC115223411","LOC115220878","LOC118765481","LOC115210839","LOC115223018","LOC115218403","LOC115228951","LOC115232511","LOC115222154","LOC118762557","LOC115219244","LOC115213334","LOC115210675","LOC115223677","LOC115224753","LOC115210562","LOC115218581","LOC115213848","LOC115218560","LOC115212337","LOC115211991","LOC118763555","LOC115215760","LOC115218755","LOC115217137","LOC115223559","LOC118764548","LOC115217493","LOC115222185","LOC115209145","LOC115229190","LOC115220942","LOC115218137","LOC115221726","LOC115212950","LOC115231038","LOC115225874","LOC115210576","LOC115227650","LOC118763461","LOC115218280","LOC118763099","LOC115214214","LOC115211627","LOC115219536","LOC115222011","LOC115226756","LOC115222540","LOC115221105","LOC115211123"
      )#50 genes 
CB50 <- c("LOC115227037","LOC115227860","LOC115214113","LOC115220836","LOC115215075","LOC115213033","LOC115221654","LOC115213194","LOC115214132","LOC115222993","LOC115218782","LOC115212087","LOC115213084","LOC115209753","LOC115211758","LOC115213836","LOC115209311","LOC115209331","LOC115215588","LOC115212402","LOC115214226","LOC118763170","LOC115222738","LOC118763027","LOC115232232","LOC115220974","LOC118763091","LOC115210485","LOC115226112","LOC115215674","LOC118762907","LOC115222573","LOC115214172","LOC115213157","LOC115213333","LOC118765708","LOC115230630","LOC118763633","LOC115218422","LOC115222298","LOC115227414","LOC115224826","LOC115211890","LOC115217300","LOC118761157","LOC115210249","LOC115211767","LOC115210843","LOC115220808")
brain <- AddModuleScore(object = brain, features = CB50, name ="Central")
plot <- FeaturePlot(object = brain, features ="Central1", reduction = "tsne",min.cutoff = 0, max.cutoff = 2.5)
plot + ggtitle("Central brain")


#saving plots as landscape A4 8,... x 10 inches

