library(Seurat)
library(dplyr)
library(limma)
setwd("/Applications/UPenn/SingleCellSequencing")
clone_tags_pos = read.table("In_order_Clone_tags", header = TRUE,sep="\t",stringsAsFactors = FALSE)
df.genes1 = read.table("chromosome_name", header = TRUE,sep="\t",stringsAsFactors = FALSE)
pbmc.data <- Read10X(data.dir="/Applications/UPenn/SingleCellSequencing/A549_RFPPos_2891732785_GRCh38_singlecell/outs/filtered_gene_bc_matrices/GRCh38")
pbmc <- CreateSeuratObject(pbmc.data)
pbmc
pbmc[["clone_tags"]]=clone_tags_pos$tag

clone_tags_pos$tag1 = ifelse(clone_tags_pos$tag == "P1", "P1", "other")
pbmc[["clone_tags1"]]=clone_tags_pos$tag1

Idents(pbmc)=pbmc$clone_tags1


pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
pbmc <- NormalizeData(pbmc)

pbmc=FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc=ScaleData(pbmc)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:12)
pbmc <- FindClusters(pbmc, resolution = 0.1)
#Idents(pbmc)=pbmc$clone_tags
pbmc <- RunUMAP(pbmc, dims = 1:12)
Idents(pbmc)=pbmc$clone_tags1

DimPlot(pbmc, reduction = "umap",label=TRUE)
#P1_dub=FindMarkers(pbmc,ident.1 = "P1",group.by = "clone_tags",min.pct = 0.25,logfc.threshold = 0.25)
P1=FindMarkers(pbmc,ident.1 = "P1",min.pct = 0.25,logfc.threshold = 0.25)
other= FindMarkers(pbmc,ident.1 = "other",min.pct = 0.25,logfc.threshold = 0.25)

P1$gene=rownames(P1)
other$gene=rownames(other)

P1_chrom=merge(P1,df.genes1,by.x = "gene",by.y="common_name",all.x=TRUE)
other=merge(other,df.genes1,by.x = "gene",by.y="common_name",all.x=TRUE)

p1_neg=P1_chrom[P1_chrom$avg_log2FC<0 & P1_chrom$p_val_adj<0.05,]
p1_pos=P1_chrom[P1_chrom$avg_log2FC>0 & P1_chrom$p_val_adj<0.05,]
other=other[other$avg_log2FC>0 & other$p_val_adj<0.05,]
p1_pos1=p1_pos[p1_pos$chromosome_name!="MT"&p1_pos$chromosome_name!="X"&p1_pos$chromosome_name!="Y",]
other1=other[other$chromosome_name!="MT"&other$chromosome_name!="X"&other$chromosome_name!="Y",]



barplot(table(strtoi(p1_pos$chromosome_name)))

barplot(table(strtoi(other1$chromosome_name)))


a=table(strtoi(p1_pos$chromosome_name))
b=table(strtoi(other1$chromosome_name))

df=cbind(a,b)

# Create a vector of alternating colors
colors <- rep(c("lightblue", "white"), length.out = nrow(df))



# Create a bar plot with alternating colors
barplot(t(as.matrix(df)), beside = TRUE, col = colors, xlab = "Chromosome number", ylab = "Distinct # genes")
#(rect(xleft = -0.25, xright = 22, ybottom = 0, ytop = 19.3, col = "gray", border = NA,density = 0))
legend("topright", c("P1", "other"), fill = c("lightblue", "white"))
abline(h = 19.3, col = "grey")
# Add a grey box to the bar plot
