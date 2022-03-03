library(Seurat)
library(filesstrings)
library(gplots)
library(ggplot2)
library(DropletUtils)
library(patchwork)
library(dplyr)
##Rscript scRNA_DE_heatmap_plot.R raw_expression.csv cluster.csv name <min_pct> <logfc> <column used in cluster file>

args=commandArgs(TRUE)
raw_expr_matrix<-args[1]
cluster<-args[2]
name<-args[3]
min_pct<-args[4]
logfc<-args[5]
cluster_column<-args[6]


if(min_pct %in% "NA"){
	min_pct<-0.1
}

if(logfc %in% "NA"){
	logfc<-0.1
}

if(cluster_column %in% "NA"){
	cluster_column<-1
}

min_pct<-as.numeric(min_pct)
logfc<-as.numeric(logfc)
cluster_column<-as.numeric(cluster_column)
print(raw_expr_matrix)
print(cluster)
print(name)
print(min_pct)
print(logfc)
print(cluster_column)



folder_name=paste(name,"_DE_heatmap_plots",sep="")
dir.create(folder_name)

data<-read.csv(raw_expr_matrix,row.name=1,check.names=F,header=T)
cluster<-read.csv(cluster,row.names=1,check.names=F,header=T)
head(cluster)
colnames(data)<-sub("\\.","-",colnames(data))
data<-data[,rownames(cluster)]
dim(data)
print(data[1:3,1:3])
print(colnames(cluster)[cluster_column])
object1 <- CreateSeuratObject(counts = data, project = "test", min.cells = 0, min.features = 0)
object1@meta.data<-cbind(object1@meta.data[intersect(rownames(cluster),rownames(object1@meta.data)),],cluster[intersect(rownames(cluster),rownames(object1@meta.data)),])
print(head(object1@meta.data))
Idents(object1)<-colnames(cluster)[cluster_column]##select which column in cluster file is used for clustering
Idents(object1)<-factor(x=Idents(object1),levels=names(sort(table(Idents(object1)),decreasing = TRUE)))
object1 <- NormalizeData(object1, normalization.method = "LogNormalize", scale.factor = 10000)
object1 <- FindVariableFeatures(object1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(object1)
object1 <- ScaleData(object1, features = all.genes)
object1.markers <- FindAllMarkers(object1, only.pos = TRUE, min.pct = min_pct, logfc.threshold = logfc)
top10<-object1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

heatmap_name<-paste(name,"_heatmap_top10_DE.pdf",sep="")
DE_table_name<-paste(name,min_pct,logfc,"DE_table.csv",sep="_")

pdf(heatmap_name,30,20)
DoHeatmap(object1, features = top10$gene)
dev.off()

write.csv(object1.markers,DE_table_name)

file.move(heatmap_name,folder_name,overwrite = TRUE)
file.move(DE_table_name,folder_name,overwrite = TRUE)


