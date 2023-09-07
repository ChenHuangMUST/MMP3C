###single cell clustering and annotation
###single cell package
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(SingleR)
library(R.utils)
###
message("time start at ",start_time)
start_time <- Sys.time()
sce<-CreateSeuratObject(counts = imputed_data, project = "imputed_sc")
###
sce@meta.data$log10GenesPerUMI <- log10(sce$nFeature_RNA) / log10(sce$nCount_RNA)
sce$mitoRatio <- PercentageFeatureSet(object = sce, pattern = "^MT-")
sce$mitoRatio <- sce@meta.data$mitoRatio / 100
###
metadata <- sce@meta.data
metadata$cells <- rownames(metadata)
metadata$sample <- substring(metadata$cells,18,19)
sce@meta.data <- metadata
###
sce <- subset(x = sce, 
              subset= (nCount_RNA >= 500) &
                (nCount_RNA <= 50000) &
                (nCount_RNA >= 250) &
                (nFeature_RNA <= 9000) &
                (log10GenesPerUMI > 0.80) & 
                (mitoRatio < 0.20))
####
save(sce,file = "D:/tcga/sc/qc/filter_seurat.Rdata")
###
sce <- NormalizeData(sce, normalization.method =  "LogNormalize",
                     scale.factor = 10000)
GetAssay(sce,assay = "RNA")
###
sce <- FindVariableFeatures(sce,
                            selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(sce), 10)
###
sce <- ScaleData(sce)
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce))
###
DimHeatmap(sce, dims = 1:12, cells = 100, balanced = TRUE)
ElbowPlot(sce)
###
sce <- JackStraw(sce, num.replicate = 100)
sce <- ScoreJackStraw(sce, dims = 1:20)
JackStrawPlot(sce, dims = 1:20)
###
sce <- FindNeighbors(sce, dims = 1:18)
library(clustree)
clustree(sce,prefix = 'res.') + coord_flip()
seq <- seq(0.1,1,by = 0.1)

for (res in seq) {
  sci <- FindClusters(sci, resolution = res)
}
clustree(sce,prefix = 'RNA_snn_res.') + coord_flip()
###
sci <- FindClusters(sce, resolution = 0.5)
sci <- RunTSNE(object = sci, dims = 1:18, do.fast = TRUE)
DimPlot(sci,reduction = "tsne",label=T)
###
sce.markers <- FindAllMarkers(object = sci, only.pos = TRUE, min.pct = 0.25,
                              thresh.use = 0.25)

end_time <- Sys.time()
message("time again ",end_time)
####
break_time2 <- Sys.time()
message("time again ",break_time)
###output file
library(data.table)
fwrite(x = as.data.frame(filtered_seurat[["RNA"]]@counts), row.names=T,file = "counts.csv")
###
sce.markers$celltype<-celltype$anno[match(sce.markers$cluster,celltype$ClusterID)]

top5 <- sce.markers %>% group_by(celltype) %>% top_n(5, avg_log2FC)

p <- DotPlot(sce, features = unique(top5$gene) ,
             assay='RNA' )  

p+ theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))
###

head(top5)
top5=top5[!duplicated(top5$gene),]
select_genes_all=split(top5$gene,top5$celltype)
select_genes_all
DotPlot(object = sce, 
        features=select_genes_all, 
        assay = "RNA") + theme(axis.text.x = element_text(angle = 45, 
                                                          vjust = 0.5, hjust=0.5))
features<-c("CD3D","CD3E","CD3G","MS4A1","CD79A","CD79B","KRT18","KRT19","EPCAM","C1QA","C1QB",
            "CD68","DCN","LUM","COL1A2","CLDN5","VWF","PECAM1","NKG7","GNLY","KLRD1")

DotPlot(object = sce, 
        features= features, 
        assay = "RNA",group.by = "celltype") + theme(axis.text.x = element_text(angle = 45, 
                                                                                vjust = 0.5, hjust=0.5))
select_genes = split(features, rep(unique(top5$celltype)[c(2,1,5,6,7,3,4)],each = 3 ))

DotPlot(object = sce, 
        features= select_genes, 
        assay = "RNA",group.by = "celltype") + theme(axis.text.x = element_text(angle = 45, 
                                                                                vjust = 0.5, hjust=0.5))
##################
#infer tumor cells
###### ########### 
library(infercnv)
library(AnnoProbe)
library(Seurat)
###gene annotation file
geneInfor=annoGene(rownames(sce),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
####
b.cells  <- row.names(sce@meta.data)[which(sce$celltype=='B_cell')]
fib.cells  <-  row.names(sce@meta.data)[sce$celltype == 'Fibroblast']
epi.cells  <-  row.names(sce@meta.data)[sce$celltype == 'Epithelial_cells']
###
geneInfor2<-geneInfor
geneInfor2$chr<-substring(geneInfor2$chr,4)
geneInfor2$chr<-as.numeric(geneInfor2$chr)
geneInfor2=geneInfor2[with(geneInfor2, order(chr, start)),]
geneInfor2$chr<-paste("chr",geneInfor2$chr,sep = "")
###sample 800 cells
fib.cells=sample(fib.cells,800)
endo.cells=sample(endo.cells,800)
###
epiMat=as.data.frame(GetAssayData(subset(sce, cells=epi.cells),slot = 'counts'))
fibMat=as.data.frame(GetAssayData(subset(sce, cells=fib.cells),slot = 'counts'))
bMat=as.data.frame(GetAssayData(subset(sce, cells=b.cells),slot = 'counts'))
###cind data
dat=cbind(epiMat,fibMat,bMat)
###reorder dat
dat<-dat[geneInfor2$SYMBOL,]
dat=dat[match( geneInfor2$SYMBOL, rownames(dat) ),] 
###group information
groupinfo=data.frame(v1=colnames(dat),
                     v2=c(rep('epi',ncol(epiMat)),
                          rep('spike-fib',500),
                          rep('ref-fib',626),
                          rep('spike-b',1000),
                          rep('ref-b',6239)))

###reorder gene pos
geneinfo$chr[geneinfo$chr %in% chr_exclude]<-"chr24"
geneinfo$chr<-as.numeric(substring(geneinfo$chr,4))
geneinfo=geneinfo[with(geneinfo, order(chr, start)),]
geneinfo$chr<-paste("chr",geneinfo$chr,sep = "")
dat=dat[match( row.names(geneinfo), rownames(dat) ),] 
row.names(groupinfo)<-groupinfo$v1
groupinfo<-groupinfo[,-1,drop = F]
###input data
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(dat),
                                    annotations_file=groupinfo,
                                    delim="\t",
                                    gene_order_file= geneinfo,
                                    ref_group_names=c("ref-b","ref-fib"),
                                    chr_exclude = "chr24")  neinfo$chr))
dat<-as.matrix(dat)
row.names(geneinfo)<-geneinfo$SYMBOL
geneinfo<-geneinfo[,-1,drop = F]

slot(infercnv_obj, "gene_order")[,"chr"] <- factor(slot(infercnv_obj, "gene_order")[,"chr"], 
                                                   levels = c("chr1", "chr2","chr3","chr4", "chr5", "chr6","chr7", "chr8", "chr9","chr10", "chr11", "chr12","chr13", "chr14", "chr15","chr16", "chr17", "chr18","chr19", "chr20", "chr21","chr22"))

###run
dir.create("infercnv2")
infercnv_all = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir= "infercnv5",  # dir is auto-created for storing outputs
                             cluster_by_groups=T, num_threads= 20,
                             denoise = T,
                             HMM = F,
                             output_format = "pdf",
                             plot_steps = F)

###
#monocle
###
library(Seurat)
library(monocle)
###delete unknowned plasma and NK cell
s.cell<-s.cell[c(1:4,6:7)]
tcells<-subset(t_cells,celltype %in% s.cell)
###
dir.create(path = "D:/tcga/sc/qc2/03_monocle/tcells")
save(tcells,file = "D:/tcga/sc/qc2/03_monocle/tcells/tcell_mp_psertome.Rdata")
###
cds <- newCellDataSet(datx,
                      phenoData = pd,
                      featureData = fd)
###
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds)
###
disp_table<-dispersionTable(cds)
disp.genes <-subset(disp_table,mean_expression >=0.1&dispersion_empirical>= 1*dispersion_fit)$gene_id
cds<-setOrderingFilter(cds,disp.genes)
plot_ordering_genes(cds)
cds<-reduceDimension(cds,max_components = 2, reduction_method = "DDRTree")
cds<-orderCells(cds)
###
plot_cell_trajectory(cds,color_by = "Pseudotime" ,show_backbone = T)
plot_cell_trajectory(cds, color_by = "seurat_clusters")
plot_cell_trajectory(cds,color_by = "Pseudotime" ,show_backbone = T)+facet_wrap("~cell")
###
df<-pData(cds)
ggplot(df, aes(Pseudotime, color = celltype, fill = celltype))+
  geom_density(bw = 0.5, linewidth = 1, alpha = 0.5)+
  theme_classic2()
