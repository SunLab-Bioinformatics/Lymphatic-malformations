library(Seurat)
library(tidyverse)
library(cowplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(qs)

setwd('~/LM/merge/')
outdir <- "./files/"
dir.create(outdir,recursive = T)
outdir2 <- "./plots/"
dir.create(outdir2,recursive = T)

#### 1.Load data ----
counts1 <-Read10X(data.dir = "./data/BFYY331/")
seurat_obj1 <- CreateSeuratObject(
  counts = counts1,
  project = "BFYY331",
  min.cells = 3,
  min.features = 200
)

counts2 <- Read10X(data.dir = "./data/BLJJ329/")
seurat_obj2 <- CreateSeuratObject(
  counts = counts2,
  project = "BLJJ329",
  min.cells = 3,
  min.features = 200
)

counts3 <- Read10X(data.dir = "./data/BHLJ1022/")
seurat_obj3 <- CreateSeuratObject(
  counts = counts3,
  project = "BHLJ1022",
  min.cells = 3,
  min.features = 200
)

counts4 <- Read10X(data.dir = "./data/BJZY1103/")
seurat_obj4 <- CreateSeuratObject(
  counts = counts4,
  project = "BJZY1103",
  min.cells = 3,
  min.features = 200
)

counts5 <- Read10X(data.dir = "./data/BMYY203/")
seurat_obj5 <- CreateSeuratObject(
  counts = counts5,
  project = "BMYY203",
  min.cells = 3,
  min.features = 200
)

counts6 <- Read10X(data.dir = "./data/CNLHL420/")
seurat_obj6 <- CreateSeuratObject(
  counts = counts6,
  project = "CNLHL420",
  min.cells = 3,
  min.features = 200
)

counts7 <- Read10X(data.dir = "./data/BLBL804/")
seurat_obj7 <- CreateSeuratObject(
  counts = counts7,
  project = "BLBL804",
  min.cells = 3,
  min.features = 200
)

counts8 <- Read10X(data.dir = "./data/BLJK804/")
seurat_obj8 <- CreateSeuratObject(
  counts = counts8,
  project = "BLJK804",
  min.cells = 3,
  min.features = 200
)

counts9 <- Read10X(data.dir = "./data/CFYY331/")
seurat_obj9 <- CreateSeuratObject(
  counts = counts9,
  project = "CFYY331",
  min.cells = 3,
  min.features = 200
)

counts10 <- Read10X(data.dir = "./data/CLJJ329/")
seurat_obj10 <- CreateSeuratObject(
  counts = counts10,
  project = "CLJJ329",
  min.cells = 3,
  min.features = 200
)

counts11 <- Read10X(data.dir = "./data/CGCJ1028/")
seurat_obj11 <- CreateSeuratObject(
  counts = counts11,
  project = "CGCJ1028",
  min.cells = 3,
  min.features = 200
) 

##
merge <- merge(seurat_obj1,
               y=c(seurat_obj2,seurat_obj3,seurat_obj4,seurat_obj5,
                   seurat_obj6,seurat_obj7,seurat_obj8,
                   seurat_obj9,seurat_obj10,seurat_obj11),
               project = "data")

#### 2.metadata ----
merge$orig.ident <- factor(merge$orig.ident, levels = c("CNLHL420", "BLBL804","BLJK804",
                                                        "BFYY331","BLJJ329","BHLJ1022","BJZY1103","BMYY203",
                                                        "CFYY331", "CLJJ329", "CGCJ1028"))
merge$donor <- ''
merge$donor[merge$orig.ident%in% c("CFYY331", "BFYY331")] <- "FYY"
merge$donor[merge$orig.ident%in% c("CLJJ329", "BLJJ329")] <- "LJJ"
merge$donor[merge$orig.ident%in% c("CNLHL420")] <- "LHL"
merge$donor[merge$orig.ident%in% c("BLBL804")] <- "LBL"
merge$donor[merge$orig.ident%in% c("BLJK804")] <- "LJK"
merge$donor[merge$orig.ident%in% c("BHLJ1022")] <- "HLJ"
merge$donor[merge$orig.ident%in% c("BJZY1103")] <- "JZY"
merge$donor[merge$orig.ident%in% c("BMYY203")] <- "MYY"
merge$donor[merge$orig.ident%in% c("CGCJ1028")] <- "GCJ"

merge$group <- ""
merge$group[merge$orig.ident%in% c("BFYY331","BLJJ329","BHLJ1022","BJZY1103","BMYY203")] <- "LM"
merge$group[merge$orig.ident%in% c("CNLHL420", "BLBL804","BLJK804")] <- "Normal"
merge$group[merge$orig.ident%in% c("CFYY331", "CLJJ329", "CGCJ1028")] <- "PE"
merge$group <- factor(merge$group, levels = c("Normal","LM","PE"))

merge <- PercentageFeatureSet(merge,
                              pattern = "^MT-",
                              col.name = "percent.mt")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merge <- CellCycleScoring(merge, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

##
genes_list <- read.csv('~/Collection/files/mt_rb_hb_hs_dig_list.csv')
merge = AddModuleScore(merge, features = list(genes_list$rb[nchar(genes_list$rb)>0]), name = 'rb_score',seed = 2023)
merge = AddModuleScore(merge, features = list(genes_list$hs[nchar(genes_list$hs)>0]), name = 'hs_score',seed = 2023)
merge = AddModuleScore(merge, features = list(genes_list$hb[nchar(genes_list$hb)>0]), name = 'hb_score',seed = 2023)
merge = AddModuleScore(merge, features = list(genes_list$dis[nchar(genes_list$dis)>0]), name = 'dis_score',seed = 2023)
colnames(merge@meta.data) <- gsub('e1', 'e', colnames(merge@meta.data))
VlnPlot(merge,group.by = 'orig.ident',
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "rb_score","hs_score", 
                     "hb_score",'dis_score','S.Score','G2M.Score'),
        ncol = 3,pt.size = 0)
ggsave(paste0(outdir2,'merged_all_index.pdf'),
       last_plot(), height= 10, width= 10)
Idents(merge) <- factor(Idents(merge), levels = c("CNLHL420", "BLBL804","BLJK804",
                                                  "BFYY331","BLJJ329","BHLJ1022","BJZY1103","BMYY203",
                                                  "CFYY331", "CLJJ329", "CGCJ1028"))
saveRDS(merge,paste0(outdir,'merged_all_ori.rds'))

merge1 <- subset(merge, subset = percent.mt <= 7 & nFeature_RNA >= 500 & nFeature_RNA <= 3500 )
VlnPlot(merge1,group.by = 'orig.ident',
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "rb_score","hs_score", 
                     "hb_score",'dis_score','S.Score','G2M.Score'),
        ncol = 3,pt.size = 0)
ggsave(paste0(outdir2,'merged_cut_index.pdf'),
       last_plot(), height= 10, width= 10)
qs::qsave(merge1,paste0(outdir, 'merged_cut_ori.qs'))

#### 3.cluster ----
set.seed(2023)
case <- 'CCA-SCT_'

seurat_obj <- qs::qread('./files/merged_cut_ori.qs')
seurat_obj<-NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
split_seurat <- SplitObject(seurat_obj, split.by = "orig.ident")
split_seurat <- lapply(
  X = split_seurat,
  FUN = function(x) {
    x <- SCTransform(
      x,
      vars.to.regress = c("percent.mt"),
      vst.flavor = "v2",
      verbose = T
    )
  }
)
features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 2000)
write.table(sort(features), file.path(outdir, paste0(case, "2000_features.txt")),
  quote = F, col.names = F, row.names = F)

exclude_genes <- read.table('~/Collection/files/exclude_genes.txt') |> pull()
HVGs <- features[!(features %in% exclude_genes)]
immunoglobulin_genes <- features[grep("^IG[HKL]", features)]
HVGs <- HVGs[!(HVGs %in% immunoglobulin_genes)]
write.table(sort(HVGs), file.path(outdir, paste0(case, "HVGs.txt")),
  quote = F, col.names = F, row.names = F)

split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = HVGs)
immune.anchors <- FindIntegrationAnchors(object.list = split_seurat, normalization.method = "SCT",
                                         anchor.features = HVGs, dims = 1:20)
seurat_obj <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
DefaultAssay(seurat_obj) <- "integrated"
seurat_obj <- seurat_obj |>
  RunPCA() |>
  RunTSNE(dims = 1:20) |>
  RunUMAP(dims = 1:20) |>
  FindNeighbors(dims = 1:20) |>
  FindClusters(resolution = 1.5)
seurat_obj$seurat_clusters <- factor(seurat_obj$seurat_clusters,
                                     levels=c(0 : length(table(seurat_obj$seurat_clusters))-1 ))
DefaultAssay(seurat_obj) <- 'RNA'
Idents(seurat_obj) <- 'seurat_clusters'
qs::qsave(seurat_obj, file.path(outdir,paste0(case, "20_2000_15.qs")))

#### 4.name ----
seurat_obj <- qs::qread('./files/seurat_obj_20_2000_15.qs')
DimPlot(seurat_obj, cols = type_colors,
        label = T, label.size = 3,label.box = T, repel = T,
        reduction = "umap") + theme_bw() +
  labs( x= "UMAP 1",y= "UMAP 2") +
  theme(panel.grid = element_blank(),aspect.ratio = 1)+NoLegend()

case = 'named_'
cell_type <- data.frame(ClusterID = 0:39,
                        cell_type = 'Unknown')
cell_type[cell_type$ClusterID %in% c(4,11,13,22,26,15), 2] = 'Monocytes'
cell_type[cell_type$ClusterID %in% c(3,8,6,7,9,18,23,24), 2] = 'CD4+ T cells'
cell_type[cell_type$ClusterID %in% c(2,10,28,37,21,12,16,17), 2] = 'CD8+ T cells'
cell_type[cell_type$ClusterID %in% c(1,14,27), 2] = 'NK cells'
cell_type[cell_type$ClusterID %in% c(20), 2] = 'MAIT'
cell_type[cell_type$ClusterID %in% c(19), 2] = 'Gamma delta T cells'
cell_type[cell_type$ClusterID %in% c(0,5,32,34,38,39), 2] = 'B cells'
cell_type[cell_type$ClusterID %in% c(29), 2] = 'Dendritic cells'
cell_type[cell_type$ClusterID %in% c(36), 2] = 'pDC'
cell_type[cell_type$ClusterID %in% c(31), 2] = 'Neutrophils'
cell_type[cell_type$ClusterID %in% c(25), 2] = 'Macrophages'
cell_type[cell_type$ClusterID %in% c(30,33), 2] = 'Platelets'
seurat_obj@meta.data$cell_type <- "NA"
for (i in 1:nrow(cell_type)) {
  seurat_obj@meta.data[which(seurat_obj@meta.data$seurat_clusters == cell_type$ClusterID[i]), 'cell_type'] <-
    cell_type$cell_type[i]
}

seurat_obj$orig.ident <- factor(seurat_obj$orig.ident, 
                                levels = c("BLBL804" ,"BLJK804" ,"CNLHL420",
                                           "BFYY331" ,"BLJJ329","BHLJ1022","BJZY1103","BMYY203",
                                           "CFYY331" ,"CLJJ329", "CGCJ1028" ))
seurat_obj$donor <- ""
seurat_obj$donor[seurat_obj$orig.ident%in% c('CNLHL420')] <- "N1"
seurat_obj$donor[seurat_obj$orig.ident%in% c('BLBL804')] <- "N2"
seurat_obj$donor[seurat_obj$orig.ident%in% c('BLJK804')] <- "N3"
seurat_obj$donor[seurat_obj$orig.ident%in% c('BFYY331')] <- "LM1"
seurat_obj$donor[seurat_obj$orig.ident%in% c('BLJJ329')] <- "LM2"
seurat_obj$donor[seurat_obj$orig.ident%in% c('BHLJ1022')] <- "LM3"
seurat_obj$donor[seurat_obj$orig.ident%in% c('BJZY1103')] <- "LM4"
seurat_obj$donor[seurat_obj$orig.ident%in% c('BMYY203')] <- "LM5"
seurat_obj$donor[seurat_obj$orig.ident%in% c('CFYY331')] <- "PE1"
seurat_obj$donor[seurat_obj$orig.ident%in% c('CLJJ329')] <- "PE2"
seurat_obj$donor[seurat_obj$orig.ident%in% c('CGCJ1028')] <- "PE3"
seurat_obj$donor <- factor(seurat_obj$donor, 
                           levels = c('N1',  'N2', 'N3',
                                      'LM1',  'LM2', 'LM3', 'LM4',  
                                      'PE1', 'PE2', 'PE3' ))
seurat_obj$group <- factor(seurat_obj$group, 
                           levels = c("Normal" ,"LM" ,"PE"))
seurat_obj$cell_type <- factor(seurat_obj$cell_type, 
                               levels = c('B cells',  'CD4+ T cells', 'CD8+ T cells', 'MAIT', 'Gamma delta T cells', 'NK cells',
                                          'Monocytes', 'Macrophages', 'Dendritic cells', 'Neutrophils', 'pDC',  'Platelets',  'Unknown'))
Idents(seurat_obj) <- 'cell_type'
qs::qsave(seurat_obj, file.path(outdir, paste0(case, "seurat.qs")))