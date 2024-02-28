#### 0.前期配置 ----
library(Seurat)
library(tidyverse)
#library(RCurl)
library(cowplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(qs)
#library(harmony)
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

setwd('~/LM/project2305/merge/results/230430/')
outdir <- "./files/"
dir.create(outdir,recursive = T)
outdir2 <- "./plots/"
dir.create(outdir2,recursive = T)
case <- 'merged_all_'

#### 1.读入数据----
#### LM
counts1 <-
  Read10X(data.dir = "~/LM/project2305/data/order/RNA/BFYY331/")
seurat_obj1 <- CreateSeuratObject(
  counts = counts1,
  project = "BFYY331",
  min.cells = 3,
  min.features = 200
) 
head(seurat_obj1@meta.data)

counts2 <-
  Read10X(data.dir = "~/LM/project2305/data/order/RNA/BLJJ329/")
seurat_obj2 <- CreateSeuratObject(
  counts = counts2,
  project = "BLJJ329",
  min.cells = 3,
  min.features = 200
)
head(seurat_obj2@meta.data)

counts3 <-
  Read10X(data.dir = "~/LM/project2305/data/order/RNA/BHLJ1022/")
seurat_obj3 <- CreateSeuratObject(
  counts = counts3,
  project = "BHLJ1022",
  min.cells = 3,
  min.features = 200
)
head(seurat_obj3@meta.data)

counts4 <-
  Read10X(data.dir = "~/LM/project2305/data/order/RNA/BJZY1103/")
seurat_obj4 <- CreateSeuratObject(
  counts = counts4,
  project = "BJZY1103",
  min.cells = 3,
  min.features = 200
) 
head(seurat_obj4@meta.data)

counts5 <-
  Read10X(data.dir = "~/LM/project2305/data/order/RNA/BMYY203/")
seurat_obj5 <- CreateSeuratObject(
  counts = counts5,
  project = "BMYY203",
  min.cells = 3,
  min.features = 200
) 
head(seurat_obj5@meta.data)

####Normal
counts6 <-
  Read10X(data.dir = "~/LM/project2305/data/order/RNA/CNLHL420/")
seurat_obj6 <- CreateSeuratObject(
  counts = counts6,
  project = "CNLHL420",
  min.cells = 3,
  min.features = 200
) 
head(seurat_obj6@meta.data)

counts7 <-
  Read10X(data.dir = "~/LM/project2305/data/order/RNA/BLBL804/")
seurat_obj7 <- CreateSeuratObject(
  counts = counts7,
  project = "BLBL804",
  min.cells = 3,
  min.features = 200
) 
head(seurat_obj7@meta.data)

counts8 <-
  Read10X(data.dir = "~/LM/project2305/data/order/RNA/BLJK804/")
seurat_obj8 <- CreateSeuratObject(
  counts = counts8,
  project = "BLJK804",
  min.cells = 3,
  min.features = 200
) 
head(seurat_obj8@meta.data)

####PE
counts9 <-
  Read10X(data.dir = "~/LM/project2305/data/order/RNA/CFYY331/")
seurat_obj9 <- CreateSeuratObject(
  counts = counts9,
  project = "CFYY331",
  min.cells = 3,
  min.features = 200
) 
head(seurat_obj9@meta.data)

counts10 <-
  Read10X(data.dir = "~/LM/project2305/data/order/RNA/CLJJ329/")
seurat_obj10 <- CreateSeuratObject(
  counts = counts10,
  project = "CLJJ329",
  min.cells = 3,
  min.features = 200
) 
head(seurat_obj10@meta.data)

counts11 <-
  Read10X(data.dir = "~/LM/project2305/data/order/RNA/CGCJ1028/")
seurat_obj11 <- CreateSeuratObject(
  counts = counts11,
  project = "CGCJ1028",
  min.cells = 3,
  min.features = 200
) 
head(seurat_obj11@meta.data)


#### 2.合并数据 ----
merge <- merge(seurat_obj1,
               y=c(seurat_obj2,seurat_obj3,seurat_obj4,seurat_obj5,
                   seurat_obj6,seurat_obj7,seurat_obj8,
                   seurat_obj9,seurat_obj10,seurat_obj11),
               project = "data11_2305")

#### 3.添加分组 ----
merge$orig.ident <- factor(merge$orig.ident, levels = c("CNLHL420", "BLBL804","BLJK804",
                                                        "BFYY331","BLJJ329","BHLJ1022","BJZY1103","BMYY203",
                                                        "CFYY331", "CLJJ329", "CGCJ1028"))

merge$group <- ""
merge$group[merge$orig.ident%in% c("BFYY331","BLJJ329","BHLJ1022","BJZY1103","BMYY203")] <- "LM"
merge$group[merge$orig.ident%in% c("CNLHL420", "BLBL804","BLJK804")] <- "Normal"
merge$group[merge$orig.ident%in% c("CFYY331", "CLJJ329", "CGCJ1028")] <- "PE"
merge$group <- factor(merge$group, levels = c("Normal","LM","PE"))

#质量控制参数
merge <- PercentageFeatureSet(merge,
                                   pattern = "^MT-",
                                   col.name = "percent.mt")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merge <-
  CellCycleScoring(
    merge,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = FALSE)

# 打分
genes_list <- read.csv('~/heart/202303/files/mt_rb_hb_hs_dig_list.csv')
merge = AddModuleScore(merge, features = list(genes_list$rb[nchar(genes_list$rb)>0]), name = 'rb_score',seed = 2023)
merge = AddModuleScore(merge, features = list(genes_list$hs[nchar(genes_list$hs)>0]), name = 'hs_score',seed = 2023)
merge = AddModuleScore(merge, features = list(genes_list$hb[nchar(genes_list$hb)>0]), name = 'hb_score',seed = 2023)
merge = AddModuleScore(merge, features = list(genes_list$dis[nchar(genes_list$dis)>0]), name = 'dis_score',seed = 2023)
colnames(merge@meta.data) <- gsub('e1', 'e', colnames(merge@meta.data))
#merge@meta.data$MALAT1 <- merge@assays$RNA@counts['MALAT1',]
VlnPlot(merge,group.by = 'orig.ident',
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,pt.size = 0)
VlnPlot(merge,group.by = 'orig.ident',
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "rb_score","hs_score", 
               "hb_score",'dis_score','S.Score','G2M.Score'),
  ncol = 3,pt.size = 0)
 
ggsave(file.path(outdir2,paste0(case,'指标.pdf')),
       ggplot2::last_plot(),
       height= 10,
       width= 10)
Idents(merge) <- factor(Idents(merge), levels = c("CNLHL420", "BLBL804","BLJK804",
                                                  "BFYY331","BLJJ329","BHLJ1022","BJZY1103","BMYY203",
                                                  "CFYY331", "CLJJ329", "CGCJ1028"))
#保存
qs::qsave(merge,file.path(outdir, paste0(case, 'ori.qs')))
saveRDS(merge,file.path(outdir, paste0(case, 'ori.rds')))

#### 4.质量过滤 ----
case <- 'merged_cut_'
merge1 <- subset(merge, subset =
                       percent.mt <= 7 & 
                       nFeature_RNA >= 500 & 
                       nFeature_RNA <= 3500 )
dim(merge1) 
table(merge1$orig.ident)

VlnPlot(merge1,group.by = 'orig.ident',
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "rb_score","hs_score", 
                     "hb_score",'dis_score','S.Score','G2M.Score'),
        ncol = 3,pt.size = 0)
ggsave(file.path(outdir2,paste0(case,'指标.pdf')),
       ggplot2::last_plot(),
       height= 10,
       width= 10)

qs::qsave(merge1,file.path(outdir, paste0(case, 'ori.qs')))

#### 5.降维聚类 ----
merge2 <- merge1 
case <- 'null_'

merge2<-NormalizeData(merge2, normalization.method = "LogNormalize", scale.factor = 10000)
merge2 <- FindVariableFeatures(merge2, selection.method = "vst", nfeatures = 2000)
merge2 <- ScaleData(merge2)
merge2 <- RunPCA(merge2)
ElbowPlot(merge2, ndims = 50)
ggsave(file.path(outdir2,paste0(case,'ElbowPlot.pdf')),
       ggplot2::last_plot(),
       height= 5,
       width= 5)
set.seed(2023)
merge2 <- merge2 |>
  RunTSNE(dims = 1:20) |>
  RunUMAP(dims = 1:20) |>
  FindNeighbors(dims = 1:20) |>
  FindClusters(resolution = 0.5)
qs::qsave(merge2, file.path(outdir,paste0(case, "20_05.qs")))
features <- merge2@assays$RNA@var.features
HVGs <- features
write.table(
  sort(HVGs),
  file.path(outdir, paste0(case, "20000HVGs.txt")),
  quote = F,
  col.names = F,
  row.names = F)

#### 6.查看结果 ----
seurat_obj <- merge2
DimPlot(seurat_obj, cols = my36colors,  
        label = T, label.size = 3,label.box = T, repel = T,
        reduction = "umap", group.by = "seurat_clusters",) + theme_bw() +
  labs( x= "UMAP 1",y= "UMAP 2",title = "seurat_clusters") +
  theme(panel.grid = element_blank(),aspect.ratio = 1)
ggsave(file.path(outdir2,paste0(case,'seurat_clusters.pdf')),
       ggplot2::last_plot(),
       height= 5,
       width= 5)

VlnPlot(seurat_obj,group.by = 'seurat_clusters',
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "rb_score","hs_score", 
                     "hb_score",'dis_score','S.Score','G2M.Score'),
        ncol = 3,pt.size = 0)
ggsave(file.path(outdir2,paste0(case,'指标.pdf')),
       ggplot2::last_plot(),
       height= 10,
       width= 20)

#data
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident",
              cols = my36colors) + theme_bw() +
  labs( x= "UMAP 1",y= "UMAP 2",title = "orig.ident") +
  theme(panel.grid = element_blank(),aspect.ratio = 1)

p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "donor",
              cols = my36colors) + theme_bw() +
  labs( x= "UMAP 1",y= "UMAP 2",title = "donor") +
  theme(panel.grid = element_blank(),aspect.ratio = 1)

p3 <-  DimPlot(seurat_obj, reduction = "umap", group.by = "group",
               cols = my36colors) + theme_bw() +
  labs( x= "UMAP 1",y= "UMAP 2",title = "group") +
  theme(panel.grid = element_blank(),aspect.ratio = 1)

p4 <- DimPlot(seurat_obj, cols = my36colors,
              reduction = "umap", group.by = "seurat_clusters",) + theme_bw() +
  labs( x= "UMAP 1",y= "UMAP 2",title = "seurat_clusters") +
  theme(panel.grid = element_blank(),aspect.ratio = 1)

#score
p5 <- DimPlot(seurat_obj, cols = my36colors,
              reduction = "umap", group.by = "seurat_clusters",) + theme_bw() +
  labs( x= "UMAP 1",y= "UMAP 2",title = "seurat_clusters") +
  theme(panel.grid = element_blank(),aspect.ratio = 1)

p6 <- FeaturePlot(seurat_obj,features= 'hs_score', max.cutoff = 5)+
  scale_colour_distiller(palette = "YlOrRd",direction = 1)+
  #guides(color=guide_colorbar(title = 'z score'))+
  theme(aspect.ratio = 1,axis.title = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = -5))

p7 <- FeaturePlot(seurat_obj,features= 'rna_MALAT1',
                  slot = 'data',min.cutoff = 3,max.cutoff = 8)+
  scale_colour_distiller(palette = "YlOrRd",direction = 1)+
  #guides(color=guide_colorbar(title = 'z score'))+
  theme(aspect.ratio = 1,axis.title = element_text(size = 10),
        plot.title = element_text(size = 10,face = "italic",vjust = -5))

p8 <- FeaturePlot(seurat_obj,features= 'dis_score',max.cutoff = 8)+
  scale_colour_distiller(palette = "YlOrRd",direction = 1)+
  #guides(color=guide_colorbar(title = 'z score'))+
  theme(aspect.ratio = 1,axis.title = element_text(size = 10),
        plot.title = element_text(size = 10,vjust = -5))



p <-p1 + p2 + p3 + p4 +plot_layout(nrow = 2, ncol = 2,guides = "keep")
p
ggsave(
  file.path(outdir2, paste0(case, "umap_data.pdf")),
  plot = p,
  height = 8,
  width = 12
)

p <-p5 + p6 + p7 + p8 +plot_layout(nrow = 2, ncol = 2,guides = "keep")
ggsave(
  file.path(outdir2, paste0(case, "umap_score.pdf")),
  plot = p,
  height = 8,
  width = 12
)


#### 7.lisi ----
lisi_score <- lisi::compute_lisi(
  X = Embeddings(seurat_obj, reduction = "pca"),
  meta_data = seurat_obj@meta.data,
  label_colnames =  "orig.ident"
) |>
  dplyr::mutate(type = "null") |>
  rownames_to_column("ID")
head(lisi_score)
round(mean(lisi_score$orig.ident), 2) #2.04
write_csv(lisi_score, file = file.path(outdir, paste0(case, "lisi_score.csv")))

#### 8.marker ----
case <- 'pan_markers_'
outdir2 <- './plots/markers/'
dir.create(outdir2, recursive = T)

features <- c('MS4A1','CD79B','CD79A' ,#2 BC
              'LEF1','TCF7','CCR7',#1 CD4_Naive #3 CD8_Naive
              'CD4',"CD40LG",#CD4+T,
              "CD8B",'CD8A',#CD8+
              'TRGV9','TRDV2',#8 γδT
              'GNLY','KLRD1','NKG7',#4 NK
              'CD14','VCAN','S100A8',#0 Mono
              'FCGR3A','CD68','CSF1R',#12 MAC
              "CD1C","CLEC10A",'HLA-DQA1',#14 DC
              'JCHAIN','IGHA1' ,'MZB1', #16,17 'Plamsma'
              "PPBP", "PF4", "GP1BB", #13 'Platelets'
              "PPBP","CCL5","CLDN5","PLA2G12A",#megakaryocytes
              "CXCL8","FCGR3B","MNDA","CD44",#Neutrophils,
              "LILRA4","IL3RA","TCF4","CLEC4C",#pDC
              "STMN1","HMGB2","HIST1H4C"#Cycle B cells
)
p <- DotPlot(seurat_obj, features = features[!duplicated(features)])  + coord_flip()+ 
  scale_colour_distiller(palette = "YlOrRd", direction = 1)+ 
  theme(axis.text.x = element_text(angle = 0, hjust =0.5),
        axis.text.y = element_text(face = 'italic'),
        axis.title = element_blank())
p
ggsave(file.path(outdir2,paste0(case,'blood_Dotplot.pdf')),
       p,
       height=10,
       width=12)

features <- c('MS4A1','CD79B', 
                      'CD4', 'CD8B',
                      'GNLY','NKG7',
                      'TRGV9','TRDV2',
                      'CD14','VCAN',
                      'FCGR3A',"CD1C",
                      'JCHAIN','IGHA1' ,
                      "PPBP", "PF4",
                      "CXCL8","FCGR3B",
                      "LILRA4","IL3RA",
                      "STMN1","HMGB2")
features <- features[features %in% rownames(seurat_obj)]

for(i in 1:length(features)){
  assign(paste0("p", i), FeaturePlot(seurat_obj,features= features[i],
                                     min.cutoff = 0,max.cutoff = 3 )+
           scale_colour_distiller(palette = "YlOrRd",direction = 1)+theme_void()+
           theme(aspect.ratio = 1,plot.title = element_text(size = 10,face = "italic",
                                                            hjust = 0.5, vjust = -1)))
}

p <-p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9+ p10+
  plot_layout(nrow = 2, ncol = 5,guides = "collect")# 
#p
ggsave(
  file.path(outdir2, paste0(case,"FeaturePlot1.pdf")),
  plot = p,
  height = 5,
  width = 15
)

p <-p11 + p12 + p13 + p14 + p15 + p16 + p17 + p18 + p19+ p20+ p21 + p22 +
  plot_layout(nrow = 3, ncol = 5,guides = "collect")# 
#p
ggsave(
  file.path(outdir2, paste0(case,"FeaturePlot2.pdf")),
  plot = p,
  height = 8,
  width = 15
)


#### 9.寻找高变基因----
case <- 'Markers_'
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
top.markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
top.markers2 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
#保存
saveRDS(markers, file.path(outdir, paste0('allmarkers', ".rds")))
saveRDS(top.markers, file.path(outdir, paste0('topmarkers', ".rds")))
saveRDS(top.markers2, file.path(outdir, paste0('topmarkers2', ".rds")))

subobj <- subset(seurat_obj, downsample = 500)
table(subobj$seurat_clusters)
p0 <- DoHeatmap(subobj, features = top.markers$gene, 
                angle = 0,
                group.colors = my36colors,label = F,
                disp.min=-2.5, disp.max=2.5)+
  scale_fill_gradientn(colors = c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
                                  colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)))

pdf(file.path(outdir2, paste0(case, "top20_heatmap.pdf")), width=40, height=20)
print(p0)
dev.off()

# #TOP10
features <- top.markers$gene[!duplicated(top.markers$gene)]
p <- DotPlot(seurat_obj, features = features)  + coord_flip()+
  scale_colour_distiller(palette = "YlOrRd", direction = 1)+
  theme(axis.text.x = element_text(angle = 0, hjust =0.5),
        axis.text.y = element_text(face = 'italic'),
        axis.title = element_blank())
ggsave(file.path(outdir2,paste0(case,'top20_Dotplot.pdf')),
       p,
       height=25,
       width=15)

#TOP5
features <- top.markers2$gene[!duplicated(top.markers2$gene)]
p <- DotPlot(seurat_obj, features = features)  + coord_flip()+ 
  scale_colour_distiller(palette = "YlOrRd", direction = 1)+ 
  theme(axis.text.x = element_text(angle = 0, hjust =0.5),
        axis.text.y = element_text(face = 'italic'),
        axis.title = element_blank())
ggsave(file.path(outdir2,paste0(case,'top5_Dotplot.pdf')),
       p,
       height=15,
       width=15)
#TOP100
subobj <- subset(seurat_obj, downsample = 500)
table(markers$cluster)
features <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) %>%
  pull()

p0 <- DoHeatmap(subobj, features = features,
                group.colors = my36colors,label = F,
                group.bar.height = 0.02, angle = 0,
                disp.min=-2.5, disp.max=2.5)+
  labs(fill = c('z score'))+
  theme(axis.text.y = element_blank())+
  scale_fill_gradientn(colors = c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
                                  colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)))

pdf(file.path(outdir2, paste0(case, "heatmap_100.pdf")), width=25, height=20)
print(p0)
dev.off()

#### 10.细胞群相似性 ----
library(corrplot)
library(pheatmap)
library(AUCell)
case <- 'cor_'
outdir2 <- './plots/cor/'
dir.create(outdir2, recursive = T)

my_palette <- c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
                colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50))

av_expr <-
  AverageExpression(
    seurat_obj,
    assays = "RNA",
    group.by = "seurat_clusters"
  )[["RNA"]] %>%
  as.data.frame()
saveRDS(av_expr, file.path(outdir,paste0(case, 'av_expr_seurat_clusters.rds')))

corr<-cor(av_expr,method="spearman")
saveRDS(corr, file.path(outdir,paste0(case, 'corr_seurat_clusters_all.rds')))

library(pheatmap)
p <- pheatmap(corr, shape = "circle",cellwidth = 10,cellheight =10,
              #breaks = my_breaks, 
              #scale_fill_gradientn = list(midpoint = my_midpoint),
              color = my_palette, main = 'spearman'
)
p
ggsave(file.path(outdir2,paste0(case, 'cluster_all.pdf')),
       p,width = 8,height =8)
#: HVGs
HVGs <- read.table('./files/null_20000HVGs.txt') |> pull()
av_expr <- av_expr[HVGs, ]
corr<-cor(av_expr,method="spearman")
saveRDS(corr, file.path(outdir,'corr_seurat_clusters_HVGs.rds'))

library(pheatmap)
p <- pheatmap(corr, shape = "circle",cellwidth = 10,cellheight =10,
              #breaks = my_breaks, 
              #scale_fill_gradientn = list(midpoint = my_midpoint),
              color = my_palette, main = 'spearman: HVGs'
)
p
ggsave(file.path(outdir2,'cluster_HVGs.pdf'),
       p,width = 8,height =8)
