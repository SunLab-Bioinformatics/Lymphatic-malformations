library(Seurat)
library(ggplot2)
library(GSEABase)
library(msigdbr)
library(dplyr)
library(stringr)
library(AUCell)
library(clusterProfiler)
library(ggraph)
library(ggpubr)
library(tidydr)
library(clusterProfiler)
library(org.Hs.eg.db)
source('~/Collection/code/plot.R')
source('~/Collection/code/scrna-seq.R')

setwd('~/LM/Fig4/')
outdir <- "./files/"
dir.create(outdir,recursive = T)
outdir2 <- "./plots/"
dir.create(outdir2,recursive = T)

seurat_obj <- qs:qread('~/LM/merge/files/impor_seurat.qs')
seurat_obj <- subset(seurat_obj, idents = c('CD4+ T cells', 'CD8+ T cells', 'MAIT', 'Gamma delta T cells'))
qs::qsave(seurat_obj, '~/LM/merge/files/TC.qs')

seurat_obj <- subset(seurat_obj, idents = 'CD8+ T cells')
qs::qsave(seurat_obj, '~/LM/merge/files/CD8+ T cells.qs')

#### 1.Fig4A-B ----
LM_N <- readRDS('./deg/LM_N/CD8+ T cells/deg.rds')
case <- 'LM vs Normal '
## UP
deg <- LM_N[LM_N$change %in% 'UP', ] |> pull('gene')
bp1 <- enrichGO(
  deg,
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL',
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2)
term1 <- bp1@result

up_pathway <- c(
  'positive regulation of innate immune response',
  'immune response-regulating signaling pathway',
  'leukocyte migration',
  'cytokine-mediated signaling pathway',
  'Fc receptor mediated stimulatory signaling pathway',
  'blood coagulation',
  'interleukin-10 production',
  'tumor necrosis factor production'
)
df1 <- term1[term1$Description %in% up_pathway, ]
df1$labelx=rep(0,nrow(df1))
df1$labely=seq(nrow(df1),1)
p <- ggplot(data = df1, aes(x = -log10(pvalue), y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=1, fill= "#ee6470", width = 0.8) +
  geom_text(aes(x=labelx, y=labely,
                label = df1$Description),
            size=3.5, hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12),
        plot.margin=unit(rep(1,4),'lines'))+
  xlab("-log10(pvalue)")+
  ggtitle(case)+
  scale_x_continuous(expand = c(0,0))
ggsave(file.path(outdir2, paste0(case, 'UP_Barplot.pdf')), p, height=3, width=4)
## DOWN
deg <- LM_N[LM_N$change %in% 'DOWN', ] |> pull('gene')
bp2 <-  enrichGO(
  deg,
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL',
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2)
term2 <- bp2@result
down_pathway <- c(
  'T cell differentiation',
  'regulation of T cell activation',
  'T cell selection',
  'thymic T cell selection',
  'T cell receptor signaling pathway',
  'regulation of T cell differentiation',
  'regulation of lymphocyte apoptotic process',
  'regulation of cell-cell adhesion'
)
df2 <- term2[term2$Description %in% down_pathway, ]
df2$labelx=rep(0,nrow(df2))
df2$labely=seq(nrow(df2),1)
p <- ggplot(data = df2,
            aes(x = -log10(pvalue), y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=1, fill= "#00a6e1", width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = df2$Description),
            size=3.5, hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12),
        plot.margin=unit(rep(1,4),'lines'))+
  xlab("-log10(pvalue)")+
  ggtitle(case)+
  scale_x_continuous(expand = c(0,0))
ggsave(file.path(outdir2, paste0(case, 'DOWN_Barplot.pdf')), p, height=3, width=4)

#### 2.Fig4C ----
## subtype
seurat_obj <- qs:qread('~/LM/merge/files/CD8+ T cells.qs')
seurat_obj$orig.ident <- as.character(seurat_obj$orig.ident)
split_seurat <- SplitObject(seurat_obj, split.by = "orig.ident")
split_seurat <- lapply(
  X = split_seurat,
  FUN = function(x) {
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  }
)
features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 2000)
exclude_genes <- read.table('~/Collection/files/exclude_genes.txt') |> pull()
HVGs <- features[!(features %in% exclude_genes)]
TCR_genes <- features[grep("^TR[ABGD]", features)]
HVGs <- HVGs[!(HVGs %in% TCR_genes)][1:1000]
DefaultAssay(seurat_obj) <- "integrated"
seurat_obj@assays$integrated@var.features <- HVGs
set.seed(2023)
seurat_obj <- seurat_obj |>
  RunPCA() |>
  RunTSNE(dims = 1:10) |>
  RunUMAP(dims = 1:10) |>
  FindNeighbors(dims = 1:10) |>
  FindClusters(resolution = 0.5)

DefaultAssay(seurat_obj) <- 'RNA'
Idents(seurat_obj) <- 'seurat_clusters'
features <- c(
  'SELL', 'CCR7', 'CTLA4',
  'CX3CR1','NKG7','GZMB',
  'DUSP2','TIGIT', 'CCL5',
  'GZMK','GPR183','LGALS3', 
  'ZNF683','CXCR3','',
  'LTB','RTKN2','TCF7')
DotPlot(seurat_obj, features = features[!duplicated(features)])  + coord_flip()+
  scale_colour_distiller(palette = "YlOrRd", direction = 1)+
  theme(axis.text.x = element_text(angle = 0, hjust =0.5),
        axis.text.y = element_text(face = 'italic'),
        axis.title = element_blank())

case = 'named_'
type_colors <- c('#009E73','#F0E442','#D55E00','#56B4E9','#0072B2','#CC79A7')
cell_type <- data.frame(ClusterID = 0:11, cell_type = 'CD8.c01')
cell_type[cell_type$ClusterID %in% c(9), 2] = 'CD8.c02'
cell_type[cell_type$ClusterID %in% c(5), 2] = 'CD8.c03'
cell_type[cell_type$ClusterID %in% c(2), 2] = 'CD8.c04'
cell_type[cell_type$ClusterID %in% c(1,10), 2] = 'CD8.c05'
cell_type[cell_type$ClusterID %in% c(6), 2] = 'CD8.c06'

seurat_obj@meta.data$cell_type <- "NA"
for (i in 1:nrow(cell_type)) {
  seurat_obj@meta.data[which(seurat_obj@meta.data$seurat_clusters == cell_type$ClusterID[i]), 'cell_type'] <-
    cell_type$cell_type[i]
}
new.cluster.ids <- c('CCR7+Tn','RTKN2+Tm','GZMK+Tm','GZMB-Tem','GZMB+Tem','ZNF683+Tex')
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj$type2 <- Idents(seurat_obj)

seurat_obj$type3 <- paste0(seurat_obj$type1, '(', seurat_obj$type2, ')')
table(seurat_obj$type3)
seurat_obj$type3 <- factor(seurat_obj$type3, levels = c(
  'CD8.c01(CCR7+Tn)', 'CD8.c02(RTKN2+Tm)', 'CD8.c03(GZMK+Tm)', 
  'CD8.c04(GZMB-Tem)', 'CD8.c05(GZMB+Tem)', 'CD8.c06(ZNF683+Tex)'
))
seurat_obj$cell_type <- seurat_obj$type3
Idents(seurat_obj) <- 'cell_type'

DimPlot(seurat_obj, pt.size = 3,raster = TRUE, raster.dpi = c(1024, 1024))+
  theme_void() +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"),
                               ends = 'last', type = "closed")) + #修改为箭头坐标系
  guides(color = guide_legend(override.aes = list(size = 5))) +#修改图例散点大小
  labs( x= "UMAP 1",y= "UMAP 2")+
  theme(plot.margin = margin(5.5,15,5.5,5.5),
        panel.grid = element_blank(),aspect.ratio = 1) +
  scale_color_manual(values = type_colors) +
  scale_fill_manual(values = type_colors)
ggsave(file.path(outdir2,paste0(case,'UMAP_seurat_types_lowR.pdf')),
       ggplot2::last_plot(), height=3, width=5)

qs::qsave(seurat_obj, file.path(outdir, paste0(case, "seurat.qs")))

#### 3.Fig4D ----
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
top.markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

case <- 'Markers_'
subobj <- subset(seurat_obj, downsample = 500)
DefaultAssay(subobj) <- 'integrated'

pdf(file.path(outdir2, paste0(case, "top5_heatmap.pdf")), width=6, height=4.5)
DoHeatmap(subobj, features = top.markers$gene, 
                angle = 0, group.colors = type_colors,label = F,
                disp.min=-2.5, disp.max=2.5)+theme(axis.text.y = element_text(face = "italic"))+
  labs(fill = c('z score'))+
  theme(axis.text.y = element_text(colour = 'black'))+
  scale_fill_gradientn(colors = c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
                                  colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)))
dev.off()

#### 4.Fig4E ----
group_colors <- c( "#228B22", "#ffbc14", "firebrick3")
case = 'CD8_'

seurat_obj <- qs::qread('./files/named_seurat.qs')
prop <- as.data.frame(prop.table(table(Idents(seurat_obj), seurat_obj$orig.ident), margin = 2) * 100)
colnames(prop) <- c('types','donor','prop')
prop$group <- ''
prop$group[prop$donor%in%c("BLBL804" ,"BLJK804" ,"CNLHL420")] <- 'Normal'
prop$group[prop$donor%in%c("CFYY331" ,"CGCJ1028" ,"CLJJ329")] <- 'PE'
prop$group[prop$donor%in%c("BMYY203","BFYY331"  ,"BLJJ329","BHLJ1022","BJZY1103")] <- 'LM'
prop$donor <- factor(prop$donor, c("BLBL804" ,"BLJK804" ,"CNLHL420",
                                   "BMYY203","BFYY331"  ,"BLJJ329","BHLJ1022","BJZY1103",
                                   "CFYY331" ,"CGCJ1028" ,"CLJJ329"))
prop$group <- factor(prop$group,levels = c('Normal','LM','PE'))
prop$types <- factor(prop$types,levels = celltypes)

p_box <- ggboxplot(prop, x = 'types', y = 'prop',
                   color ='group',repel = F,
                   palette = group_colors,add = "jitter", 
                   width = 0.5, size = 0.5, legend = 'right') +
  labs(x = ' ', y = 'Percentage of cells (%)')+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=1),
        axis.ticks = element_blank())

library(rstatix)
stat_t <- t_test(group_by(prop, types), prop~group)
stat_t <- add_significance(stat_t, 'p')
stat_t.test <-  add_xy_position(stat_t, x = 'types', dodge = 0.8)
p_box + stat_pvalue_manual(stat_t.test, label = 'p.signif', tip.length = 0.01,
                           bracket.nudge.y =-2,
                           color = 'black',label.size = 5,label.y = c(29, 35, 40),
                           hide.ns = "p",
                           bracket.size = 0.2)
ggsave(file.path(outdir2, paste0(case, "box_group_type.pdf")),last_plot(),width =5.5,height = 5)

#### 5.Fig4F ----
seurat_obj <- qs::qread('./files/named_seurat.qs')
DefaultAssay(seurat_obj) <- 'RNA'
case <- 'GZMB_'
markers <- FindMarkers(seurat_obj,
                       ident.1 = "CD8.c05(GZMB+Tem)",
                       ident.2 = "CD8.c04(GZMB-Tem)",
                       min.pct = 0,
                       logfc.threshold = 0)
markers$gene <- rownames(markers)
logFC.cutoff <-0.25
pvalue.cutoff <- 0.05
markers$change <-
  as.factor(ifelse(
    markers$p_val < pvalue.cutoff &
      abs(markers$avg_log2FC) > logFC.cutoff,
    ifelse(markers$avg_log2FC > logFC.cutoff , 'UP', 'DOWN'),
    'NOT'
  ))
saveRDS(markers, file = file.path(outdir,paste0(case, "deg.rds")))

##
library(ggrepel)
library(ggrastr)

adata <- readRDS('./files/GZMB_deg.rds')
genes <- c('GZMB','CX3CR1','KLRD1','LGALS1','TBX21','CCL4','TGFB1','KLF2','RUNX3','TCF7',
           'CD28','MYC','GZMK','IL7R','CD27','LTB','CXCR3', 'CXCR4','JUNB','TIGIT')
adata %>% glimpse()
adata$label <- ""
adata$label[adata$gene %in% genes] <- adata[adata$gene %in% genes, ]$gene
p <- catvolcano(adata, x = avg_log2FC, y = p_val,
                log2FC = 0.25, p_value = 0.05, text = genes)
ggsave(file.path(outdir2, paste0(case, "DEG_volcano.pdf")), p, height = 4, width = 5, dpi = 600)

#### 6.Fig4G ----
seurat_obj <- qs::qread('./files/named_seurat.qs')
markers <- FindAllMarkers(seurat_obj, assay = "integrated", only.pos = T)
saveRDS(markers, "./files/markers.rds")

data <- as(as.matrix(seurat_obj@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat_obj@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
data <- data[,rownames(pd)]
cds <- newCellDataSet(
  data,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10)) 
saveRDS(cds, "./files/cds_ori.rds")

top_markers <- markers %>% 
  filter(p_val <= 0.05) %>% 
  arrange("avg_log2FC") %>% 
  group_by(cluster) %>% 
  top_n(10, wt = avg_log2FC) %>% 
  pull(gene)
ordergene <- top_markers
cds <- setOrderingFilter(cds, ordergene)
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree', verbose = F)
cds <- orderCells(cds)
saveRDS(cds, file.path(outdir, "cds.rds"))

p <- plot_cell_trajectory(cds, color_by = "Pseudotime",
                           cell_size = 0.1, size = 1, show_backbone = TRUE)
ggsave(file.path(outdir2, 'Pseudotime.pdf'), p, width = 3.5, height = 4)

p <- plot_cell_trajectory(cds, color_by = "cell_type",
                           cell_size = 0.1, size = 1, show_branch_points = F, show_backbone = TRUE) +
  scale_color_manual(values = type_colors)
ggsave(file.path(outdir2, 'cell_type.pdf'), p, width = 3.5, height = 4)

#### 7.Fig4H ----
all_deg_CD8 <- readRDS('./DEG/DEG_LM_N.rds')
all_deg_CD8 <- all_deg_CD8 |> filter( p_val_adj < 0.001 & abs(avg_log2FC) >0.3)
ordergene2 <- all_deg_CD8[all_deg_CD8$change %in% c('DOWN', 'UP'), ]$gene
TFs <- c("TBX21", "EOMES", "PRDM1", "STAT4", "GATA3", "FOXO1", "TCF7", "MYB",  "RUNX3",
         "ZBTB7A", 'JUND', 'CEBPB', 'PITX1')
ordergene2 <- c(ordergene2, TFs)
ordergene2 <- ordergene2[!duplicated(ordergene2)]
Time_diff <- differentialGeneTest(cds[ordergene2,],
                                  cores = 1,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
rownames(Time_diff) <- Time_diff$gene_short_name
pdf(file=file.path(outdir2, "non-branched_heatmap_label.pdf"), width= 5, height= 5)
plot_pseudotime_heatmap(
  cds[row.names(subset(Time_diff, qval < 0.01)),],
  num_clusters = 5, cores = 1,
  use_gene_short_name = T,
  show_rownames = T)
dev.off()

##
plotdf=pData(cds)
ggplot(plotdf, aes(x=Pseudotime,y=cell_type,fill=cell_type))+
  geom_density_ridges(scale=1) +
  scale_y_discrete("")+
  scale_fill_manual(values = type_colors) +
  theme_minimal()+  
  theme(axis.title = element_blank(),
        axis.text = element_text(face = 'bold'),
        panel.grid = element_blank())+
  scale_x_continuous(position = "top")+
  scale_y_discrete(position = "right")+NoLegend()
ggsave(file.path(outdir2, "celltype_tmp1.pdf"),
       width = 20,height = 5,units = "cm")

#### 8.Fig4I ----
s_genes <- c('GZMB','PRDM1', "TBX21", "STAT4",  "RUNX3", "ZBTB7A", 'JUND')
p <- plot_genes_in_pseudotime(cds[s_genes,], color_by = "cell_type")+
  scale_color_manual(values = type_colors)
ggsave(file.path(outdir2, 'genes_s3_pseudotime.pdf'), p, width = 6, height = 6)

#### 9.Fig4J ----
seurat_obj <- qs::qread('./files/named_seurat.qs')
seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                              species = "human",
                              outdir = "./TF/results/")

TF_G <- read.table('./TF/results/tfs_target.tsv',sep = ',',header = T) 
TF_G <- TF_G[!duplicated(TF_G[c('tf')]), ]
regulonAUC <- importAUCfromText("auc_mtx.csv")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
cellInfo <- seurat_obj@meta.data
seurat_obj[["scenic"]] <- CreateAssayObject(counts = t(getAUC(regulonAUC)))

genes <- c("TBX21", "EOMES", "PRDM1", "STAT4", "GATA3", "FOXO1", "TCF7", "MYB",  "RUNX3",
           "ZBTB7A", 'JUND', 'CEBPB', 'PITX1')
regulonActivity_byCell <- getAUC(regulonAUC)
seurat_obj$cell_id <- colnames(seurat_obj)
annotation_col <- as.data.frame(seurat_obj@meta.data[, c("orig.ident",'cell_type','group')])
rownames(annotation_col) <- rownames(seurat_obj@meta.data)

annotation_colors <- list(
  group=c('Normal'='#228B22',
          'LM'='#ffbc14',
          'PE'="firebrick3"),
  cell_type = c('CD8.c01(CCR7+Tn)'='#009E73',
                'CD8.c02(RTKN2+Tm)'='#F0E442',
                "CD8.c03(GZMK+Tm)"='#D55E00',  
                "CD8.c04(GZMB-Tem)"='#56B4E9',
                "CD8.c05(GZMB+Tem)"='#0072B2',
                "CD8.c06(ZNF683+Tex)"='#CC79A7'),
  orig.ident = c("BLBL804" = "#53A85F",
                 "BLJK804" = "#23452F",
                 "CNLHL420" = "#91D0BE",
                 'BFYY331'="#F1BB72",
                 'BLJJ329'="#D6E7A3",
                 'BHLJ1022'="#c5942e",
                 'BJZY1103'="#fabf04",
                 'BMYY203'="#FFFF00",
                 'CFYY331'="#F3B1A0",
                 'CLJJ329'="#E95C59",
                 'CGCJ1028'="#E59CC4"))

colnames(regulonActivity_byCell) <- gsub(',','',colnames(regulonActivity_byCell))
colnames(regulonActivity_byCell) <- TF_G$symbol
scaled <- scale(t(regulonActivity_byCell), center = T, scale = T)
scaled[scaled >= 1.3] <- 1.3
scaled[scaled <= -1.3] <- -1.3
my.breaks <- c(seq(-2, 0, by=0.01), seq(0.1, 2, by=0.01))
my.colors <- c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(length(my.breaks)/2),
               colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(length(my.breaks)/2))
pdf(file.path(outdir2, "TF_s_heatmap2.pdf"), width = 5,height = 2)
pheatmap(scaled[tf.seleted, ],
         annotation_col = annotation_col, annotation_colors = annotation_colors,
         color = my.colors,
         treeheight_row = 10, treeheight_col = 10,
         cellwidth = 0.005, cellheight = 5,
         fontsize = 5, border_color = "white", angle_col = "45",
         cluster_rows = T, cluster_cols = T,
         show_rownames = T, show_colnames = F)
dev.off()

#### 10.Fig4K ----
seurat_obj0 <- qs::qread('./files/named_seurat.qs')
seurat_obj <- subset(seurat_obj0, idents = 'CD8.c04(GZMB-Tem)')
Idents(seurat_obj) <- 'group'

TF_G <- read.table('./TF/results/tfs_target.tsv',sep = ',',header = T) 
TF_G <- TF_G[!duplicated(TF_G[c('tf')]), ]
regulonAUC <- importAUCfromText('./TF/results/auc_mtx.csv')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonAUC <- regulonAUC[colnames(seurat.obj), ]
cellTypes <- seurat.obj@meta.data[,c('group',"orig.ident")]
table(cellTypes$orig.ident)
rss <- calcRSS(AUC=t(getAUC(regulonAUC)), cellAnnotation=cellTypes[rownames(regulonAUC), 'group'])
rss=na.omit(rss) 
rownames(rss) <- gsub(',','',rownames(rss))
rownames(rss) <- TF_G$symbol
rssPlot <- plotRSS(rss,zThreshold = 0.5)
plotly::ggplotly(rssPlot$plot)

data <- as.data.frame(rss)
data$TFs <- rownames(data)
data$gene_name = rownames(data)
data$rss <- data$Normal
data$hits <- if_else(data$LM < data$Normal,
                     "down", "no")
data <- data[order(data$Normal,decreasing = T),]
data$rank <- 1:length(data$gene_name)
data <- data[,c('rss','hits','rank')]

pointed_out <- data |> filter(rownames(data) %in% 
                                c("TBX21 (260g)", "EOMES (371g)", "PRDM1 (294g)",
                                  "FOXO1 (22g)", "MYB (95g)",  "RUNX3 (54g)",
                                  "ZBTB7A (94g)", 'JUND (8g)', 'CEBPB (500g)', 'PITX1 (40g)'))
ggplot(data, aes(x=rank, y=rss, color=hits))+
  geom_point(size=3,color="grey50",alpha =0.4)+
  geom_point(data = pointed_out,  stroke = 0.5, size=4, shape=16, color="#DC050C",alpha =0.6)+
  ggrepel::geom_text_repel(data=pointed_out, aes(label=rownames(pointed_out)), color="black", 
                           size=4, fontface="italic", segment.size=0.5, nudge_x=60, 
                           direction="y", hjust=0,ylim = c(0.06, 0.50))+ theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title = element_text(colour = 'black', size = 15),
        axis.text = element_text(colour = 'black', size = 12))+
  labs(x='Rank', y='RSS')
ggsave(file.path(outdir2,'TF_rank_N.pdf'), last_plot(), width = 4,height = 5)

#### 11.Fig4L ----
tf.seleted <- c("FOXO1 (22g)","STAT4 (12g)","MYB (95g)","PRDM1 (294g)","TBX21 (260g)",
                "EOMES (371g)","ZBTB7A (94g)","PITX1 (40g)","JUND (8g)","RUNX3 (54g)","CEBPB (500g)")
seurat_obj3 <- subset(seurat_obj, idents = 'CD8.c05(GZMB-Tem)')
regulonAUC <- importAUCfromText("./TF/results/auc_mtx.csv")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
regulonAUC <- regulonAUC[colnames(seurat_obj3),]
cellInfo <- seurat_obj3@meta.data
seurat_obj3[["scenic"]] <- CreateAssayObject(counts = t(getAUC(regulonAUC)))
Idents(seurat_obj3) <- 'group'
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), Idents(seurat_obj3)),
                                     function(cells)
                                       rowMeans(t(getAUC(regulonAUC))[, cells]))
rownames(regulonActivity_byCellType) <- gsub(',','',rownames(regulonActivity_byCellType))
rownames(regulonActivity_byCellType) <- TF_G$symbol
pdf(file.path(outdir2, "TF_gene_group_sTF_GZMB.pdf"), width = 2,height = 2)
p <- pheatmap(
  regulonActivity_byCellType[tf.seleted,],
  scale = "row",
  color = c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
            colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)),
  treeheight_row = 10, treeheight_col = 10,
  border_color = "white",
  cellheight = 5, cellwidth = 10,
  fontsize = 5,
  angle_col = "45",
  main = 'CD8.c05(GZMB+Tem)',
  cluster_cols = F, cluster_rows = T)
dev.off()
