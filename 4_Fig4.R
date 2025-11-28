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

DefaultAssay(seurat_obj) <- 'RNA'
Idents(seurat_obj) <- 'cell_type'
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
ggsave(file.path(outdir2,paste0(case,'markers_types.pdf')),
       ggplot2::last_plot(), height=5, width=5)

#### 3.Fig4D ----
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


