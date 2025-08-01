library(Seurat)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(ggplot2)
library(UpSetR)
library(qs) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(ggrastr)
set.seed(777)

setwd('~/LM/Fig2/')
outdir <- "./files/"
dir.create(outdir,recursive = T)
outdir2 <- "./plots/"
dir.create(outdir2,recursive = T)

type_colors <- c('#53A85F', '#E59CC4','#F1BB72', '#E95C59', '#F3B1A0', '#8C549C', 
                 '#f46f20', "#6a73cf",'#0eb0c8', '#476D87')
celltypes <- c('B cells', 'CD4+ T cells', 'CD8+ T cells', 'MAIT', 'Gamma delta T cells', 'NK cells',
               'CD14+ Monocytes', 'CD14+CD16+ Monocytes', 'CD16+ Monocytes', 'Dendritic cells')
names(type_colors) <- celltypes

#### 1.Fig2A ----
DEG <- readRDS("./DEG/DEG_LM_N.rds")|>
  filter(abs(avg_log2FC) > 0.25 & p_val_adj < 0.05)
DEG_split <- split(DEG, DEG$cluster)
for(i in 1:length(table(DEG$cluster))){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}
upsetData=fromList(DEG_split)
case = 'all_'

ct_colors <- c('#6a73cf','#0eb0c8','#f46f20','#476D87','#8C549C','#F3B1A0','#53A85F','#E59CC4','#E95C59','#F1BB72')
pdf(file=file.path(outdir2, paste0(case, 'upset.pdf')),onefile = FALSE,width=6,height=5)
upset(upsetData,
      nsets = length(DEG_split),            
      nintersects = 24,             
      order.by = "freq",       
      show.numbers = "yes",           
      number.angles = 0,             
      point.size = 3,                 
      matrix.color="darkgrey",             
      line.size = 0.8,                     
      mainbar.y.label = "Differentially expressed genes(DEGs)",
      sets.x.label = "Set Size",
      sets.bar.color = ct_colors,
      main.bar.color ="black", 
      queries = list(
        list(query=intersects, params=list(celltypes), color="#BEBADA", active=T),
        list(query=intersects, params=list("CD14+ Monocytes"), color="#BEBADA", active=T),
        list(query=intersects, params=list("CD16+ Monocytes"), color="#BEBADA", active=T),
        list(query=intersects, params=list('CD14+CD16+ Monocytes'), color="#BEBADA", active=T),
        list(query=intersects, params=list('NK cells'), color="#BEBADA", active=T),
        list(query=intersects, params=list('CD14+ Monocytes', 'CD14+CD16+ Monocytes'), color="#BEBADA", active=T),
        list(query=intersects, params=list('CD16+ Monocytes', 'CD14+CD16+ Monocytes'), color="#BEBADA", active=T),
        list(query=intersects, params=list('CD14+ Monocytes', 'CD14+CD16+ Monocytes', 'CD16+ Monocytes'), color="#BEBADA", active=T),
        list(query=intersects, params=list('CD14+ Monocytes', 'CD14+CD16+ Monocytes', 'CD16+ Monocytes', "Dendritic cells"), color="#BEBADA", active=T)
      )) 
dev.off()

intersectGenes=Reduce(intersect,DEG_split)
write.table(file=file.path(outdir, paste0(case, 'intergenes.txt')),
            intersectGenes,sep=",",quote=F,col.names=T,row.names=F)

#### 2.Fig2B ----
marker_condition <- readRDS("./DEG/DEG_LM_N.rds")
logFC.cutoff <- 0.25
pvalue.cutoff <- 0.05
marker_condition$logP <- -log10(marker_condition$p_val_adj)
marker_condition$type <- marker_condition$cluster |> as.character()
marker_condition$type[marker_condition$change == 'NOT'] <- 'NOT'
marker_condition$type <- factor(marker_condition$type, 
                                levels = c('B cells','CD4+ T cells','CD8+ T cells','MAIT','Gamma delta T cells',
                                           'NK cells','CD14+ Monocytes','CD14+CD16+ Monocytes','CD16+ Monocytes',
                                           'Dendritic cells','NOT'))
type_colors <- c('#53A85F','#E59CC4','#F1BB72','#E95C59','#F3B1A0','#8C549C','#f46f20',"#6a73cf",'#0eb0c8','#476D87','grey')
genes <- c('FCGR1A','CD163','S100A8','HLA-DRB1')
marker_condition <- marker_condition %>% mutate(label = if_else(gene %in% genes & change != 'NOT', gene, ""))
Volcano_paired <-
  rasterise(ggscatter(marker_condition,
                      x = "avg_log2FC", y = "logP", color = "type",
                      palette = type_colors,
                      size = 0.5, font.label = 8, repel = T,
                      xlab = "log2 Foldchange", ylab = "-log10 (adjust pvalue)"), dpi = 600) +
  ggrepel::geom_text_repel(aes(label = label), max.overlaps = Inf,
                           color = "black", size = 3, fontface = "italic") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-logFC.cutoff, logFC.cutoff), linetype = "dashed") +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = "right",
        plot.title = element_text(size = 12, color = "black", hjust = 0.5),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank()) +
  guides(color=guide_legend(title = "cell type",override.aes = list(size=4,alpha=1)))+
  scale_x_continuous(limits = c(-max(abs(marker_condition$avg_log2FC)), max(abs(marker_condition$avg_log2FC))))
ggsave(paste0(outdir2,"types_up_volcano.pdf"),
       Volcano_paired, height = 6, width = 6, dpi = 600)

#### 3.Fig2C-D ----
library(enrichplot)
type_colors <- c('#53A85F', '#f46f20',"#6a73cf",'#0eb0c8','#E59CC4','#F1BB72','#476D87','#F3B1A0','#E95C59','#8C549C')

## UP
DEG <- readRDS("./DEG/DEG_LM_N.rds") |> filter(avg_log2FC >= 0.25 & p_val_adj < 0.05) 
DEG_split <- split(DEG, DEG$cluster)
for(i in 1:length(table(DEG$cluster))){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}
data <- DEG_split
data <- lapply(X = data, FUN = function(x) {
  x <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  x <- x[,-1]
})
GO_pathway <- compareCluster(data, 
                             fun='enrichGO', ont= 'BP',
                             OrgDb='org.Hs.eg.db' ,
                             pvalueCutoff = 0.05, qvalueCutoff = 1)
GO_results <- pairwise_termsim(GO_pathway)
saveRDS(GO_results, file.path(outdir, 'GO_results_UP.rds'))

df <- GO_results@compareClusterResult
terms <- c(
  'defense response to virus',
  'response to type I interferon',
  'cytokine-mediated signaling pathway',
  'response to lipopolysaccharide',
  'immune response-regulating signaling pathway',
  'interferon-mediated signaling pathway',
  'neutrophil chemotaxis',
  'granulocyte chemotaxis',
  'leukocyte mediated cytotoxicity',
  'positive regulation of NF-kappaB transcription factor activity',
  'defense response to bacterium',
  'positive regulation of leukocyte activation',
  'type I interferon-mediated signaling pathway',
  'response to tumor necrosis factor'
)
terms[!(terms %in% GO_results@compareClusterResult$Description)]
df2 <- df |> filter(Description %in% terms)
p1 <- df2 |> ggplot(aes(x = Description, y = Cluster, 
                        color = -log10(pvalue), size = Count)) +
  geom_point(shape = 15) +
  theme_bw()+
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        panel.spacing.x = unit(0, "pt")) +
  coord_fixed() +
  guides(x = guide_axis(angle = 30)) +
  scale_color_gradient2(low = "#0571b0", mid = "white", high = "#ca0020")
ggsave(file.path(outdir2, 'Pathway_heatmap_UP.pdf'),
       p1, height=5, width=6)

## DOWN
DEG <- readRDS("./DEG/DEG_LM_N.rds") |> filter(avg_log2FC <= -0.25 & p_val_adj < 0.05) 
DEG_split <- split(DEG, DEG$cluster)
for(i in 1:length(table(DEG$cluster))){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}
data <- DEG_split
data <- lapply(X = data, FUN = function(x) {
  x <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  x <- x[,-1]
})
GO_pathway <- compareCluster(data, 
                             fun='enrichGO', ont= 'BP',
                             OrgDb='org.Hs.eg.db' ,
                             pvalueCutoff = 0.05, qvalueCutoff = 1)
GO_results <- pairwise_termsim(GO_pathway)
saveRDS(GO_results, file.path(outdir, 'GO_results_DOWN.rds'))

terms <- c(
  'MHC class II protein complex assembly',
  'peptide antigen assembly with MHC protein complex',
  'positive regulation of leukocyte cell-cell adhesion',
  'positive regulation of lymphocyte activation',
  'positive regulation of leukocyte activation',
  'lymphocyte mediated immunity',
  'lymphocyte proliferation',
  'cell killing',
  'immunoglobulin production',
  'cellular defense response',
  'regulation of cell killing',
  'leukocyte mediated cytotoxicity',
  'regulation of cellular senescence',
  'positive regulation of cellular senescence'
)
terms[!(terms %in% GO_results@compareClusterResult$Description)]
df2 <- df |> filter(Description %in% terms)
## 绘图
p2 <- df2 |> ggplot(aes(x = Description, y = Cluster,
                        color = -log10(pvalue), size = Count)) +
  geom_point(shape = 15) +
  theme_bw()+
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        panel.spacing.x = unit(0, "pt")) +
  coord_fixed() +
  guides(x = guide_axis(angle = 30)) +
  scale_color_gradient2(low = "#ca0020", mid = "white", high = "#2166AC")
ggsave(file.path(outdir2, 'Pathway_heatmap_DOWN.pdf'),
       p2, height=5, width=6)

#### 4.Fig2E ----
library(Seurat)
library(ggplot2)
library(GSVA)
library(GSEABase)
library(msigdbr)
library(dplyr)
library(stringr)
library(AUCell)
library(clusterProfiler)
library(ggraph)
library(ggpubr)
library(ggrain)

group_colors <- c( "#228B22", "#ffbc14")
seurat_obj <- qs:qread('~/LM/merge/files/named_seurat.qs')
seurat_obj <- subset(seurat_obj, idents=c('B cells','CD4+ T cells','CD8+ T cells','MAIT','Gamma delta T cells','NK cells',
                                          'CD14+ Monocytes','CD14+CD16+ Monocytes','CD16+ Monocytes','Dendritic cells'))
qs::qsave(seurat_obj, '~/LM/merge/files/impor_seurat.qs')

seurat_obj <- subset(seurat_obj, group%in%c('LM','Normal'))
seurat_obj$group <- factor(seurat_obj$group, levels = c('Normal','LM'))
DefaultAssay(seurat_obj) <- 'RNA'

cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@data) 
category <- "H"
Dataset <- msigdbr(species = "Homo sapiens", category = category)

feature_list <- c('HALLMARK_TGF_BETA_SIGNALING', 
                  'HALLMARK_INFLAMMATORY_RESPONSE')
dataset <- Dataset %>% filter(Dataset$gs_name %in% feature_list)
geneSets <- lapply(unique(dataset$gs_name), 
                   function(x){dataset$gene_symbol[dataset$gs_name == x]})
names(geneSets) <- unique(dataset$gs_name)

names(geneSets) <- names(geneSets)|>
  str_replace_all("HALLMARK_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
s_term <- c('Tgf beta signaling', 'Inflammatory response')
for(i in s_term){
  cells_AUC <- AUCell_calcAUC(geneSets[[i]], cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
  aucs <- as.numeric(getAUC(cells_AUC))
  seurat_obj$AUC  <- aucs
  
  meta.data <- seurat_obj@meta.data[, c("group", "AUC")]
  
  p <- ggplot(meta.data, aes(group, AUC, fill = group)) +
    geom_rain(alpha = .4, 
              boxplot.args.pos = list( width = 0.05, position = position_nudge(x = 0.13)),
              violin.args.pos  = list( width = 0.70, position = position_nudge(x = 0.2),side = "r")) + 
    theme_classic() +
    theme(
      aspect.ratio = 0.6,
      axis.text = element_text(size = 6,color="black"), 
      axis.title = element_text(size = 8),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
      plot.title = element_text(size = 8, vjust = 1, hjust = 0),
      legend.text = element_text(size = 6))+
    scale_color_manual(values = group_colors)+
    scale_fill_manual(values = group_colors)+
    labs(title =  i |> str_replace_all("HALLMARK_", " ") |>
           str_replace_all("_", " ")|>
           str_wrap(width = 40)|>
           str_to_sentence(), x= '', y='AUC score') +
    guides(fill = 'none', color = 'none')+ 
    ggsignif::geom_signif(comparisons = list(c("LM","Normal")),
                          test ="wilcox.test",
                          map_signif_level =c("***"=0.001, "**"=0.01, "*"=0.05),
                          size =0.5, 
                          textsize =3)+coord_flip()
  ggsave(file.path(outdir2,paste0(i,'_score_rain.pdf')), p, height = 4,width = 4)
}

#### 5.Fig2F ----
av_expr <-
  AverageExpression(
    seurat_obj,
    assays = "RNA",
    group.by = "group"
  )[["RNA"]] %>%
  as.data.frame()
cytokine[!(cytokine %in% rownames(av_expr))]
av_expr2 <- av_expr[cytokine, ]

cytokine <- c(
  'IL6', 'IFNG','TNFRSF1A','CSF2RA', 'CSF2RB', 
  'IRF7', 'IL1R2', 'ACSL1', 'ISG15', 'IFITM2', 'IFITM1',
  'CSF2', 'TNF', 'IL1B','LTA')
av_expr2 <- av_expr[cytokine, ]
p1 <- pheatmap::pheatmap(
  av_expr2,
  scale = "row",
  border_color = "white",
  cellwidth = 10,
  cellheight = 10,
  fontsize = 10,
  cluster_rows = F,
  cluster_cols = F,
  number_color = "white",
  color = c(colorRampPalette(colors = c("#FDDBC7","#F4A582","#B2182B"))(10)))
ggsave(file.path(outdir2, "cytokine_heatmap.pdf"), p1, height=3, width=2)