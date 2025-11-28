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

setwd('~/LM/Fig6/')
outdir <- "./files/"
dir.create(outdir,recursive = T)
outdir2 <- "./plots/"
dir.create(outdir2,recursive = T)

type_colors <- c("#2b9b81","#92c051")
group_colors <- c( "#228B22", "#ffbc14", "firebrick3")

#### 1.Fig6A ----
seurat_obj <- qs::qread('./files/named_seurat.qs')
Idents(seurat_obj) <- factor(Idents(seurat_obj), levels = c('NK CD56highCD16low', 'NK CD56lowCD16high'))
DimPlot(seurat_obj, pt.size = 5, raster = TRUE, raster.dpi = c(1024,1024))+
  theme_void() +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), ends = 'last', type = "closed")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs( x= "UMAP 1",y= "UMAP 2")+
  theme(plot.margin = margin(5.5,15,5.5,5.5),
        panel.grid = element_blank(),aspect.ratio = 1) +
  scale_color_manual(values = type_colors) +
  scale_fill_manual(values = type_colors)
ggsave(file.path(outdir2,paste0(case,'UMAP_celltype.pdf')), last_plot(), height=2, 

#### 2.Fig6B ----
features <- c('NKG7', 'PRF1','FCGR3A','FGFBP2','CCL5','GZMB','IL2RG',
              'NCAM1', 'XCL1','NCR1','IL7R','CCR7','IL2RB','IL2RA')
DotPlot(seurat_obj,
        features= features,)+
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,face = 'italic'),
        axis.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.line = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.ticks = element_blank())+
  scale_color_distiller(palette= "RdBu")
ggsave(file.path(outdir2,'marker_Dotplot2.pdf'), last_plot(), height=3, width=6)

#### 3.Fig6C ----
case <- 'NK2_'
markers <- FindMarkers(seurat_obj,
                       ident.1 = "NK CD56bright",
                       ident.2 = "NK CD56dim",
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
deg <- readRDS(file.path(outdir,paste0(case, "deg.rds")))
deg_up <- deg[deg$change=="UP",]
deg_down <- deg[deg$change=="DOWN",]

bp1 <-
  enrichGO(
    deg_up$gene,
    OrgDb =org.Hs.eg.db, 
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)
term1 <- bp1@result
up_pathway <- c('leukocyte migration',
                'lymphocyte mediated immunity',
                'positive regulation of cytokine production',
                'cell chemotaxis',
                'type II interferon production',
                'regulation of lymphocyte mediated immunity',
                'cellular response to interleukin-4')
df1 <- term1[term1$Description %in% up_pathway, ]
df1$labelx=rep(0,nrow(df1))
df1$labely=seq(nrow(df1),1)

bp2 <-
  enrichGO(
    deg_down$gene,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)
term2 <- bp2@result
down_pathway <- c('I-kappaB kinase/NF-kappaB signaling',
                  'cell killing',
                  'cytolysis',
                  'leukocyte mediated cytotoxicity',
                  'natural killer cell chemotaxis',
                  'cellular response to tumor necrosis factor',
                  'natural killer cell mediated cytotoxicity')
df2 <- term2[term2$Description %in% down_pathway, ]
df2$labelx=rep(0,nrow(df2))
df2$labely=seq(nrow(df2),1)

data <- rbind(df1[df1$Description%in%up_pathway,],df2[df2$Description%in%down_pathway,])
a <- as.numeric(length(up_pathway))
b <- as.numeric(length(down_pathway))
data$change <- c(rep('NK CD56bright',a),rep('NK CD56dim',b))
data$group <- c(rep(1,a),rep(-1,b))
data$pvalue <- as.numeric(data$pvalue)
data$LogP <- -log10(data$pvalue)
data$LogP[8:14] <-  -data$LogP[8:14]

ggplot(data,aes(x=reorder(Description, LogP),y= LogP,fill=change))+
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.2)) +
  theme_bw()+coord_flip()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        panel.border = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.line.x = element_line(color='black'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none')+
  geom_text(data = data[which(data$group>0),],aes(x=Description, y=-0.1, label=Description),
            hjust=1, size=4)+
  geom_text(data = data[which(data$group<0),],aes(x=Description, y=0.1, label=Description),
            hjust=0, size=4)+
  scale_fill_manual(values = type_colors)+ 
  scale_x_discrete(labels=function(x) str_wrap(x, width=20), expand = expansion(mult = c(0,0)))+
  scale_y_continuous(breaks = seq(-9, 15, 3))+
  labs(x='', y='-log10(pvalue)')+NoLegend()
ggsave(file.path(outdir2, paste0(case, 'UP_DOWN_Barplot.pdf')), last_plot(), height=5, width=7)

#### 4.Fig6D ----
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
prop2 <- prop[prop$types %in% c("NK CD56highCD16low","NK CD56lowCD16high"),]
stat_t <- t_test(group_by(prop2, types), prop~group)
stat_t <- add_significance(stat_t, 'p')

stat_t.test <-  add_xy_position(stat_t, dodge = 3)
prop2$group <- factor(prop2$group,levels = c('Normal','LM','PE'))
prop2$types <- factor(prop2$types,levels = c("NK CD56highCD16low","NK CD56lowCD16high"))
p_box <- ggboxplot(prop2, x = 'group', y = 'prop',
                   facet.by = 'types', fill = 'group',
                   color = 'black',repel = F,
                   width = 0.5, size = 0.3, legend = 'right') +
  scale_fill_manual(values = group_colors) +
  labs(x = '', y = 'Percentage of cells (%)', fill = 'types')+
  theme(axis.text.x = element_text(angle=0,hjust=0.5,vjust=1),
        axis.text.x.bottom = element_text(size = 10),
        aspect.ratio = 1, axis.ticks = element_blank())+NoLegend()+ 
  facet_wrap(~ types, ncol = 2)
ggsave(file.path(outdir2, "box_face_group_type_group.pdf"), last_plot(), width = 5, height = 3)

#### 5.Fig6E ----
cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@data)
category <- "C5"
Dataset <- msigdbr(species = "Homo sapiens", category = category)
Dataset <- Dataset %>% filter(str_detect(Dataset$gs_name,"GOBP_"))
Dataset$gs_name2 <- Dataset$gs_name|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ")|>
  str_to_sentence()
geneSets <- lapply(unique(Dataset$gs_name),
                   function(x){Dataset$gene_symbol[Dataset$gs_name == x]})
names(geneSets) <- unique(Dataset$gs_name)
names(geneSets) <- unique(Dataset$gs_name2)
s_term <- c(
  'Antibody dependent cellular cytotoxicity',
  'Natural killer cell activation',
  'Natural killer cell activation involved in immune response',
  'Natural killer cell chemotaxis',
  'Natural killer cell cytokine production',
  'Natural killer cell degranulation',
  'Natural killer cell differentiation',
  'Natural killer cell mediated immunity',
  'Natural killer cell proliferation',
  'Regulation of natural killer cell activation',
  'Positive regulation of natural killer cell activation',
  'Positive regulation of natural killer cell chemotaxis',
  'Positive regulation of natural killer cell mediated cytotoxicity',
  'Positive regulation of natural killer cell mediated immunity',
  'Positive regulation of natural killer cell proliferation',
  'Positive regulation of nk t cell activation',
  'Protection from natural killer cell mediated cytotoxicity',
  'Regulation of natural killer cell chemotaxis',
  'Regulation of natural killer cell differentiation',
  'Regulation of natural killer cell mediated immunity',
  'Regulation of nk t cell activation',
  'Regulation of nk t cell proliferation',
  'Positive regulation of natural killer cell mediated cytotoxicity',
  'Natural killer cell mediated immunity',
  'Cell killing',
  'Regulation of cell killing',
  'Leukocyte cell cell adhesion',
  'Immune response regulating cell surface receptor signaling pathway',
  'Cell surface receptor signaling pathway involved in cell cell signaling',
  'Regulation of inflammatory response',
  'Response to interferon gamma',
  'Response to type i interferon',
  'Leukocyte mediated immunity',
  'Platelet activation',
  'I kappab kinase nf kappab signaling',
  'Positive regulation of cytokine production',
  'Response to tumor necrosis factor',
  'Tumor necrosis factor superfamily cytokine production',
  'Alpha beta t cell activation',
  'Regulation of t cell activation',
  'Leukocyte proliferation',
  'Regulation of lymphocyte apoptotic process',
  'Cellular response to interferon gamma',
  'Interferon gamma mediated signaling pathway',
  'Interferon gamma production',
  'Positive regulation of interferon gamma production',
  'Positive regulation of response to interferon gamma',
  'Regulation of response to interferon gamma',
  'Response to interferon gamma'
)
dataset <- Dataset %>% filter(Dataset$gs_name2 %in% s_term)
length(unique(dataset$gs_name)) #24
geneSets <- lapply(unique(dataset$gs_name), 
                   function(x){dataset$gene_symbol[dataset$gs_name == x]})
names(geneSets) <- unique(dataset$gs_name2)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
saveRDS(AUCell_socre,file.path(outdir, paste0(category,'_',"AUCscore.rds")))

AUCell_socre <- readRDS('./files/C5_AUCscore.rds')
seurat_obj$ID <- rownames(seurat_obj@meta.data)
meta <- AUCell_socre |> as.data.frame()
meta$cell_type <- seurat_obj$group
object <-
  aggregate(meta[, -length(colnames(meta))], list(meta$cell_type), FUN = mean)
rownames(object) <- object$Group.1
object <-object[,-1]
colnames(object) <- colnames(object)|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()

## Vlnplot
s_term <- c('Natural killer cell activation involved in immune response','Natural killer cell cytokine production')
data0 <- seurat_obj[, seurat_obj[["cell_type"]] == 'NK CD56lowCD16high']
data0$ID <- rownames(data2@meta.data)

for(i in s_term){
  data <- data0
  data$group <- factor(data$group, levels = c('Normal', 'LM', 'PE'))
  data$group <- factor(data$group)
  s_geneSets <- geneSets[[i]]
  cells_rankings <- AUCell_buildRankings(data@assays$RNA@data)
  cells_AUC <- AUCell_calcAUC(geneSets[[i]], cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
  aucs <- as.numeric(getAUC(cells_AUC))
  data$AUC  <- aucs
  
  meta.data <- data@meta.data[, c("group", "AUC")]
  meta.data <- meta.data %>% reshape2::melt(id.vars = c("group"))
  colnames(meta.data)[2:3] <- c("signature", "score")
  mean.score <- mean(meta.data[meta.data$group %in% 'Normal',]$score)
  median.score <- meta.data %>% group_by(group) %>%
    summarise(median_score = median(score))
  median.score <- median.score[order(median.score$median_score,decreasing = T),]
  median.score$group <- factor(median.score$group, levels = median.score$group)
  meta.data$group <- factor(meta.data$group)
  compar <- list(c('Normal', 'LM'), c('LM', 'PE'))
  
  p2 <- meta.data %>%
    ggplot(aes(x = group, y = score, fill = '', color =group)) +
    geom_violin(fill = 'NA') +
    geom_boxplot(data = meta.data, aes(x = group, y = score), width = 0.1, alpha = 0.1)+
    geom_hline(yintercept = mean.score, color = "black", linetype = "dashed") +
    scale_y_continuous(i |> str_replace_all("GOBP_", " ") |>
                         str_replace_all("_", " ")|>
                         str_wrap(width = 35)|>
                         str_to_sentence()) +
    theme_test() + scale_fill_manual(values = group_colors) +
    labs(title= names(table(data2$cell_type)))+
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 10),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 8,color = "black"),
          axis.text.x = element_text(angle = 0, vjust = 1,hjust = 0.5)) +
    stat_compare_means(comparisons = compar, label = "p.signif",
                       hide.ns = T, ref.group = ".all.",vjust = 2)+ 
    scale_color_manual(values = group_colors)
  ggsave(file.path(outdir2,paste0(i,'_violin_celltype.pdf')), p2, height = 3,width = 4)
}

#### 6.Fig6F ----
seurat.obj <- qs::qread('./files/named_seurat.qs')
seurat.obj <- hp_run_pyscenic(x = seurat.obj,
                              species = "human",
                              outdir = "./TF/results/")

#### 7.Fig6G ----
##
tf <- 'IRF1(+)'
DefaultAssay(seurat.obj) <- "scenic"
p1 <- FeaturePlot(
  seurat.obj,split.by = 'group',
  features = tf,
  reduction = 'umap',
  min.cutoff = 0,
  cols = c("lightgrey" , "#DE1F1F"),
  pt.size = 1,
  order = T)

DefaultAssay(seurat.obj) <- "RNA"
p2 <- FeaturePlot(
  seurat.obj,
  features = tf,
  reduction = "umap",
  pt.size = 1,
  order = T,
  min.cutoff = 0)+theme(plot.title = element_text(face = "bold.italic"))

p <- p1 + p2 + plot_layout(ncol = 2)
ggsave(file.path(outdir2, 'IRF1_TF.pdf'), p, width = 8, height = 3.5)

## Vlnplot_gene
meta.data <- seurat.obj@meta.data[,c('group','cell_type')] 
meta.data$G <- as.numeric(seurat.obj@assays$RNA['IRF1',])
colnames(meta.data)[2:3] <- c("signature", "score")
compar <- list(c('LM', 'Normal'), c('PE', 'LM'))
mean.score <- mean(meta.data$score)
median.score <- meta.data %>% group_by(group) %>% summarise(median_score = median(score))

p <- meta.data %>%
  ggplot(aes(x = group, y = score, fill = group)) +
  geom_violin() +
  geom_point(data = median.score,color = "black", aes(x = group, y = median_score), shape = 3, size = 3) +
  scale_y_continuous("Normalized expression of IRF1")+
  theme_test() + scale_fill_manual(values = type_colors) +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 10,color = "black" ),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8,color = "black"),
        axis.text.x = element_text(angle = 0)) +
  stat_compare_means(comparisons = compar,label.y = 4, label = "p.signif", vjust = 2 )+
  scale_fill_manual(values = group_colors)
ggsave(file.path(outdir2,paste0('G_IRF1_group_violin.pdf')), p, height = 3,width = 4)

## Vlnplot_TF
seurat.obj@meta.data$G <- as.numeric(seurat.obj@assays$scenic['IRF1(+),',])
meta.data <- seurat.obj@meta.data[, c("group", "G")]
meta.data <- meta.data %>% reshape2::melt(id.vars = c("group"))
colnames(meta.data)[2:3] <- c("signature", "score")
mean.score <- mean(meta.data$score)
median.score <- meta.data %>% group_by(group) %>% summarise(median_score = median(score))
median.score <- median.score[order(median.score$median_score,decreasing = T),]
median.score$group <- factor(median.score$group, levels = median.score$group)
meta.data$group <- factor(meta.data$group)
compar <- list(c('LM', 'Normal'), c('PE', 'LM'))

p <- meta.data %>%
  ggplot(aes(x = group, y = score, fill = '', color =group)) +
  geom_violin(fill = 'NA') +
  geom_boxplot(data = meta.data, aes(x = group, y = score), width = 0.1, alpha = 0.1)+
  scale_y_continuous(paste0('Transcription factor activity of IRF1(+)'))+
  theme_test() + scale_fill_manual(values = type_colors) +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 10,color = "black" ),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 8,color = "black"),
        axis.text.x = element_text(angle = 0)) +
  stat_compare_means(comparisons = compar,label.y = 0.5, label = "p.signif", vjust = 2 )+
  scale_color_manual(values = group_colors)
ggsave(file.path(outdir2,paste0('TF_IRF1_group_violin.pdf')), p, height = 3,width = 4)

#### 8.Fig6H ----
library(monocle)
seurat.obj <- qs::qread('./files/named_seurat.qs')
markers <- FindAllMarkers(seurat.obj, assay = "integrated", only.pos = T)

data <- as(as.matrix(seurat.obj@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat.obj@meta.data)
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
top_markers <- markers %>% 
  filter(p_val <= 0.05) %>% 
  arrange("avg_log2FC") %>% 
  group_by(cluster) %>% 
  top_n(40, wt = avg_log2FC) %>% 
  pull(gene)
ordergene <- top_markers
cds <- setOrderingFilter(cds, ordergene)
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree', verbose = F)
cds <- orderCells(cds)
saveRDS(cds, file.path(outdir, "cds.rds"))

p1 <- plot_cell_trajectory(cds, color_by = "Pseudotime",
                           cell_size = 0.1, show_branch_points = F, size = 1, show_backbone = TRUE)
ggsave(file.path(outdir2, 'Pseudotime.pdf'), p1, width = 3.5, height = 4)
p5 <- plot_cell_trajectory(cds, color_by = "cell_type",
                           cell_size = 0.1, show_branch_points = F, size = 1, show_backbone = TRUE)+
  scale_color_manual(values = type_colors)
ggsave(file.path(outdir2, 'cell_type.pdf'), p5, width = 3.5, height = 4)

#### 9.Fig6I ----
## heatmap
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
diff <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~cell_type", cores = 40)
diff_genes <- subset(diff, qval < 0.001)
diff_genes <- diff_genes[order(diff_genes$qval, decreasing = F),]
diff_genes <- diff_genes[1:1000,]
ordergene <- diff_genes$gene_short_name
ordergene <- unique(ordergene)
Time_diff <- differentialGeneTest(cds[ordergene,],  cores = 1, fullModelFormulaStr = "~sm.ns(Pseudotime)")

rownames(Time_diff) <- Time_diff$gene_short_name
pdf(file=file.path(outdir2, "non-branched_heatmap2.pdf"), width= 5, height= 5)
plot_pseudotime_heatmap(
  cds[row.names(subset(Time_diff, qval < 0.05)),],
  num_clusters = 5, cores = 1,
  use_gene_short_name = T,
  show_rownames = F) +
  theme(plot.margin=unit(rep(1,4),'lines'))
dev.off()

## density
library(ggridges)
library(RColorBrewer)
library(scales)

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
ggsave(file.path(outdir2, "celltype_tmp1.pdf"), width = 15, height = 3, units = "cm")

#### 10.Fig6J ----
s_genes <- c('KLRC1', 'IRF7', 'IL32', 'ID2')
for(i in s_genes){
  p <- plot_genes_in_pseudotime(cds[i,], color_by = "group")+
    scale_color_manual(values = group_colors)
  df <- p$data
  df$group <- factor(df$group, levels = c('Normal', 'LM', 'PE'))
  ggplot(df, aes(Pseudotime, adjusted_expression)) + 
    geom_point(aes(colour = group), size = 0.5 )+
    geom_smooth(aes(colour = group),method = "loess", se = FALSE)+
    labs(title=i, y='Adjusted Expression')+theme_bw()+
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.text = element_text(color = "black", size = 6),
          axis.title.x = element_text(color = "black", size = 8),
          plot.title = element_text(color = "black", size = 8, face = 'italic'))+
    scale_color_manual(values = group_colors)
  
  ggsave(file.path(outdir2, paste0(i,"_geomplot.pdf")), last_plot(), height = 3, width = 3)
}
