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
source('~/Collection/code/plot.R')
source('~/Collection/code/scrna-seq.R')

setwd('~/LM/Fig3/')
outdir <- "./files/"
dir.create(outdir,recursive = T)
outdir2 <- "./plots/"
dir.create(outdir2,recursive = T)

type_colors <- c('#f46f20', "#6a73cf",'#0eb0c8', '#476D87')
group_colors <- c( "#228B22", "#ffbc14", "firebrick3")
#### 1.Fig3A ----
seurat_obj <- qs:qread('~/LM/merge/files/impor_seurat.qs')
celltypes <- c('CD14+ Monocytes','CD14+CD16+ Monocytes','CD16+ Monocytes','Dendritic cells')
Mye <- subset(seurat_obj, idents = c('CD14+ Monocytes','CD14+CD16+ Monocytes','CD16+ Monocytes','Dendritic cells'))
seurat_obj$cell_type <- factor(seurat_obj$cell_type, levels = names(table(Idents(seurat_obj))))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
qs::qsave(seurat_obj, '~/LM/merge/files/Mye.qs')

DimPlot(seurat_obj)+
  theme_void() +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), ends = 'last', type = "closed")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs( x= "UMAP 1",y= "UMAP 2")+
  theme(plot.margin = margin(5.5,15,5.5,5.5),
        panel.grid = element_blank(),aspect.ratio = 1) +
  scale_color_manual(values = type_colors) +
  scale_fill_manual(values = type_colors)
ggsave(paste0(outdir2,'Mye_Dimplot.pdf'), last_plot(), height=5, width=5)

#### 2.Fig3B ----
cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@data) 
category <- "C5"
Dataset <- msigdbr(species = "Homo sapiens", category = category)
Dataset <- Dataset %>% filter(str_detect(Dataset$gs_name,"GOBP_"))
geneSets <- lapply(unique(Dataset$gs_name),
                   function(x){Dataset$gene_symbol[Dataset$gs_name == x]})
names(geneSets) <- unique(Dataset$gs_name)
s_term <- c('GOBP_PHAGOCYTOSIS_RECOGNITION', 'GOBP_INFLAMMATORY_RESPONSE')
for(i in s_term){
  seurat_obj$group <- factor(seurat_obj$cell_type)#
  s_geneSets <- geneSets[[i]]
  geneSets[[i]]
  cells_AUC <- AUCell_calcAUC(geneSets[[i]], cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
  aucs <- as.numeric(getAUC(cells_AUC))
  seurat_obj$AUC  <- aucs
  
  meta.data <- seurat_obj@meta.data[, c("cell_type", "AUC")]
  meta.data <- meta.data %>% reshape2::melt(id.vars = c("cell_type"))
  colnames(meta.data)[2:3] <- c("signature", "score")
  p2 <- meta.data %>%
    ggplot(aes(x = cell_type, y = score, fill = '', color =cell_type)) +
    geom_violin(fill = 'NA') +
    geom_boxplot(data = meta.data,
                 aes(x = cell_type, y = score),
                 width = 0.1,  alpha = 0.1)+
    geom_hline(yintercept = mean(meta.data$score),
               linetype = "dashed", size = 0.5, color = "darkgrey")+
    scale_y_continuous(i |>
                         str_replace_all("GOBP_", " ") |>
                         str_replace_all("_", " ")|>
                         str_wrap(width = 40)|>
                         str_to_sentence()) +
    theme_test() + scale_fill_manual(values = type_colors) +
    theme(legend.position = "none",
          axis.title.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.text.y = element_text(color = 'black'),
          axis.text.x = element_blank()) +
    stat_compare_means(comparisons = list(c('CD14+ Monocytes', 'CD16+ Monocytes')),
                       label = "p.signif", vjust = 1.5 )+
    scale_color_manual(values = type_colors)
  ggsave(file.path(outdir2,paste0(i,'_violin_celltype.pdf')), p2, height=2.5, width= 3.5)
}

#### 3.Fig3C ----
library(clusterProfiler)
library(org.Hs.eg.db)

case <- 'LM vs Normal '
for (i in c('CD14+ Monocytes','CD14+CD16+ Monocytes','CD16+ Monocytes','Dendritic cells')) {
  path <- paste0('./deg/LM_N/',i,'/deg.rds')
  LM_N <- readRDS(path)|> filter(abs(avg_log2FC) > 0.25 & p_val < 0.05)
  ## UP
  deg <- LM_N[LM_N$change %in% 'UP', ] |> pull('gene')
  bp1 <-
    enrichGO(
      deg,
      OrgDb = org.Hs.eg.db,
      keyType = 'SYMBOL',
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2)
  term1 <- bp1@result
  write.table(term1, file.path(outdir, paste0(i, '_', case,'UP_GOBP.csv')), quote = F, sep = ",", row.names = F)
  ## DOWN
  deg <- LM_N[LM_N$change %in% 'DOWN', ] |> pull('gene')
  bp2 <-
    enrichGO(
      deg,
      OrgDb = org.Hs.eg.db,
      keyType = 'SYMBOL',
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2)
  term2 <- bp2@result
  write.table(term2, file.path(outdir, paste0(i, '_', case,'DOWN_GOBP.csv')), quote = F, sep = ",",row.names = F)
}

## load data
UP_A <- read.csv(paste0(outdir, 'CD14+ Monocytes_LM vs Normal UP_GOBP.csv'))
UP_A$change <- 'UP'
UP_A$celltype <- 'CD14+ Monocytes'
DOWN_A <-  read.csv(paste0(outdir, 'CD14+ Monocytes_LM vs Normal DOWN_GOBP.csv'))
DOWN_A$change <- 'DOWN'
DOWN_A$celltype <- 'CD14+ Monocytes'

UP_B <- read.csv(paste0(outdir, 'CD14+CD16+ Monocytes_LM vs Normal UP_GOBP.csv'))
DOWN_B <-  read.csv(paste0(outdir, 'CD14+CD16+ Monocytes_LM vs Normal DOWN_GOBP.csv'))
UP_B$change <- 'UP'
UP_B$celltype <- 'CD14+CD16+ Monocytes'
DOWN_B$change <- 'DOWN'
DOWN_B$celltype <- 'CD14+CD16+ Monocytes'

UP_C <- read.csv(paste0(outdir, 'CD16+ Monocytes_LM vs Normal UP_GOBP.csv'))
DOWN_C <- read.csv(paste0(outdir, 'CD16+ Monocytes_LM vs Normal DOWN_GOBP.csv'))
UP_C$change <- 'UP'
UP_C$celltype <- 'CD16+ Monocytes'
DOWN_C$change <- 'DOWN'
DOWN_C$celltype <- 'CD16+ Monocytes'

UP_D <- read.csv(paste0(outdir, 'Dendritic cells_LM vs Normal UP_GOBP.csv'))
DOWN_D <- read.csv(paste0(outdir, 'Dendritic cells_LM vs Normal DOWN_GOBP.csv'))
UP_D$change <- 'UP'
UP_D$celltype <- 'Dendritic cells'
DOWN_D$change <- 'DOWN'
DOWN_D$celltype <- 'Dendritic cells'

UP_data <- rbind(UP_A, UP_B, UP_C, UP_D)
UP_data <- na.omit(UP_data)
UP_data$pvalue <- as.numeric(UP_data$pvalue)

DOWN_data <- rbind(DOWN_A, DOWN_B, DOWN_C, DOWN_D)
DOWN_data <- na.omit(DOWN_data)
DOWN_data$pvalue <- as.numeric(DOWN_data$pvalue)

up_term <- c(
  'phagocytosis',
  'regulation of phagocytosis',
  'positive regulation of inflammatory response',
  'regulation of inflammatory response',
  'defense response to bacterium', 
  'response to insulin',
  'response to type II interferon',
  'positive regulation of secretion',
  'cytokine-mediated signaling pathway',
  'integrated stress response signaling',
  'cellular response to tumor necrosis factor',
  'cellular response to interleukin-1',
  'chronic inflammatory response',
  'positive regulation of NF-kappaB transcription factor activity',
  'innate immune response-activating signaling pathway',
  'regulation of innate immune response')
down_term <- c(
  'antigen processing and presentation',
  'MHC class II protein complex assembly',
  'immunoglobulin mediated immune response',
  'positive regulation of T cell activation',
  'positive regulation of cell adhesion',
  'type II interferon production',
  'cellular response to interleukin-4',
  'ubiquitin-dependent ERAD pathway',
  'cell maturation',
  'positive regulation of regulatory T cell differentiation',
  'tolerance induction',
  'cell-cell adhesion mediated by integrin',
  'cytolysis',
  'leukocyte mediated cytotoxicity',
  'positive regulation of interleukin-12 production')
UP_data2 <- UP_data |> filter(Description %in% up_term)
DOWN_data2 <- DOWN_data |> filter(Description %in% down_term)

UP_data2$pvalue <- as.numeric(UP_data2$pvalue)
UP_data2$Count <- as.numeric(UP_data2$Count)
UP_data2$Description <- factor(UP_data2$Description, levels = up_term)
DOWN_data2$pvalue <- as.numeric(DOWN_data2$pvalue)
DOWN_data2$Count <- as.numeric(DOWN_data2$Count)
DOWN_data2$Description <- factor(DOWN_data2$Description, levels = down_term)

p1 <- UP_data2 %>%
  catdotplot(
    x = celltype,
    y = Description,
    size = Count ,
    color = -log10(pvalue) ,
    title = 'UP-terms: LM vs Normal') +
  coord_fixed() +
  theme(
    axis.text.x  = element_text(angle = 30, vjust = 1, hjust = 1),
    panel.grid = element_line(size = 0.2, color = "lightgrey"))+
  scale_color_gradient(low = "#ffdde1", high = "#ee6470")
ggsave(file.path(outdir2, paste0(case, 'UP terms.pdf')), p1, height=4, width=5)

DOWN_data2$pvalue[DOWN_data2$pvalue< 1.399863e-8] <- 1.399863e-8
p2 <- DOWN_data2 %>%
  catdotplot(
    x = celltype,
    y = Description,
    size = Count ,
    color = -log10(pvalue) ,
    title = 'DOWN-terms: LM vs Normal') +
  coord_fixed() +
  theme(
    axis.text.x  = element_text(angle = 30, vjust = 1, hjust = 1),
    panel.grid = element_line(size = 0.2, color = "lightgrey"))+
  scale_color_gradient(low = "#a0ddff", high = "#00a6e1")
ggsave(file.path(outdir2, paste0(case, 'DOWN terms.pdf')), p2, height=4, width=5)

#### 4.Fig3D ----
## CD14+CD16+ Monocytes
gsea.input <- readRDS('./deg/LM_N/CD14+CD16+ Monocytes/deg.rds')
category <- "H"
genesets <- msigdbr(species = "Homo sapiens", category = category)
genesets <- subset(genesets, select = c("gs_name", "gene_symbol"))
rownames(gsea.input) <- gsea.input$gene
gsea.input <- gsea.input[order(gsea.input$avg_log2FC, decreasing = T),]
genelist <- structure(gsea.input$avg_log2FC,names=rownames(gsea.input))
res <- GSEA( genelist, TERM2GENE = genesets , pvalueCutoff=1, eps=0)
gsea_result<-as.data.frame(res@result)
gsea_result <- gsea_result[order(gsea_result$pvalue, decreasing = F),]
gsea_result$Description<- gsub('_',' ',gsea_result$Description)
gsea_result$Description<- str_to_sentence(gsea_result$Description)
saveRDS(res, file = paste0(outdir,"H_gsea.rds"))

res@result$Description <- res@result$Description|>
  str_replace_all("HALLMARK_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
p <- res |> cat_gseaplot(
  'HALLMARK_INFLAMMATORY_RESPONSE',
  subplots = c(1, 2),
  pvalue_table = T,
  title = 'HALLMARK_INFLAMMATORY_RESPONSE'|> 
    str_replace_all("HALLMARK_", "") |>
    str_replace_all("_", " ") |>
    str_to_sentence())
ggsave(file.path(outdir2, paste0(i,'.pdf')), p, height = 2.5, width = 3.5)

## CD14+ Monocytes
gsea.input <- readRDS('./deg/LM_N/CD14+ Monocytes/deg.rds')
category <- "C5"
genesets <- msigdbr(species = "Homo sapiens", category = category)
genesets <- subset(genesets, select = c("gs_name", "gene_symbol"))
rownames(gsea.input) <- gsea.input$gene
gsea.input <- gsea.input[order(gsea.input$avg_log2FC, decreasing = T),]
genelist <- structure(gsea.input$avg_log2FC,names=rownames(gsea.input))
res <- GSEA( genelist, TERM2GENE = genesets , pvalueCutoff=1, eps=0)
gsea_result<-as.data.frame(res@result)
gsea_result <- gsea_result[order(gsea_result$pvalue, decreasing = F),]
gsea_result$Description<- gsub('_',' ',gsea_result$Description)
gsea_result$Description<- str_to_sentence(gsea_result$Description)
saveRDS(res, file = paste0(outdir,"C5_gsea.rds"))

res@result$Description <- res@result$Description|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
p <- res |> cat_gseaplot(
  'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_ANTIGEN',
  subplots = c(1, 2),
  pvalue_table = T,
  title = 'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_ANTIGEN'|> 
    str_replace_all("GOBP_", "") |>
    str_replace_all("_", " ") |>
    str_to_sentence())
ggsave(file.path(outdir2, paste0(i,'.pdf')), p, height = 2.5, width = 3.5)

#### 5.Fig3F ----
seurat_obj <- qs::qread('~/LM/merge/files/Mye.qs')
seurat_obj <- subset(seurat_obj, idents = "CD14+ Monocytes")
seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                              species = "human",
                              outdir = "./TF/results/")
seurat_obj$orig.ident <- factor(seurat_obj$orig.ident, 
                                levels = c("BLBL804" ,"BLJK804" ,"CNLHL420",
                                           "BFYY331" ,"BLJJ329","BHLJ1022","BJZY1103","BMYY203",
                                           "CFYY331" ,"CLJJ329", "CGCJ1028" ))
adata <- read.table( "./TF/results/tfs_target.tsv", sep = ',', header = T)
tf.seleted <- c('REL','RFX5')
adata2 <- adata |> filter(tf %in% tf.seleted)
deg <- readRDS('./deg/LM_N/CD14+ Monocytes/deg.rds')
deg_down <- deg |> filter(p_val <= 0.05 & avg_log2FC <= -0.10)
inter_G <- intersect(deg_down$gene, adata2$target_gene)
TF_G1 <- deg[unique(adata2$target_gene),]
colnames(TF_G1)[7] <- 'target_gene'
TF_G2 <- merge(adata2, TF_G1, by = 'target_gene',)
TF_G3 <- TF_G2[,c("target_gene", "symbol",  "tf",  "p_val", "avg_log2FC", "change")]
write.table(TF_G3, file = file.path(outdir, 'TFs DEGs.csv'), quote = F, sep = ",", row.names = F)

##
regulonAUC <- importAUCfromText("auc_mtx.csv")
regulonAUC <-
  regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
colnames(regulonAUC) <- gsub(',','', colnames(regulonAUC))
cellInfo <- seurat_obj@meta.data
seurat_obj[["scenic"]] <- CreateAssayObject(counts = t(getAUC(regulonAUC)))
seurat_obj2 <- subset(seurat_obj,subset= group !='PE')
DefaultAssay(seurat_obj2) <- "scenic"
tf.seleted <- c('REL(+)', 'RFX5(+)')
expr <- t(as.data.frame(seurat_obj2@assays$scenic@data[tf.seleted,]))
features <- c('REL', 'RFX5')
colnames(expr) <- features
write.table(expr, file=file.path(outdir, "selected TF activity.csv"),
            quote = F, sep = ",", row.names = T)
df <- cbind(expr, seurat_obj2@meta.data)
for(i in 1:length(features)){
  ggviolin(df, "group",  features[i], 
           color = "group",fill = "group",
           width = 1, palette = c("#228B22", "#ffbc14"),
           add.params = list(fill = "white"))+
    geom_boxplot(width = 0.1, fill = "white", color = 'black',outlier.shape = NA) +
    theme(axis.title.x = element_blank()) + NoLegend()+ labs(y = paste0(features[i],' TF activity'))+
    stat_compare_means(comparisons = list(c('Normal', 'LM')), label = "p.signif", vjust = 2 )
  ggsave(file.path(outdir2, paste0(features[i], '_activity_group.pdf')),
         ggplot2::last_plot(), width = 3,height = 3)
}
