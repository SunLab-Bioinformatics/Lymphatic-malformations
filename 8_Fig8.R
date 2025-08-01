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

setwd('~/LM/Fig8/')
outdir <- "./files/"
dir.create(outdir,recursive = T)
outdir2 <- "./plots/"
dir.create(outdir2,recursive = T)

#### 1.Fig8A ----
library(ggVennDiagram)
library(ggvenn)
library(readxl)

DEG_LM_N <- readRDS("~/LM/Fig2/DEG/DEG_LM_N.rds")
DEG_up_LM <- DEG_LM_N |> filter(change == 'UP' & p_val < 0.05) |> 
  filter(cluster %in% c('CD14+ Monocytes','CD14+CD16+ Monocytes','CD16+ Monocytes','Dendritic cells')) |>
  pull(gene) |> unique()

DEG_PE_N <- readRDS("~/LM/Fig2/DEG/DEG_PE_N.rds")
DEG_up_PE <- DEG_PE_N |> filter(change == 'UP' & p_val < 0.05) |>
  filter(cluster %in% c('CD14+ Monocytes','CD14+CD16+ Monocytes','CD16+ Monocytes','Dendritic cells')) |>
  pull(gene) |> unique()

deg_LM <- intersect(DEG_up_LM, DEG_up_PE)
gene_list <- list("LM vs Normal" = DEG_up_LM,
                  "PE vs Normal" = DEG_up_PE)
ggvenn(gene_list,show_percentage = F,
       stroke_color = "white",
       fill_color = c("#156077","#f46f20"),
       fill_alpha = 0.8,
       set_name_color = "black",
       text_size = 10,
       text_color = "black")
ggsave(paste0(outdir2,"venn_LM_PE.pdf"), last_plot(), height = 5, width = 6)

#### 2.Fig8B ----
library(ggraph)
library(tidydr)
library(ggthemes)
library(AUCell)

seurat_obj <- qs:qread('~/LM/merge/files/Mye.qs')
cells_ranking <- AUCell_buildRankings(seurat_obj@assays$RNA@counts)
cells_AUC <- AUCell_calcAUC(deg_LM, cells_ranking, aucMaxRank=nrow(cells_ranking)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
seurat_obj$AUC  <- AUCell_socre
seurat_obj$AUC[seurat_obj$AUC> 0.50] <- 0.50
seurat_obj$AUC[seurat_obj$AUC< 0.20] <- 0.20

a <- c("Normal","LM","PE")
data0 <- data.frame(seurat_obj@meta.data, seurat_obj@reductions$umap@cell.embeddings)
for (i in 1:3) {
  data <- subset(data0, group == a[i])
  assign(paste0("p",i), ggplot(data,  aes(UMAP_1, UMAP_2, color=AUC)) + 
           geom_point(size=0.1) + labs(title = a[i])+
           scale_color_gradient(low = "#68217A",  high = "#fbe183")+ 
           theme_dr(xlength = 0.2, ylength = 0.2,
                    arrow = grid::arrow(length = unit(0.1, "inches"), ends = 'last', type = "closed")) +
           guides(color=guide_legend(title = "SMEP score",override.aes = list(size = 5),reverse = T))+
           theme(
             aspect.ratio = 1, panel.grid = element_blank(),
             plot.title = element_text(hjust=0.5,vjust = -1)))
  ggsave(paste0(outdir2, a[i], '_SMEP_score_UMAP.pdf'), last_plot(),width = 4, height = 4)
}

p <- p1 + p2 + p3 + plot_layout(guides = 'collect')
ggsave(file.path(outdir2, 'SMEP_score_UMAP.pdf'), last_plot(),width = 12, height = 4)

#### 3.Fig8E ----
immune_gene <- read_xlsx("./files/immune_genes.xlsx",col_names = F) |> pull(...1)
vessel_gene <- read_xlsx("./files/vessel_genes.xlsx",col_names = F)|> pull(...1)
inter_dis <- intersect(immune_gene, vessel_gene)
disease_list <- list("vessel" = vessel_gene,
                     "immune" = immune_gene)
ggvenn(disease_list,show_percentage = F,
       stroke_color = "white",
       fill_color = c("#4da0a0","#9b3a74"),
       fill_alpha = 0.8,
       set_name_color = "black",
       text_size = 10,
       text_color = "black")
ggsave(paste0(outdir2,"venn_disease.pdf"), last_plot(), height = 5, width = 6)

#### 4.Fig8F ----
cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@counts) 
cells_AUC <- AUCell_calcAUC(inter_dis, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
seurat_obj$AUC  <- AUCell_socre

data <- data.frame(seurat_obj@meta.data, seurat_obj@reductions$umap@cell.embeddings)
data_N <- subset(data, group == "Normal")
compar <- list(c('Normal', 'LM'),
               c('Normal', 'PE'),
               c('LM', 'PE'))
ggviolin(data, x='group',y="AUC",
         xlab = '', ylab = "Shared_VID score",
         color ='group', add = 'boxplot',
         add.params = list(fill='white', width=0.1))+
  scale_color_manual(values = c("#238b3a","#f9b920","#cd2125"))+
  geom_hline(yintercept = median(data_N$AUC),
             linetype = "dashed", size = 0.5, color = "black")+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5), axis.ticks.x = element_blank(),
        legend.position = 'right') + NoLegend() +
  stat_compare_means(comparisons = compar, label = "p.signif",)
ggsave(file.path(outdir2, 'Shared_VID score vlnplot.pdf'), last_plot(),width = 6, height = 4)

#### 5.Fig8H ----
library(Nebulosa)
library(scCustomize)

features <- c("S100A8","FCGR1A","CD163")

for (i in features) {
  outdir3 <- paste0(outdir2,'coexp/')
  dir.create(outdir3, recursive = T)

  p1 <- Plot_Density_Joint_Only(seurat_obj, features = c("CD14","FCGR3A",i),
                                custom_palette = BlueAndRed()) +
    theme(plot.title = element_text(face = "italic"))
  ggsave(paste0(outdir3,i,"_coexp.pdf"), p1, height = 4.5, width = 5)
}

#### 6.Fig8I ----
library(msigdbr)

cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@counts) 
category <- "C5"
Dataset <- msigdbr(species = "Homo sapiens", category = category)
Dataset <- Dataset %>% filter(str_detect(Dataset$gs_name,"GOBP_"))

dataset <- Dataset
geneSets <- lapply(unique(dataset$gs_name), 
                   function(x){dataset$gene_symbol[dataset$gs_name == x]})
names(geneSets) <- unique(dataset$gs_name)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
colnames(AUCell_socre) <- colnames(AUCell_socre)|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
terms_GO <- c('Regulation of inflammatory response',
              'Acute inflammatory response',
              'Chronic inflammatory response',
              'Cytokine production involved in inflammatory response',
              'Leukocyte migration involved in inflammatory response',
              'Production of molecular mediator involved in inflammatory response')
object <- AUCell_socre[ ,colnames(AUCell_socre) %in%terms_GO]
a <- t(as.data.frame(seurat_obj@assays[["RNA"]]@data[features,]))

## SMEP
cells_ranking <- AUCell_buildRankings(seurat_obj@assays$RNA@counts)
cells_AUC <- AUCell_calcAUC(deg_LM, cells_ranking, aucMaxRank=nrow(cells_ranking)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
seurat_obj$SMEP  <- AUCell_socre

## Shared_VID
cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@counts) 
cells_AUC <- AUCell_calcAUC(inter_dis, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
seurat_obj$Shared_VID  <- AUCell_socre

object <- cbind(object, a, seurat_obj@meta.data[,c("SMEP","Shared_VID")])
condition <- colnames(object)

library(corrplot)
cor_obj <- cor(object)
col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))

pdf(paste0(outdir2,"cor_all_mye_test.pdf"), width = 14,height = 14)
p <- corrplot(cor_obj, method = "color", col = col(200),
              diag = FALSE, addCoef.col = "black", 
              type = "upper", 
              tl.pos="td", tl.cex=1, tl.col="black",tl.offset=0.5,tl.srt = 45)
dev.off()
