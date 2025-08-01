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

setwd('~/LM/Fig7/')
outdir <- "./files/"
dir.create(outdir,recursive = T)
outdir2 <- "./plots/"
dir.create(outdir2,recursive = T)

#### 1.Fig7A ----
cellchat.1 <- readRDS("./results/Normal/cellchat.rds")
cellchat.2 <- readRDS("./results/LM/cellchat.rds")

cellchat.1 <- netAnalysis_computeCentrality(cellchat.1)
cellchat.2 <- netAnalysis_computeCentrality(cellchat.2)
object.list <- list(Normal = cellchat.1, LM = cellchat.2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = "structural")
cellchat <- netClustering(cellchat, type = "structural")

sources = c(1:10)
targets = c(6)
pdf(file = paste0(outdir2, 'MHC-I NK.pdf'),width=6, height= 4.5)
p <- netVisual_bubble(cellchat, sources.use = sources, targets.use = targets,  
                      signaling = c("MIF"),
                      comparison = c(1, 2), angle.x = 90)
dev.off()

#### 2.Fig7B ----
seurat_obj <- readRDS('~/LM/Fig6/files/named_seurat.qs')
selected_genes <- c('KLRC1', 'KLRB1', 'KLRD1', 'KLRK1', 'KLRC2', 'KLRC3','KIR2DL1', 'KIR3DL1') 
selected_genes <- selected_genes[!duplicated(selected_genes[selected_genes %in% rownames(seurat_obj)])]
VlnPlot(seurat_obj, features = selected_genes,
        cols = c("#d8443c","#9f5691","#2b9b81","#fe9b00","#92c051",
                 "#e6a2a6","#fbe183","#1f6e9c","#aa7aa1"),
        group.by = 'group',flip = T,stack = T)+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 10,angle = 0,hjust = 0.5),
        plot.title = element_text(hjust = 0.5,size = 10))+NoLegend()
ggsave(file.path(outdir2, "LR gene expr group.pdf"), last_plot(), height=4, width=3)

#### 3.Fig7C ----
seurat_obj0 <- readRDS('~/LM/Fig6/files/named_seurat.qs')
seurat_obj <- seurat_obj0[,seurat_obj0[['cell_type']]== "NK CD56lowCD16high"]
category <- "C5"
case = 'all_'

Dataset <- msigdbr(species = "Homo sapiens", category = category)
Dataset <- Dataset %>% filter(str_detect(Dataset$gs_name,"GOBP_"))
Dataset$gs_name2 <- Dataset$gs_name|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ")|>
  str_to_sentence()
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
  'Response to interferon gamma')
dataset <- Dataset %>% filter(Dataset$gs_name2 %in% s_term)
length(unique(dataset$gs_name))
geneSets <- lapply(unique(dataset$gs_name), 
                   function(x){dataset$gene_symbol[dataset$gs_name == x]})
names(geneSets) <- unique(dataset$gs_name2)
cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@data)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
saveRDS(AUCell_socre, file.path(outdir, paste0(case,category,"AUCscore.rds")))

method = "pearson"
pathway_cor_test <- data.frame(
  feature_x=NULL,
  feature_y=NULL,
  p_value=NULL,
  estimate=NULL,
  num=NULL,
  method=NULL)
for (feature_x in colnames(expr)) {
  for (feature_y in colnames(AUCell_socre)) {
    cor.test.res <-
      cor.test(adata[, feature_x],
               adata[, feature_y],
               method = method)
    p.value <- cor.test.res$p.value
    estimate <- cor.test.res$estimate
    pathway_cor_test <-
      rbind(
        pathway_cor_test,
        data.frame(
          feature_x = feature_x,
          feature_y = feature_y,
          p_value = p.value,
          estimate = estimate,
          num = nrow(adata),
          method = method))
  }
}
saveRDS(pathway_cor_test, file.path(outdir, paste0(case,category,"gene_pathway_cor.rds")))

pathway_cor_test <- readRDS(file.path(outdir, paste0(case,category,"gene_pathway_cor.rds")))
selected_genes <-  c('KLRK1','KLRC2','KLRC3','KLRB1','KLRC1') 
pathway_cor_test$text=ifelse(pathway_cor_test$p_value<0.001,"***",
                             ifelse(pathway_cor_test$p_value<0.01,"**",
                                    ifelse(pathway_cor_test$p_value<0.05,"*","")))
pathway_cor_test$estimate=as.numeric(pathway_cor_test$estimate)
pathway_cor_test$feature_x <- factor(pathway_cor_test$feature_x, levels = selected_genes)
pathway_cor_test2 <- pathway_cor_test |> 
  filter( abs(estimate) >0.1 & estimate < 0.9  &  p_value< 0.05) |>
  arrange(desc(estimate))
pathway_cor_test2$LogP <- -log10(pathway_cor_test2$p_value)
pathway_cor_test2$LogP[pathway_cor_test2$LogP==Inf] <- 100
pathway_cor_test2$LogP[pathway_cor_test2$LogP > 100] <- 100

selected_features <- c('Antibody dependent cellular cytotoxicity',
                       'Positive regulation of natural killer cell activation',
                       'Positive regulation of natural killer cell proliferation',
                       'Natural killer cell proliferation',
                       'Regulation of natural killer cell activation',
                       'Protection from natural killer cell mediated cytotoxicity',
                       'Natural killer cell activation involved in immune response',
                       'Positive regulation of natural killer cell mediated cytotoxicity',
                       'Natural killer cell degranulation')
  
pdf(file = file.path(outdir2,paste0(case, "dotplot_cor_s.pdf")), width=4, height=2.5)
pathway_cor_test2 %>%
  filter(feature_x %in% selected_genes, feature_y %in% selected_features) %>%
  mutate(feature_x = fct_relevel(feature_x, selected_genes), feature_y = fct_relevel(feature_y, selected_features)) %>%
  catdotplot(x = feature_x, y = feature_y, size = LogP, color = estimate, title = '') +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6, face = "italic"),
        panel.grid = element_line(size = 0.2, color = "lightgrey"))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
dev.off()

#### 4.Fig7D ----
cellchat <- readRDS(paste0(outdir, 'merged_cellchat_PE_LM.rds'))
pdf(file = paste0(outdir2, "DiffInteraction_hratmap.pdf"),width= 7,height= 4)
netVisual_heatmap(cellchat, color.use = type_colors)
dev.off()

#### 5.Fig7E ----
pdf(file = paste0(outdir2,  "B cells :Signaling changes2.pdf"),width=20,height=20)
gg <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "B cells")+theme(aspect.ratio = 1)
patchwork::wrap_plots(plots = list(gg))
dev.off()

#### 6.Fig7F ----
pdf(file = paste0(outdir2, 'BC DOWN Change Chord.pdf'),width=5, height= 5)
netVisual_chord_gene(object.list[[1]], sources.use = sources, targets.use = targets, 
                     color.use = type_colors,
                     slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, show.legend = T,
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()

#### 7.Fig7G ----
library(CellChat)
library(patchwork)
library(Seurat)
library(qs)
library(future)
options(future.globals.maxSize = 100000 * 1024^5)
getOption("future.globals.maxSize")
options(future.plan="multiprocess",mc.cores = parallel::detectCores() - 1L)

seurat_obj <- qread('~/LM/merge/files/all_subtypes_RNA.qs')
case <- "all"
outdir <- paste0("./cellchat/files/")
outdir2 <- paste0("./cellchat/plots/")

expr <- seurat_obj@assays$RNA@data
data.input <- expr
meta <- as.data.frame(Idents(seurat_obj))
colnames(meta) <- "labels"
unique(meta$labels)
cellchat <-createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 25) 
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) 
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
saveRDS(cellchat, file = paste0(outdir, "./cellchat/cellchat.rds"))

## plot
LR_data <- p[["data"]][,1:4]
LR_data <- LR_data[!is.na(LR_data$ligand), ] 

LR_data$receptor <- gsub('CD74_CD44', 'CD44', LR_data$receptor)
LR_data$receptor <- gsub('CD74_CXCR4', 'CXCR4', LR_data$receptor)
LR_data$receptor <- gsub('ITGAL_ITGB2', 'ITGB2', LR_data$receptor)
LR_data$receptor <- gsub('ITGAM_ITGB2', 'ITGB2', LR_data$receptor)
LR_data$receptor <- gsub('ITGAX_ITGB2', 'ITGB2', LR_data$receptor)
LR_data_name <- data.frame(c(names(table(LR_data$ligand)), names(table(LR_data$receptor))),
                           c(1:length(table(LR_data$ligand)),1:length(table(LR_data$receptor))))

data <- data.frame(sou=LR_data$ligand,
                   x1=rep(2,length(LR_data$ligand)),
                   net_y1='',
                   tar=LR_data$receptor,
                   x2=rep(3,length(LR_data$receptor)),
                   net_y2='')
colnames(LR_data_name) <- c('sou', 'number')
match_indices <- match(data$sou, LR_data_name$sou)
data$net_y1 <- LR_data_name$num[match_indices]
colnames(LR_data_name) <- c('tar', 'number')
match_indices <- match(data$tar, LR_data_name$tar)
data$net_y2 <- LR_data_name$num[match_indices]
data <- unique(data)
write.table(data, file = file.path(outdir, "data_dotplot.csv"), quote = F, sep = ",", row.names = T)

##
features_1 <- c(LR_data_name$genes[1:9])
p1 <- DotPlot(object = seurat_obj[, seurat_obj[["cell_type"]] == c('NK CD56highCD16low', 'NK CD56lowCD16high',         
                                                                   'CD14+ Monocytes',  'CD14+CD16+ Monocytes',    
                                                                   'CD16+ Monocytes',  'Dendritic cells' )], 
              features = features_1, assay = "RNA") +
  guides(color = guide_colorbar(title = 'Average  Expression')) +
  scale_color_gradientn(colours = c('white','grey',"black"))+ 
  coord_flip() +
  theme(axis.title = element_blank(),
        legend.position = "left",
        plot.margin = unit(c(0.5,0,0.5,0.5), 'cm'),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.2),
        panel.border = element_rect(color = "black",fill=NA, size = 0.1),
        panel.grid.major = element_line(size = 0.1,color = 'black',linetype = 2),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.1)) +
  scale_x_discrete(position = "top") + scale_y_discrete(position = "right")
ggsave(file.path(outdir2, 'left_target.pdf'), p1, width = 4.5, height = 5.5)

##
features_2 <- c(LR_data_name$genes[10:19])
p3 <- DotPlot(object = seurat_obj, features = features_2, assay = "RNA") +
  guides(color = guide_colorbar(title = 'Average  Expression')) +
  scale_color_gradientn(colours = c('white','grey',"black"))+ 
  coord_flip() +
  theme(axis.title = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0), 'cm'),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.2),
        panel.border = element_rect(color = "black",fill=NA, size = 0.1),
        panel.grid.major = element_line(size = 0.1,color = 'black',linetype = 2),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.1)) +
  scale_y_discrete(position = "right")  + 
  NoLegend()
ggsave(file.path(outdir2, 'right_source.pdf'), p3, width = 7, height = 6)

##
p2 <- ggplot(data)+
  geom_segment(aes(x1,net_y1,xend=x2,yend=net_y2),
               size=0.5,color=c('black','black','black','black','#CC0033',
                                'black','black','black','black','black',
                                'black','black','black'))+
  geom_point(aes(x=x1,y=net_y1),size=1, fill="#3cb346", color="#3cb346",
             stroke=1, shape = 21)+
  geom_point(aes(x=x2,y=net_y2),size=1, fill="#44c1f0", color="#44c1f0",
             stroke=1, shape = 21)+
  scale_y_continuous(limits = c(1, 10),expand = expansion(add=c(0.5,0.7)))+
  scale_x_continuous(expand = expansion(0,0.1))+
  theme_void()
ggsave(file.path(outdir2, 'mid_line.pdf'), p2, width = 2, height = 5)


