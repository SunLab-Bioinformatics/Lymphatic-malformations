library(scRepertoire)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(qs)

setwd('~/LM/Fig5/')
outdir <- "./files/"
dir.create(outdir,recursive = T)
outdir2 <- "./plots/"
dir.create(outdir2,recursive = T)

my_theme <- theme(
  panel.background = element_rect(fill = "white"),
  plot.title = element_text(size = 8, face = "bold"),
  axis.title = element_text(size = 8),
  axis.line = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  legend.position = "right",
  legend.background = element_rect(fill = "white"),
  axis.text.y = element_text(size = 6, colour = 'black'),
  axis.text.x = element_blank(),axis.ticks.x = element_blank(),
  axis.title.x =  element_blank(),aspect.ratio = 1)

group_colors <- c("#228B22", "#ffbc14", "firebrick3")
#### 1.Fig5A ----
S1 <- read_csv('./TCR/BLBL0804/filtered_contig_annotations.csv')
S2 <-  read_csv('./TCR/BLJK0804/filtered_contig_annotations.csv')
S3 <- read_csv('./TCR/CNLHL420/filtered_contig_annotations.csv')
S4 <-  read_csv('./TCR/BFYY331/filtered_contig_annotations.csv')
S5 <-  read_csv('./TCR/BLJJ329/filtered_contig_annotations.csv')
S6 <-  read_csv('./TCR/BHLJ1022/filtered_contig_annotations.csv')
S7 <-  read_csv('./TCR/BJZY1103/filtered_contig_annotations.csv')
S8 <-  read_csv('./TCR/BMYY0203/filtered_contig_annotations.csv')
S9 <-  read_csv('./TCR/CFYY331/filtered_contig_annotations.csv')
S10 <-  read_csv('./TCR/CLJJ329/filtered_contig_annotations.csv')
S11 <-  read_csv('./TCR/CGCJ028/filtered_contig_annotations.csv')

combined <- combineTCR(
  contig_list <- list(S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11),
  samples = c("BLBL804" ,"BLJK804" ,"CNLHL420",
              "BFYY331" ,"BLJJ329","BHLJ1022","BJZY1103","BMYY203",
              "CFYY331" ,"CLJJ329", "CGCJ1028"),
  ID = c("Normal","Normal","Normal", 
         "LM", "LM", "LM", "LM", "LM",
         "PE","PE","PE"))
for(i in 1:11){
  combined[[i]]$sample <- factor(combined[[i]]$sample, levels = c("BLBL804" ,"BLJK804" ,"CNLHL420",
                                                                  "BFYY331" ,"BLJJ329","BHLJ1022","BJZY1103","BMYY203",
                                                                  "CFYY331" ,"CLJJ329", "CGCJ1028"))
  combined[[i]]$ID <- factor(combined[[i]]$ID, levels = c("Normal","LM","PE"))
  combined[[i]]$barcode <- gsub("_Normal_","_",combined[[i]]$barcode)
  combined[[i]]$barcode <- gsub("_LM_","_",combined[[i]]$barcode)
  combined[[i]]$barcode <- gsub("_PE_","_",combined[[i]]$barcode)
}

seurat <- qs::qread('~/LM/merge/files/CD4CD8_sub.qs')
DefaultAssay(seurat) <- 'RNA'
seurat$barcode_ori <- rownames(seurat@meta.data)
rownames(seurat@meta.data)  <- paste(seurat$orig.ident, colnames(seurat), sep = '_')
rownames(seurat@meta.data)  <- gsub('.{2}$', '', rownames(seurat@meta.data) )
rownames(seurat@meta.data)  <- gsub('-1_', '-1', rownames(seurat@meta.data) )
qs::qsave(seurat, file.path(outdir, 'seurat_ori.qs'))

seurat1 <- combineExpression(
  combined,
  seurat,
  cloneCall = "gene",
  proportion = FALSE,
  cloneTypes=c(Single=1, Small=5, Medium=20, Large=500)
)
seurat1@meta.data$cloneType <- factor(seurat1@meta.data$cloneType, 
                                      levels = c("Large (20 < X <= 500)", 
                                                 "Medium (5 < X <= 20)", 
                                                 "Small (1 < X <= 5)",
                                                 "Single (0 < X <= 1)", 
                                                 NA))
rownames(seurat1@meta.data)  <- seurat1$barcode_ori
qs::qsave(seurat1, file.path(outdir, 'combined_seurat.qs'))

colorblind_vector <- colorRampPalette(rev(c('#ffb480','#ec6800','#b44c97','#65318e')))
DimPlot(seurat1, group.by = "cloneType") + 
  scale_color_manual(values=colorblind_vector(4), na.value="#d3d6db") + 
  theme(plot.title = element_blank())+ 
  theme(aspect.ratio = 1)
ggsave(file.path(outdir2,'UMAP_cloneType.pdf'), last_plot(), height= 4, width= 8)

#### 2.Fig5B ----
seurat_obj <- qs:qread(paste0(outdir, 'combined_seurat.qs'))

pdf(file = paste0(outdir2, "prop_group.pdf"),width= 5,height= 3)
data <- seurat_obj@meta.data
ggpiestats(data, 'cloneType', 
           results.subtitle = F,
           #factor.levels = c('4', '6', '8'),
           slice.label = 'percentage',
           perc.k = 2,
           direction = 1,
           title = 'ALL')+
  scale_fill_manual(values = type_colors)+theme(plot.title = element_text(vjust = -5,hjust = 0.5))
dev.off()

#### 3.Fig5C ----
seurat <- qs::qread(file.path(outdir, 'combined_seurat.qs'))
occupiedscRepertoire(seurat, x.axis = "group",label = F,proportion = T)+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.text = element_text(colour = 'black',size = 8),
        axis.title = element_text(colour = 'black',size = 10))+
  scale_fill_manual(values = colorblind_vector(5))
ggsave(file.path(outdir2,'bar_pro_occupiedscRepertoire_group.pdf'),
       ggplot2::last_plot(),
       height=4,
       width=7)

#### 4.Fig5D ----
prop <- as.data.frame(prop.table(table(Idents(seurat_obj), seurat_obj$orig.ident), margin = 2) * 100)
colnames(prop) <- c('types','donor','prop')
prop$group <- ''
prop$group[prop$donor%in%c("BLBL804" ,"BLJK804" ,"CNLHL420")] <- 'Normal'
prop$group[prop$donor%in%c("BMYY203","BFYY331"  ,"BLJJ329","BHLJ1022","BJZY1103")] <- 'LM'
prop$group[prop$donor%in%c("CFYY331","CLJJ329"  ,"CGCJ1028")] <- 'PE'
prop$group <- factor(prop$group,levels = c('Normal','LM','PE'))
prop$types <- factor(prop$types,levels = celltypes)

celltypes <- c('Large (20 < X <= 500)','Medium (5 < X <= 20)', 'Small (1 < X <= 5)', 'Single (0 < X <= 1)')
prop_LM_N <- prop|> filter(types %in% 'Single (0 < X <= 1)')
p <- prop_LM_N %>%
  ggplot(aes(x = group, y = prop, fill = group)) +
  geom_boxplot() +geom_point( size=1.5)+
  ggsignif::geom_signif(textsize = 2, comparisons = list(c("Normal", "LM"),c("LM", "PE"))) +
  labs(title = "TCR",y = 'Percentage of unique clones (%)') + theme_classic()+
  theme(aspect.ratio = 1.5,
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  scale_fill_manual(values = group_colors) +
  NoLegend()
ggsave(file.path(outdir2, paste0(case, "Single_group.pdf")), p,width = 3,height = 3)

#### 5.Fig5E ----
seurat <- subset(seurat, group%in%c('Normal','LM'))
occupiedscRepertoire(seurat, x.axis = "ident",label = F,proportion = T)+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        axis.text = element_text(colour = 'black',size = 8),
        axis.title = element_text(colour = 'black',size = 10))+
  scale_fill_manual(values = colorblind_vector(5))
ggsave(file.path(outdir2,'bar_pro_occupiedscRepertoire_type.pdf'), last_plot(), height=4, width=7)

#### 6.Fig5F ----
CD8_GZMB <- subset(seurat, cell_type%in%c('CD8.c05(GZMB+Tem)'))
CD8_GZMB1 <- subset(seurat, cell_type%in%c('CD8.c04(GZMB-Tem)'))
CD4_GZMA <- subset(seurat, cell_type%in%c('CD4.c04(GZMA+Tem)'))

colorblind_vector <- colorRampPalette(rev(c('#ffb480','#ec6800','#b44c97','#65318e')))
p1=occupiedscRepertoire(CD8_GZMB, x.axis = "group",
                        proportion = T,label = F)+
  labs(title = 'CD8.c05(GZMB+Tem)')+
  theme(axis.text = element_text(colour = 'black',size = 8),
        axis.title = element_text(colour = 'black',size = 10),
        plot.title = element_text(colour = 'black',size = 10))+
  scale_fill_manual(values = colorblind_vector(4))+
  NoLegend()
p2=occupiedscRepertoire(CD8_GZMB1, x.axis = "group",
                        proportion = T,label = F)+
  labs(title = 'CD8.c04(GZMB-Tem)')+
  theme(axis.text = element_text(colour = 'black',size = 8),
        axis.title = element_text(colour = 'black',size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(colour = 'black',size = 10))+
  scale_fill_manual(values = colorblind_vector(4))+
  NoLegend()
p3=occupiedscRepertoire(CD4_GZMA, x.axis = "group",
                        proportion = T,label = F)+
  labs(title = 'CD4.c04(GZMA+Tem)')+
  theme(axis.text = element_text(colour = 'black',size = 8),
        axis.title = element_text(colour = 'black',size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(colour = 'black',size = 10))+
  scale_fill_manual(values = colorblind_vector(4))
p4=p1+p2+p3
ggsave(file.path(outdir2,'bar_pro_occupiedscRepertoire_stype.pdf'), last_plot(), height=4, width=8)

#### 7.Fig5G ----
seurat_obj$cloneType2 <- ''
seurat_obj$cloneType2[seurat_obj$cloneType%in%c("Single (0 < X <= 1)")] <- 'unique'
seurat_obj$cloneType2[seurat_obj$cloneType%in%c('Large (20 < X <= 500)','Medium (5 < X <= 20)', 
                                                'Small (1 < X <= 5)')] <- 'multiple'
seurat_obj$cloneType2 <- factor(seurat_obj$cloneType2,levels = c('unique','multiple'))
data <- seurat_obj@meta.data |> filter(cell_type %in% c('CD8.c04(GZMB-Tem)','CD8.c05(GZMB+Tem)')) 
data$type_id <- paste(data$cell_type,data$orig.ident,sep = '_')

df <- as.data.frame(prop.table(table(data$cloneType2, data$type_id), margin = 2) * 100)
split_col <- strsplit(as.character(df$Var2), "_")
df$cell_type <- sapply(split_col, function(x) x[1])
df$orig.ident <- sapply(split_col, function(x) x[2])
df$type <- NULL
df$group <- ''
df$group[df$orig.ident%in%c("BLBL804" ,"BLJK804" ,"CNLHL420")] <- 'Normal'
df$group[df$orig.ident%in%c("BMYY203","BFYY331"  ,"BLJJ329","BHLJ1022","BJZY1103")] <- 'LM'
df$group[df$donor%in%c("CFYY331","CLJJ329"  ,"CGCJ1028")] <- 'PE'
df$group <- factor(df$group,levels = c('Normal','LM','PE'))
colnames(df)<- gsub('Var1','clonetype',colnames(df))

p <- ggboxplot(df, x = "cell_type", y = "Freq",
               color = "clonetype", palette = c("#f26115","#65318e"),add = "jitter")+
  stat_compare_means(aes(group = clonetype), label = "p",method = 'wilcox.test')+
  rotate_x_text(angle = 0,vjust = 0.5,color = "black", size = 8)+
  labs(x='',y="Percentage of clonetypes (%)")+ theme(legend.position = 'top', axis.text.x = element_text(angle = 0))
ggsave(file.path(outdir2, "GZMB_cloneType2_group.pdf"), p,width = 4,height = 4.5)

#### 8.Fig5H ----
seurat_obj0 <- qs:qread(paste0(outdir, 'combined_seurat.qs'))
seurat_obj <- seurat_obj0 %>% subset(.,cell_type == 'CD8.c05(GZMB+Tem)')

s_Gene <- c('GZMB','GNAS','XIST','JUND','DIP2A','CCR7','PDE3B',
            'GZMA','GZMH','GZMB','GNLY','PRF1','GZMK','NKG7','ZNF683','CX3CR1','CCL4','CCL5','IFNG',
            'KLRG1','KLRB1','TIGIT','IKZF2','SELL','IL7R','TBX21','EOMES','CD244',
            'TRBV10-2', 'TRBV5-4', 'TRGV4', 'TRAV26-2',
            'TRAV29DV5','TRBV27','TRBV2','TRBV9','TRAV26-1','TRGV10','TRAV26-1',
            'TRAV38-1','TRAC',
            'TRG-AS1', 'TRGV3','TRAPPC10','TRGC2','TRBC1','TRAV26-2')
s_Gene <- unique(s_Gene)
expr_byCell <- as.data.frame(seurat_obj@assays$RNA@counts)
expr_byCell <- expr_byCell[s_Gene,]
seurat_obj$cell_id <- colnames(seurat_obj)
annotation_col <- as.data.frame(seurat_obj@meta.data[, c('group',"orig.ident",'cloneType')])
rownames(annotation_col) <- rownames(seurat_obj@meta.data)
annotation_colors <- list(orig.ident = c("BLBL804" = "#53A85F", "BLJK804" = "#23452F", "CNLHL420" = "#91D0BE",
                                         'BFYY331'="#F1BB72", 'BLJJ329'="#D6E7A3", 'BHLJ1022'="#c5942e", 'BJZY1103'="#fabf04",
                                         'BMYY203'="#FFFF00", 'CFYY331'="#F3B1A0", 'CLJJ329'="#E95C59", 'CGCJ1028'="#E59CC4"),
                          group=c('Normal'='#228B22', 'LM'='#ffbc14', 'PE'="firebrick3"),
                          cloneType = c('Single (0 < X <= 1)'="#ffb480", 'Small (1 < X <= 5)'="#ec6800",
                                        "Medium (5 < X <= 20)"="#b44c97", "Large (20 < X <= 500)"="#65318e"))
scaled <- scale(t(expr_byCell), center = T, scale = T)
scaled <- t(scaled)
scaled[scaled >= 2] <- 2
scaled[scaled <= -2] <- -2

pdf(file.path(outdir2, "G_heatmap.pdf"), width = 15,height = 10)
pheatmap(scaled,
         annotation_col = annotation_col, 
         annotation_colors = annotation_colors, 
         color = c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
                   colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)),
         treeheight_row = 10, treeheight_col = 10,
         border_color = "white",
         cellwidth = 0.2, cellheight = 5,
         fontsize = 5, angle_col = "45",
         cluster_rows = T, cluster_cols = F,
         show_rownames = T, show_colnames = F,
         gaps_col = c(1521,3195),
         clustering_method = "average")
dev.off()

#### 9.Fig5I ----
quantContig_output <- quantContig(
  combined,
  cloneCall = "gene",
  scale = T,
  exportTable = T)
quantContig_output$group <- c(rep('Normal',3),rep('LM',5),rep('PE',3))
quantContig_output$group <- factor(quantContig_output$group, levels = c('Normal', 'LM','PE'))
quantContig_output %>%
  ggplot(aes(x = group, y = contigs, fill = group)) +
  geom_boxplot() + geom_point(size=1.5) +
  ggsignif::geom_signif(textsize = 2, comparisons = list(c("Normal", "LM"), c("LM", "PE"))) +
  labs(title = "TCR",y = 'Diversity in T cells') + theme_classic()+
  theme(aspect.ratio = 1.5,
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  scale_fill_manual(values = group_colors) +
  NoLegend()
ggsave(file.path(outdir2, "Diversity_group.pdf"), p,width = 3,height = 3)

#### 10.Fig5J ----
s_type <- c('CD4.c01(SELL+Tn)', 'CD4.c02(ISG+Tn)', 'CD8.c01(CCR7+Tn)')

data <- seurat_obj@meta.data
data$type_id <- paste(data$cell_type,data$orig.ident,sep = '_')

type_split <- split(data, data$type_id)
quantContig_celltype <- quantContig(type_split, cloneCall = "gene", scale = T, exportTable = T)
df <- quantContig_celltype
split_col <- strsplit(as.character(df$values), "_")
df$cell_type <- sapply(split_col, function(x) x[1])
df$orig.ident <- sapply(split_col, function(x) x[2])
df$type <- NULL
df$group <- ''
df$group[df$orig.ident%in%c("BLBL804" ,"BLJK804" ,"CNLHL420")] <- 'Normal'
df$group[df$orig.ident%in%c("BMYY203","BFYY331"  ,"BLJJ329","BHLJ1022","BJZY1103")] <- 'LM'
df$group[df$orig.ident%in%c("CFYY331","CLJJ329"  ,"CGCJ1028")] <- 'PE'
df$group <- factor(df$group,levels = c('Normal','LM','PE'))

for(i in 1:length(s_type)){
  df2 <- df[df$cell_type %in% s_type[i],]
  assign(paste0("p", i),
         ggboxplot(df2, x = 'group', y = 'contigs',
                   color = "group", palette = group_colors,add = "jitter") + 
           scale_fill_manual(values = group_colors) +
           labs(x = '', y = 'Diversity in T cells', fill = 'cell_type')+
           theme(axis.text.x = element_text(angle=0,hjust=0.5,vjust=1),
                 axis.text.x.bottom = element_text(size = 10),
                 aspect.ratio = 1, axis.ticks = element_blank())+NoLegend()+ 
           facet_wrap(~ cell_type, ncol = 3)+
           ggsignif::geom_signif(textsize = 3,vjust = 1.4, comparisons = list(c("Normal", "LM"), c("LM", "PE"))))
}
p <- p1+p2+p3
ggsave(file.path(outdir2, paste0(case, "Diversity_s_celltype_group.pdf")), p,width = 10,height = 3.5)



