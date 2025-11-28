library(Seurat)
library(tidyverse)
library(cowplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(qs)
library(tidydr)

setwd('~/LM/Fig1/')
outdir <- "./files/"
dir.create(outdir,recursive = T)
outdir2 <- "./plots/"
dir.create(outdir2,recursive = T)

group_colors <- c( "#228B22", "#ffbc14", "firebrick3")
#### 1.Fig1C/E ----
case = 'named_'
seurat_obj <- qs:qread('~/LM/merge/files/named_seurat.qs')
df = data.frame(clu=names(table(seurat_obj$cell_type)),
                per=sprintf("%1.2f%%", 100*table(seurat_obj$cell_type)/length(seurat_obj$cell_type)))
seurat_obj$per = df[match(seurat_obj$cell_type,df$clu),2]
seurat_obj$type_per = paste0(seurat_obj$cell_type," (",seurat_obj$per,")")
seurat_obj$type_per <- factor(seurat_obj$type_per, 
                              levels = c(
                                'B cells (12.33%)',  'Plasma cells (0.38%)', 
                                'CD4+ T cells (26.43%)',  'CD8+ T cells (21.69%)','MAIT (2.47%)',  'Gamma delta T cells (2.49%)', 'NK cells (10.34%)',
                                'CD14+ Monocytes (15.09%)', 'CD14+CD16+ Monocytes (2.91%)', 'CD16+ Monocytes (1.45%)', 'Dendritic cells (0.95%)',
                                'Neutrophils (0.79%)', 'pDC (0.36%)', 'Platelets (1.63%)',  'Unknown (0.69%)' ))
df1=seurat_obj@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(cluster=seurat_obj$type_per)
df2=sapply( split(df1[, 1:2], df1$cluster), function(x){ apply(x,2,mean) }) %>%
  t() %>% as.data.frame()
df2$cluster=rownames(df2)
p1 <- DimPlot(seurat_obj,pt.size = 2, group.by = 'type_per', 
              raster = TRUE, raster.dpi = c(1024, 1024))+
  theme_void() +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"),
                               ends = 'last', type = "closed")) + 
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs( x= "UMAP 1",y= "UMAP 2", title = '')+
  theme(plot.margin = margin(5.5,15,5.5,5.5),
        panel.grid = element_blank(),aspect.ratio = 1) +
  scale_color_manual(values = type_colors) +
  scale_fill_manual(values = type_colors)
ggsave(paste0(outdir2,'named_UMAP_seurat_cluster_per_lowR.pdf'),
       p1, height=5, width=7)

p2 <- DimPlot(seurat_obj,pt.size = 3, group.by = 'cell_type', split.by = 'group', 
              raster = TRUE, raster.dpi = c(1024, 1024))+
  theme(plot.margin = margin(5.5,15,5.5,5.5),
        panel.grid = element_blank(),aspect.ratio = 1) +
  scale_color_manual(values = type_colors) +
  labs( x= "UMAP 1",y= "UMAP 2", title = '')+
  scale_fill_manual(values = type_colors)+NoLegend()
ggsave(file.path(outdir2,paste0(case,'split_group_lowR.pdf')),
       p2, height=3.5, width=8)

#### 2.Fig1D ----
library(ComplexHeatmap)

merged <- subset(seurat_obj, idents = 'Unknown', invert = T)
features <- c('MS4A1','CD79B','CD79A' ,
              'JCHAIN','MZB1' ,'IGHA1' ,
              'CD4',"CD40LG",
              "CD8B",'CD8A','CCL5',
              'SLC4A10','TRAV1-2','KLRB1',
              'TRDV2', 'TRGV9','TRDC', 
              'GNLY','KLRD1','NKG7',
              'S100A8','CD14','VCAN',
              'IFI30','FCN1','AIF1',
              'CD68','CSF1R','FCGR3A',
              "CD1C","CLEC10A",'HLA-DQA1',
              'CSF3R', 'MNDA', 'FCGR3B',
              "LILRA4","TCF4","CLEC4C",
              "PPBP", "PF4", "GP1BB"
)
gene_cell_exp <- AverageExpression(merged,
                                   features = features,
                                   group.by = 'cell_type',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'cell_type'
df$cell_type <- factor(df$cell_type, 
                       levels = c('B cells',  'Plasma cells',  'CD4+ T cells', 'CD8+ T cells', 'MAIT', 'Gamma delta T cells', 'NK cells',
                                  'CD14+ Monocytes', 'CD14+CD16+ Monocytes', 'CD16+ Monocytes', 'Dendritic cells', 
                                  'Neutrophils', 'pDC',  'Platelets'))
top_anno = HeatmapAnnotation(df = df,
                             border = F, annotation_height = c(0.05),
                             show_annotation_name = F,show_legend = F,
                             gp = gpar(col = "white"),
                             col = list(cell_type = c('B cells'= '#53A85F', 
                                                      'Plasma cells'= '#585658',
                                                      'CD4+ T cells'= '#E59CC4',
                                                      'CD8+ T cells'= '#F1BB72',
                                                      "MAIT"= '#E95C59',
                                                      "Gamma delta T cells"= '#F3B1A0',
                                                      'NK cells'= '#8C549C',
                                                      'CD14+ Monocytes'= '#f46f20',
                                                      'CD14+CD16+ Monocytes'= "#6a73cf",
                                                      'CD16+ Monocytes'= '#0eb0c8',
                                                      "Dendritic cells"= '#476D87',
                                                      'Neutrophils'='#91D0BE',
                                                      'Platelets'= '#AA9A59',
                                                      'pDC'= '#E39A35')))
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
marker_exp[marker_exp > 1.5] <- 1.5
marker_exp[marker_exp < -1.5] <- -1.5
colnames(marker_exp) <- factor(colnames(marker_exp), levels = c(
  'B cells',  'Plasma cells',  'CD4+ T cells', 'CD8+ T cells', 'MAIT', 'Gamma delta T cells', 'NK cells',
  'CD14+ Monocytes', 'CD14+CD16+ Monocytes', 'CD16+ Monocytes', 'Dendritic cells', 
  'Neutrophils', 'pDC',  'Platelets'
))
col_cluster <- setNames(c(rep('#53A85F',3),rep('#585658',3), 
                          rep('#E59CC4',2),
                          rep('#F1BB72',3), rep('#E95C59',3),
                          rep('#F3B1A0',3),rep('#8C549C',3),
                          rep('#f46f20',3),rep('#6a73cf',3),
                          rep('#0eb0c8',3),
                          rep('#476D87',3),rep('#91D0BE',3),
                          rep('#AA9A59',3),rep('#E39A35',3)),
                        rownames(marker_exp))
row_info = rowAnnotation(foo = anno_text(rownames(marker_exp), 
                                         location = 0, 
                                         just = "left",
                                         gp = gpar(fill = col_cluster, fontface="italic", col = "white"),
                                         width = max_text_width(rownames(marker_exp))*1.1))
s_pathway <- c("B cell activation","B cell receptor signaling pathway","B cell differentiation","regulation of hemopoiesis",
               "mitotic cell cycle phase transition","response to virus","T-helper cell differentiation","T cell receptor signaling pathway",
               "T cell differentiation","type II interferon production","regulation of T cell activation","regulation of immune effector process",                   
               "regulation of cell killing","leukocyte mediated cytotoxicity","gamma-delta T cell activation","gamma-delta T cell differentiation",                    
               "leukocyte mediated cytotoxicity","natural killer cell mediated immunity","natural killer cell mediated cytotoxicity",
               "natural killer cell activation","phagocytosis","positive regulation of cytokine production","mononuclear cell migration",
               "myeloid leukocyte activation","response to type II interferon","myeloid leukocyte differentiation",
               "leukocyte activation involved in inflammatory response","leukocyte activation involved in immune response",
               "regulation of inflammatory response","leukocyte cell-cell adhesion","antigen processing and presentation of exogenous antigen",
               "MHC class II protein complex assembly","leukocyte chemotaxis","positive regulation of phagocytosis",
               "neutrophil chemotaxis","dendrite development","dendrite morphogenesis","interferon-alpha production",
               "hemostasis","blood coagulation","platelet aggregation")
row_info2 = rowAnnotation(foo = anno_text(s_pathway,
                                          location = 0.01,
                                          just = "left",
                                          gp = gpar(fill = col_cluster, col = "white"),
                                          width = max_text_width(s_pathway)*1.1))
pdf(file.path(outdir2,paste0(case,'markers_heatmap.pdf')), width=10, height=7)
p <- Heatmap(marker_exp,
             cluster_rows = F,cluster_columns = F,
             show_column_names = F,show_row_names = T,
             column_title = NULL,
             col = c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
                     colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)),
             heatmap_legend_param = list(at = c(-1.5, 0, 1.5),title = "z scale"),
             border = 'white',
             rect_gp = gpar(col = "white", lwd = 1),
             row_names_gp = gpar(fontsize = 10),
             column_names_gp = gpar(fontsize = 10),
             top_annotation = top_anno)+
  row_info+row_info2
dev.off()

#### 3.Fig1F ----
seurat_obj <- qs:qread('~/LM/merge/files/named_seurat.qs')
seurat_obj <- subset(seurat_obj, idents = 'Unknown', invert = T)
qs::qsave(seurat_obj, '~/LM/merge/files/named_seurat_no_unknown.qs')
case <- 'propor_all_'

prop.table(table(Idents(seurat_obj), seurat_obj$orig.ident), margin = 2) * 100
prop <- as.data.frame(prop.table(table(Idents(seurat_obj), seurat_obj$orig.ident), margin = 2) * 100)
colnames(prop) <- c('types','donor','prop')
prop$group <- ''
prop$group[prop$donor%in%c("BLBL804" ,"BLJK804" ,"CNLHL420")] <- 'Normal'
prop$group[prop$donor%in%c("CFYY331" ,"CGCJ1028" ,"CLJJ329")] <- 'PE'
prop$group[prop$donor%in%c("BMYY203","BFYY331"  ,"BLJJ329","BHLJ1022","BJZY1103")] <- 'LM'
prop$group <- factor(prop$group,levels = c('Normal','LM','PE'))
prop$types <- factor(prop$types,levels = celltypes)
write.table(prop, file.path(outdir, paste0(case, "box_group_type.csv")),
            quote = F, sep = ",", row.names = F)

#### plot
library(rstatix)
prop <- read.csv(file.path(outdir, paste0(case, "box_group_type.csv")))
##
prop1 <- prop[prop$types %in% c('CD4+ T cells', 'CD8+ T cells'),]
stat_t <- t_test(group_by(prop1, types), prop~group)
stat_t <- add_significance(stat_t, 'p')
stat_t.test <-  add_xy_position(stat_t, dodge = 3)
p_box <- ggboxplot(prop1, x = 'group', y = 'prop',
                   facet.by = 'types', fill = 'group',
                   color = 'black',repel = F,
                   width = 0.5, size = 0.3, legend = 'right') +
  scale_fill_manual(values = group_colors) +
  labs(x = '', y = 'Percentage of cells (%)', fill = 'types')+
  theme(axis.text.x = element_text(angle=0,hjust=0.5,vjust=1),
        axis.text.x.bottom = element_text(size = 10),
        aspect.ratio = 1, axis.ticks = element_blank())+NoLegend()+ 
  facet_wrap(~ types, ncol = 2)+
  stat_pvalue_manual(stat_t.test, label = 'p', tip.length = 0.01,
                     bracket.nudge.y =-3, vjust = 1.5,
                     color = 'black',label.size = 2.5,label.y = c(29, 35, 40),
                     bracket.size = 0.2)
ggsave(file.path(outdir2, paste0(case, "face_group_type_group1.pdf")),
       ggplot2::last_plot(),width = 5,height = 2.5)
##
prop1 <- prop[prop$types %in% c('B cells', 'CD16+ Monocytes'),]
stat_t <- t_test(group_by(prop1, types), prop~group)
stat_t <- add_significance(stat_t, 'p')
stat_t.test <-  add_xy_position(stat_t, dodge = 3)
p_box <- ggboxplot(prop1, x = 'group', y = 'prop',
                   facet.by = 'types', fill = 'group',
                   color = 'black',repel = F,
                   width = 0.5, size = 0.3, legend = 'right') +
  scale_fill_manual(values = group_colors) +
  labs(x = '', y = 'Percentage of cells (%)', fill = 'types')+
  theme(axis.text.x = element_text(angle=0,hjust=0.5,vjust=1),
        axis.text.x.bottom = element_text(size = 10),
        aspect.ratio = 1, axis.ticks = element_blank())+NoLegend()+ 
  facet_wrap(~ types, ncol = 2)+
  stat_pvalue_manual(stat_t.test, label = 'p', tip.length = 0.01,
                     bracket.nudge.y =-3, vjust = 1.5,
                     color = 'black',label.size = 2.5,label.y = c(29, 35, 40),
                     bracket.size = 0.2)
ggsave(file.path(outdir2, paste0(case, "face_group_type_group2.pdf")),
       ggplot2::last_plot(),width = 5,height = 2.5)
##
prop1 <- prop[prop$types %in% c('CD14+ Monocytes','CD14+CD16+ Monocytes'),]
stat_t <- t_test(group_by(prop1, types), prop~group)
stat_t <- add_significance(stat_t, 'p')
stat_t.test <-  add_xy_position(stat_t, dodge = 3)
p_box <- ggboxplot(prop1, x = 'group', y = 'prop',
                   facet.by = 'types', fill = 'group',
                   color = 'black',repel = F,
                   width = 0.5, size = 0.3, legend = 'right') +
  scale_fill_manual(values = group_colors) +
  labs(x = '', y = 'Percentage of cells (%)', fill = 'types')+
  theme(axis.text.x = element_text(angle=0,hjust=0.5,vjust=1),
        axis.text.x.bottom = element_text(size = 10),
        aspect.ratio = 1, axis.ticks = element_blank())+NoLegend()+ 
  facet_wrap(~ types, ncol = 2)+
  stat_pvalue_manual(stat_t.test, label = 'p', tip.length = 0.01,
                     bracket.nudge.y =-3, vjust = 1.5,
                     color = 'black',label.size = 2.5,label.y = c(29, 35, 40),
                     bracket.size = 0.2)
ggsave(file.path(outdir2, paste0(case, "face_group_type_group3.pdf")),
       ggplot2::last_plot(),width = 5,height = 2.5)