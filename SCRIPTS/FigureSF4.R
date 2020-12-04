library(Seurat)
library(ggplot2)
library(matrixStats)
library(gridExtra)
library(RColorBrewer)
library(ggsci)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(matrixStats)

args = commandArgs(trailingOnly=TRUE)

colorer <- c("Normal_pouch"="dodgerblue2","Pouchitis"="red3", "UC_colon"="forestgreen", "UC_inflamed"="darkorange1")
colors_clusters = c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))

big_colorer <- c(act_CD8_NR4A2_hi=colors_clusters[1], 
                 act_CD8_NR4A2_lo=colors_clusters[2], 
                 activated_CD4=colors_clusters[3], 
                 CD4_fos_hi=colors_clusters[4], 
                 CD4_fos_lo=colors_clusters[5],
                 Foxp3_cells=colors_clusters[6], 
                 memory_CD4=colors_clusters[7], 
                 memory_CD8=colors_clusters[8], 
                 naive_CD8=colors_clusters[9], 
                 NK_cells=colors_clusters[10],
                 Tfh_cells=colors_clusters[11], 
                 Th17_cells=colors_clusters[12])

level_fix = c("act_CD8_NR4A2_hi", 
              "act_CD8_NR4A2_lo", 
              "activated_CD4", 
              "CD4_fos_hi", 
              "CD4_fos_lo",
              "Foxp3_cells", 
              "memory_CD4", 
              "memory_CD8", 
              "naive_CD8", 
              "NK_cells",
              "Tfh_cells", 
              "Th17_cells")

s_obj = readRDS("OBJECTS/T_cells/seurat_obj.rds")

s_obj@meta.data$MinorPopulations = factor(s_obj@meta.data$MinorPopulations, levels = level_fix)
Idents(s_obj) = s_obj@meta.data$MinorPopulations

colorer <- colors_clusters
plot_umap =
  DimPlot(s_obj, reduction = "umap", cells = sample(colnames(s_obj)), pt.size=1, cols.use=big_colorer) +
  theme(aspect.ratio = 1) + scale_color_manual(values=big_colorer)

ggsave(paste0("FIGURES/SF4/FigureSF4a.umap.png"), plot = plot_umap, width = 12, height = 10, units = "in")
ggsave(paste0("FIGURES/SF4/FigureSF4a.umap.pdf"), plot = plot_umap, width = 12, height = 10, units = "in")



#
#
#
#
#### contour of 3 different states
UMAP_1 = s_obj@reductions$umap@cell.embeddings[,1]
UMAP_2 = s_obj@reductions$umap@cell.embeddings[,2]
PouchStatus = s_obj@meta.data$pouch.status

plotter = data.frame(UMAP_1, UMAP_2, PouchStatus)


contourg=ggplot(plotter, aes(UMAP_1, UMAP_2, fill=PouchStatus, color=PouchStatus)) + 
  scale_fill_manual(values=colorer) +
  scale_color_manual(values=colorer) +
  geom_point(alpha=0.15, size=0.25) +
  scale_x_continuous(limits=c(-11,11.5)) +
  scale_y_continuous(limits=c(-7,7)) +
  geom_density_2d(contour=TRUE, n=100) + theme_bw() + facet_wrap(~PouchStatus, nrow=1) +
  theme_void()

  
ggsave(paste0("FIGURES/SF4/FigureSF4_contour.umap.png"), plot = contourg, width = 11, height = 3, units = "in")
ggsave(paste0("FIGURES/SF4/FigureSF4_contour.umap.pdf"), plot = contourg, width = 11, height = 3, units = "in")
  