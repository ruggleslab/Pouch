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




#
# All cells UMAP


s_obj <- readRDS("OBJECTS/seurat_obj.rds")

tcells <- subset(s_obj, subset = snn_res.0.4 %in% c("C02", "C03", "C05", "C06", "C07", "C08"))
bcells_gc_f <- subset(s_obj, subset = snn_res.0.4 %in% c("C01", "C10"))
bcells_cycling <- subset(s_obj, subset = snn_res.0.4 %in% c("C11"))
bcells_plasma <- subset(s_obj, subset = snn_res.0.4 %in% c("C04", "C13"))
mcells_monomac <- subset(s_obj, subset = snn_res.0.4 %in% c("C09", "C14"))
mcells_mast <- subset(s_obj, subset = snn_res.0.4 %in% c("C12"))

objs = list(tcells,	
            bcells_gc_f,	
            bcells_cycling ,
            bcells_plasma ,
            mcells_monomac,
            mcells_mast)
names(objs) = c("tcells","bcells_gc_f","bcells_cycling","bcells_plasma",
                "mcells_monomac","mcells_mast")
cells_meta <- data.frame(cells_id = NA, cells_state = NA)
for(i in 1:length(objs)){
  cells_id = colnames(objs[[i]])
  cells_state = names(objs)[i]

  cells_add <- data.frame(cells_id, cells_state)
  cells_meta <- rbind(cells_meta, cells_add)
}
cells_meta <- cells_meta[-1,]

rownames(cells_meta) <- cells_meta$cells_id
cells_meta2 <- cells_meta[colnames(s_obj),]

s_obj@meta.data$MajorPopulations <- cells_meta2$cells_state
s_obj@meta.data$MajorPopulations = factor(s_obj@meta.data$MajorPopulations, levels = c("tcells","bcells_gc_f","bcells_cycling","bcells_plasma",
                                                    "mcells_monomac","mcells_mast", "mcells_pdcs"))

if("MajorPopulations" %in% colnames(s_obj@meta.data)){
} else{
  saveRDS(s_obj, "OBJECTS/seurat_obj.rds")
}

#
#### plot cell states

Idents(s_obj) = s_obj@meta.data$MajorPopulations

colorer <- c("tcells"="forestgreen", 
             "bcells_gc_f"="darkred", "bcells_cycling"="orangered", "bcells_plasma"="lightcoral", 
             "mcells_monomac"="navy", "mcells_mast"="deepskyblue4", "mcells_pdcs"="dodgerblue2")

plot_umap =
  DimPlot(s_obj, reduction = "umap", cells = sample(colnames(s_obj)), pt.size=2, cols.use=colorer) +
  theme(aspect.ratio = 1) + scale_color_manual(values=colorer)

ggsave(paste0("FIGURES/FIGURE1/Figure1B.umap.png"), plot = plot_umap, width = 12, height = 12, units = "in")
ggsave(paste0("FIGURES/FIGURE1/Figure1B.umap.pdf"), plot = plot_umap, width = 12, height = 12, units = "in")



######
### plot different donors

Idents(s_obj) = s_obj@meta.data$orig.ident

colors_clusters = c(pal_d3("category20b")(20), pal_igv("default")(51))


plot_umap =
  DimPlot(s_obj, reduction = "umap", cells = sample(colnames(s_obj)), pt.size=2, cols.use=colors_clusters) +
  theme(aspect.ratio = 1) + scale_color_manual(values=colors_clusters)

ggsave(paste0("FIGURES/SF1/FigureSF1A.umap.png"), plot = plot_umap, width = 12, height = 12, units = "in")
ggsave(paste0("FIGURES/SF1/FigureSF1A.umap.pdf"), plot = plot_umap, width = 12, height = 12, units = "in")


#
#
### plot disease
colorer <- c("Normal_pouch"="dodgerblue2","Pouchitis"="red3", "UC_colon"="forestgreen", "UC_inflamed"="darkorange1")

Idents(s_obj) = s_obj@meta.data$pouch.status

plot_umap =
  DimPlot(s_obj, reduction = "umap", cells = colnames(s_obj), pt.size=2, cols.use=colorer) +
  theme(aspect.ratio = 1) + scale_color_manual(values=colorer)

ggsave(paste0("FIGURES/SF1/FigureSF1B.umap.png"), plot = plot_umap+facet_wrap(~s_obj@meta.data$pouch.status, nrow=1), width = 36, height = 12, units = "in")
ggsave(paste0("FIGURES/SF1/FigureSF1B.umap.pdf"), plot = plot_umap+facet_wrap(~s_obj@meta.data$pouch.status, nrow=1), width = 36, height = 12, units = "in")
