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
#
#

s_obj <- readRDS("OBJECTS/Myeloid_cells/seurat_obj.rds")
s_obj@meta.data$MinorPopulations2 <- factor(s_obj@meta.data$MinorPopulations2, 
                                           levels = c("mono_mac1", "mono_mac2", "mono_mac3", "mast_cells",
                                                      "DC", "pdcs"))
Idents(s_obj) = s_obj@meta.data$MinorPopulations2
colorer <- c("mast_cells"="navy", 
             "pdcs"="pink", "DC"="darkorange1",
             "mono_mac1"="red3", "mono_mac2"="purple", "mono_mac3"="forestgreen")

genes = c("TREM1", "CXCL10")
gradient_colors = c("gray85", "red2")
plot_umap =
  FeaturePlot(s_obj, features=genes, reduction = "umap", cells = sample(colnames(s_obj)), cols = gradient_colors, ncol = 2, min.cutoff=0)

ggsave("FIGURES/SF3/FigureSF3_TREM1.pdf", plot = plot_umap, width = 10, height = 5, units = "in")
ggsave("FIGURES/SF3/FigureSF3_TREM1.png", plot = plot_umap, width = 10, height = 5, units = "in")

#
# and heatmap of top genes
#
#
#
#
marker_list = read.table("OBJECTS/Myeloid_cells/clusters-MinorPopulations2-clust6/markers-global/markers.clust6.wilcox.all.csv", T, ',')
clusts = levels(s_obj@meta.data$MinorPopulations2)

good_genes <- data.frame(gene=NA, ct = NA)
i=1
type = clusts[i]
gener = subset(marker_list, cluster == clusts[i])
gener <- gener[order(gener$avg_logFC, gener$p_val_adj, decreasing=T),]
gener = subset(gener, avg_logFC>1.2)

genes = as.character(gener$gene)
genes2 = data.frame(gene=genes, ct = clusts[i])
print(length(genes))
good_genes <- rbind(good_genes, genes2)

#good_genes <- c(good_genes, "CD8A", "CD4", "CD3D", "CD79A", "FOS", "FOXP3", "IL17A")
good_genes <- good_genes[-1,]
good_good <- unique(good_genes$gene)

good_good <- intersect(good_good, rownames(s_obj))
print(good_good)

good_genes2 <- subset(good_genes, good_genes$gene %in% good_good)
good_genes2 <- subset(good_genes2, !is.na(gene))
good_genes2$ct <- factor(good_genes2$ct, levels = levels(s_obj@meta.data$MinorPopulations2))
good_genes2 <- good_genes2[order(good_genes2$ct),]

#tab2write
figSF3_tab = subset(marker_list, gene %in% good_good)
write.table(figSF3_tab, "FIGURES/SF3/FigureSF3_Top_markers.txt", sep='\t', row.names=F, quote=F)
#

counts <- data.frame(s_obj@assays$integrated@scale.data)[good_genes2$gene,]
#monomac1 = gsub(":", "\\.", rownames(subset(s_obj@meta.data, MinorPopulations2 == "mono_mac1")))
#monomac1 = gsub("-", "\\.", monomac1)

#monomac2 = gsub(":", "\\.", rownames(subset(s_obj@meta.data, MinorPopulations2 == "mono_mac2")))
#monomac2 = gsub("-", "\\.", monomac2)

#monomac3 = gsub(":", "\\.", rownames(subset(s_obj@meta.data, MinorPopulations2 == "mono_mac3")))
#monomac3 = gsub("-", "\\.", monomac3)

#monomac=c(monomac1, monomac2, monomac3)

#counts = counts[,monomac]

#rownames(s_obj@meta.data) <- make.names(colnames(s_obj))
#clusterofchoice = "MinorPopulations2"
cluster2 = s_obj@meta.data$MinorPopulations2
#cluster2 = factor(cluster2, levels = c("mono_mac1", "mono_mac2", "mono_mac3"))

#scounts <- aggregate(tcounts, by=list(cluster2), 'median')
#rownames(scounts) <- scounts[,1]
#scounts <- scounts[,-1]

#scounts[is.na(scounts)] <- 0

#counts=t(scale(t(counts)))

ha <- columnAnnotation(df = data.frame(Cluster=cluster2),
                    col=list(Cluster = c("mast_cells"="navy", 
                                         "pdcs"="pink", "DC"='darkorange1',
                                         "mono_mac1"="red3", "mono_mac2"="purple", "mono_mac3"="forestgreen")))

namer <- paste0("FIGURES/SF3/FigureSF3_heatmap.pdf")
pdf(namer, height = 8, width = 8)
Heatmap(data.matrix(counts), show_row_names = T, show_column_names = F,
        top_annotation = ha,
        heatmap_legend_param = list(title = "Log2 Expression\nLevel"),
        cluster_rows = T, cluster_columns = T, #row_names_side = 'left',
        row_names_gp = gpar(fontsize=10),
        row_title_gp = gpar(fontsize = 10),
        row_names_max_width = unit(10,'cm'),
        use_raster = T,
        cluster_column_slices=F,
        column_split = cluster2,
        #split = cluster2,
        #left_annotation = ha,
        col = colorRamp2(c(-2,0,2), c("blue", "white", "red")))
dev.off()

