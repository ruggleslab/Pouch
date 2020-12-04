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

#
#
#
# Heatmap of top genes
#
#
#
#
marker_list = read.table("OBJECTS/T_cells/clusters-MinorPopulations-clust12/markers-global/markers.clust12.wilcox.all.csv", T, ',')
clusts = levels(s_obj@meta.data$MinorPopulations)

good_genes <- data.frame(gene=NA, ct = NA)
for(i in 1:length(clusts)){
  type = clusts[i]
  gener = subset(marker_list, cluster == clusts[i])
  gener <- gener[order(gener$avg_logFC, gener$p_val_adj, decreasing=T),]
  
  
  genes = as.character(gener$gene[1:10])
  genes2 = data.frame(gene=genes, ct = clusts[i])
  print(length(genes))
  good_genes <- rbind(good_genes, genes2)
}
#good_genes <- c(good_genes, "CD8A", "CD4", "CD3D", "CD79A", "FOS", "FOXP3", "IL17A")
good_genes <- good_genes[-1,]
good_good <- unique(good_genes$gene)

good_good <- intersect(good_good, rownames(s_obj))
print(good_good)

good_genes2 <- subset(good_genes, good_genes$gene %in% good_good)
good_genes2 <- subset(good_genes2, !is.na(gene))
good_genes2$ct <- factor(good_genes2$ct, levels = levels(s_obj@meta.data$MinorPopulations))
good_genes2 <- good_genes2[order(good_genes2$ct),]

#tab2write
fig3a_tab = subset(marker_list, gene %in% good_good)
write.table(fig3a_tab, "FIGURES/FIGURE3/Figure3A_Top_markers.txt", sep='\t', row.names=F, quote=F)
#

counts <- data.frame(s_obj@assays$integrated@scale.data)[good_genes2$gene,]
tcounts = t(counts)

rownames(s_obj@meta.data) <- make.names(colnames(s_obj))
clusterofchoice = "MinorPopulations"
cluster2 = s_obj@meta.data[colnames(counts),clusterofchoice]
cluster2 = factor(cluster2, levels = levels(s_obj@meta.data$MinorPopulations))

#scounts <- aggregate(tcounts, by=list(cluster2), 'median')
#rownames(scounts) <- scounts[,1]
#scounts <- scounts[,-1]

#scounts[is.na(scounts)] <- 0

ha <- rowAnnotation(df = data.frame(Cluster=cluster2),
                    col=list(Cluster = big_colorer))

namer <- paste0("FIGURES/FIGURE3/Figure3A_heatmap.png")
png(namer, height = 10, width = 15, units = "in", res=150)
Heatmap(data.matrix(tcounts), show_row_names = F, show_column_names = T,
        #top_annotation = ha,
        heatmap_legend_param = list(title = "Scaled Value"),
        cluster_rows = T, cluster_columns = F, #row_names_side = 'left',
        column_names_gp = gpar(fontsize=12),
        row_title_gp = gpar(fontsize = 10),
        row_names_max_width = unit(10,'cm'),
        use_raster = T,
        cluster_row_slices=F,
        split = cluster2,
        left_annotation = ha,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
dev.off()
