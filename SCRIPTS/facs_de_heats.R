library(ggplot2)
library(reshape)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(RColorBrewer)
library(ggsci)
library(reshape)
library(Seurat)
library(ggrepel)
library(DESeq2)
library(matrixStats)

colorer <- c("Normal_pouch"="dodgerblue2","Pouchitis"="red3", "UC_colon"="forestgreen", "UC_inflamed"="darkorange1", "other"="grey85")
colors_clusters = c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))

#
#
#
#### heatmap for each de results from facs
des = c("OBJECTS/Myeloid_cells/clusters-MinorPopulations-clust5/diff-expression-orig.ident/SOX4/clusters-Facs_pop-clust1/diff-expression-orig.ident/de.clust1.orig.ident.wilcox.all.csv",
        "OBJECTS/Myeloid_cells/clusters-MinorPopulations-clust5/diff-expression-orig.ident/IL1B/clusters-Facs_pop-clust1/diff-expression-orig.ident/de.clust1.orig.ident.wilcox.all.csv",
        "OBJECTS/Myeloid_cells/clusters-MinorPopulations-clust5/diff-expression-orig.ident/APOE/clusters-Facs_pop-clust1/diff-expression-orig.ident/de.clust1.orig.ident.wilcox.all.csv",
        "OBJECTS/T_cells/clusters-MinorPopulations-clust12/diff-expression-orig.ident/FOXP3/clusters-Facs_pop-clust1/diff-expression-orig.ident/de.clust1.orig.ident.wilcox.all.csv",
        "OBJECTS/T_cells/clusters-MinorPopulations-clust12/diff-expression-orig.ident/CD8A/clusters-Facs_pop-clust1/diff-expression-orig.ident/de.clust1.orig.ident.wilcox.all.csv")

cutoffs=c(0.75,0.75,0.75,0.5,0.5)

for(de in 1:length(des)){

  de_table = read.table(des[de], T, sep=',')
  
  cells = as.character(unique(de_table$cluster))
  for(i in 1:length(cells)){
    cell_table=subset(de_table, cluster==cells[i])
    cell_table_de=subset(cell_table, abs(avg_logFC)>cutoffs[de] & p_val_adj<0.05)
    good_genes = as.character(unique(cell_table_de$gene))
    cell_table_filt = subset(cell_table, gene %in% good_genes)
    
    if(nrow(cell_table_filt)>0){
      de_mat = cast(cell_table_filt, gene ~ group1+group2, value='avg_logFC')
      rownames(de_mat) = de_mat[,1]
      de_mat=data.matrix(de_mat[,-1])
      
      de_mat[,1]=de_mat[,1]*-1
      de_mat[,3]=de_mat[,3]*-1
      colnames(de_mat) = c("Pouchitis vs. Normal Pouch", "UC Inflamed vs. Normal Pouch",
                           "Pouchitis vs. UC Inflamed")
      de_mat[is.na(de_mat)]<-0
      if(i==3){row_size=4} else{row_size=8}
      
      namer <- gsub("de.clust1.orig.ident.wilcox.all.csv", paste0(cells[i], "_FC_heatmap.pdf"), des[de])
      pdf(namer, height = 10, width = 5)
      print(Heatmap(de_mat, show_row_names = T, show_column_names = T,
                    #top_annotation = ha,
                    heatmap_legend_param = list(title = "Fold Change"),
                    cluster_rows = T, cluster_columns = F, #row_names_side = 'left',
                    column_names_gp = gpar(fontsize=12),
                    row_names_gp = gpar(fontsize=row_size),
                    row_title_gp = gpar(fontsize = 10),
                    #row_names_max_width = unit(10,'cm'),
                    use_raster = T,
                    #cluster_row_slices=F,
                    #split = cluster2,
                    #left_annotation = ha,
                    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
      )
      dev.off()
    }
    
  }
}
