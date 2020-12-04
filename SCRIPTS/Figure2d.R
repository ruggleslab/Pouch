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


s_obj = readRDS("OBJECTS/seurat_obj.rds")
s_obj2 = readRDS("OBJECTS/Myeloid_cells/seurat_obj.rds")
meta <- s_obj@meta.data

#### mono_mac1

plotter <- s_obj@meta.data
plotter$SOX4 <- s_obj@assays$integrated@data["SOX4",]
plotter$MAFA <- s_obj@assays$integrated@data["MAFA",]

plotter <- subset(plotter, MajorPopulations == "mcells_monomac")

plotter$pouch.status2 <- plotter$pouch.status
plotter$pouch.status2[plotter$SOX4 < 0.25] <- "other"
plotter$pouch.status2[plotter$MAFA < 0.25] <- "other"

ggSox = ggplot(plotter, aes(SOX4, MAFA, color=pouch.status2)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=colorer) +
  theme_bw() + facet_wrap(~pouch.status, nrow=1) +
  geom_vline(xintercept =0.25) + geom_hline(yintercept = 0.25)


plotter_int <- subset(plotter, SOX4 > 0.25 & MAFA > 0.25)
samples <- unique(plotter_int$orig.ident)
resser <- data.frame(sampleID = NA, pouchStatus = NA, perc = NA)
for(j in 1:length(samples)){
  curr <- subset(plotter_int, orig.ident == samples[j])
  perc = nrow(curr)/as.numeric(as.character(curr$cell_total[1]))
  adder = data.frame(sampleID = samples[j],
                     pouchStatus = curr$pouch.status[1],
                     perc = perc)
  resser <- rbind(resser, adder)
}
resser=resser[-1,]
foxresser = resser

ressy=pairwise.wilcox.test(x=foxresser$perc, g=foxresser$pouchStatus, p.adjust.method="fdr")$p.value
pvals=as.numeric(ressy)
comps=paste0(rep(rownames(ressy), ncol(ressy)), ";", rep(colnames(ressy), each=nrow(ressy)))
ressy1=data.frame(Comparison=comps, pvals=pvals)
ressy1$cluster="mono_mac1"

bbSox = ggplot(resser, aes(pouchStatus, perc*100, color=pouchStatus)) +
  geom_boxplot(outlier.shape=NA) +
  ggtitle("SOX4+ MAFA+ Cell Percentages") +
  ylab("Cell Percentages") +
  scale_color_manual(values = colorer) +
  geom_jitter(width=0.2) + theme_bw() +
  theme(axis.text.x=element_blank())

###
celler1 = rownames(plotter_int)
s_obj_1 = s_obj2[,celler1]
Idents(s_obj_1) = s_obj_1@meta.data$pouch.status
group_combinations = combn(levels(s_obj_1), m = 2, simplify = TRUE) 
s_obj_1@meta.data$Facs_pop = "SOX4"
saveRDS(s_obj_1, "OBJECTS/Myeloid_cells/clusters-MinorPopulations-clust5/diff-expression-orig.ident/SOX4/seurat_obj.rds")

# res = data.frame(p_val=NA, avg_logFC=NA,pct.1=NA,pct.2=NA,p_val_adj=NA, group1=NA, group2=NA, gene=NA)
# for(j in 1:ncol(group_combinations)){
#   g1 = group_combinations[1, j]
#   g2 = group_combinations[2, j]
#   
#   de_genes1 = FindMarkers(s_obj_1, ident.1 = g1, ident.2 = g2, assay = "RNA",
#                        test.use = 'wilcox', logfc.threshold = log(0), 
#                        min.pct = 0.1, only.pos = FALSE,
#                        print.bar = FALSE)
#   de_genes1$group1=g1
#   de_genes1$group2=g2
#   de_genes1$gene=rownames(de_genes1)
#   res = rbind(res, de_genes1)
# }
# res=res[-1,]
# 
# 
# cell_table_de=subset(res, abs(avg_logFC)>0.75 & p_val_adj<0.05)
# good_genes = as.character(unique(cell_table_de$gene))
# cell_table_filt = subset(res, gene %in% good_genes)
#   
# if(nrow(cell_table_filt)>0){
#   de_mat = cast(cell_table_filt, gene ~ group1+group2, value='avg_logFC')
#   rownames(de_mat) = de_mat[,1]
#   de_mat=data.matrix(de_mat[,-1])
#   
#   de_mat[,1]=de_mat[,1]*-1
#   de_mat[,3]=de_mat[,3]*-1
#   colnames(de_mat) = c("Pouchitis vs. Normal Pouch", "UC Inflamed vs. Normal Pouch",
#                        "Pouchitis vs. UC Inflamed")
#   de_mat[is.na(de_mat)]<-0
#   if(nrow(de_mat)>30){row_size=4} else{row_size=8}
#   
#   namer <- paste0("OBJECTS/Myeloid_cells/clusters-MinorPopulations-clust5/diff-expression-orig.ident/SOX4_MAFA_FC_heatmap.pdf")
#   pdf(namer, height = 10, width = 5)
#   print(Heatmap(de_mat, show_row_names = T, show_column_names = T,
#                 #top_annotation = ha,
#                 heatmap_legend_param = list(title = "Fold Change"),
#                 cluster_rows = T, cluster_columns = F, #row_names_side = 'left',
#                 column_names_gp = gpar(fontsize=12),
#                 row_names_gp = gpar(fontsize=row_size),
#                 row_title_gp = gpar(fontsize = 10),
#                 #row_names_max_width = unit(10,'cm'),
#                 use_raster = T,
#                 #cluster_row_slices=F,
#                 #split = cluster2,
#                 #left_annotation = ha,
#                 col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
#   )
#   dev.off()
# }

  


#
#
#
##deseq on this pop?
#counts = data.matrix(s_obj@assays$RNA@counts[,rownames(plotter_int)])

#counts = counts[rowSums(counts)>10,]
#counts = counts[rowVars(counts)>0.1,]

#library("BiocParallel")
#param <- BatchtoolsParam(workers=10, cluster="slurm", template="SCRIPTS/slurm.tmpl")

#dds <- DESeqDataSetFromMatrix(counts, data.frame(pouch.status=plotter_int$pouch.status), ~ pouch.status)
#dds <- DESeq(dds, parallel=T, BPPARAM = param)
#saveRDS(dds, "OBJECTS/Myeloid_cells/Monomac1_DESeq.rds")
#res <- results(dds)

### umaps
source("SCRIPTS/gene_umap.R")

markers = c("MAFA")
maf_u = gene_umap(s_obj2, markers, pts=1)

markers = c("SOX4")
sox_u = gene_umap(s_obj2, markers, pts=1)

uSox = arrangeGrob(maf_u+theme_void(), sox_u+theme_void(), ncol=1)



#
#
#

#####
#### mono_mac2

plotter <- s_obj@meta.data
plotter$IL1B <- s_obj@assays$integrated@data["IL1B",]
plotter$LYZ <- s_obj@assays$integrated@data["LYZ",]

plotter <- subset(plotter, MajorPopulations == "mcells_monomac")

plotter$pouch.status2 <- plotter$pouch.status
plotter$pouch.status2[plotter$IL1B < 0.25] <- "other"
plotter$pouch.status2[plotter$LYZ < 0.25] <- "other"

ggIL1B = ggplot(plotter, aes(IL1B, LYZ, color=pouch.status2)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=colorer) +
  theme_bw() + facet_wrap(~pouch.status, nrow=1) +
  geom_vline(xintercept =0.25) + geom_hline(yintercept = 0.25)


plotter_int <- subset(plotter, IL1B > 0.25 & LYZ > 0.25)
samples <- unique(plotter_int$orig.ident)
resser <- data.frame(sampleID = NA, pouchStatus = NA, perc = NA)
for(j in 1:length(samples)){
  curr <- subset(plotter_int, orig.ident == samples[j])
  perc = nrow(curr)/as.numeric(as.character(curr$cell_total[1]))
  adder = data.frame(sampleID = samples[j],
                     pouchStatus = curr$pouch.status[1],
                     perc = perc)
  resser <- rbind(resser, adder)
}
resser=resser[-1,]
foxresser = resser

##stats
ressy=pairwise.wilcox.test(x=foxresser$perc, g=foxresser$pouchStatus, p.adjust.method="fdr")$p.value
pvals=as.numeric(ressy)
comps=paste0(rep(rownames(ressy), ncol(ressy)), ";", rep(colnames(ressy), each=nrow(ressy)))
ressy2=data.frame(Comparison=comps, pvals=pvals)
ressy2$cluster="mono_mac2"

bbIL1B = ggplot(resser, aes(pouchStatus, perc*100, color=pouchStatus)) +
  geom_boxplot(outlier.shape=NA) +
  ggtitle("IL1B+ LYZ+ Cell Percentages") +
  ylab("Cell Percentages") +
  scale_color_manual(values = colorer) +
  geom_jitter(width=0.2) + theme_bw() +
  theme(axis.text.x=element_blank())


###
celler2 = rownames(plotter_int)
s_obj_1 = s_obj2[,celler2]
Idents(s_obj_1) = s_obj_1@meta.data$pouch.status
group_combinations = combn(levels(s_obj_1), m = 2, simplify = TRUE)
s_obj_1@meta.data$Facs_pop = "IL1B"
saveRDS(s_obj_1, "OBJECTS/Myeloid_cells/clusters-MinorPopulations-clust5/diff-expression-orig.ident/IL1B/seurat_obj.rds")

# res = data.frame(p_val=NA, avg_logFC=NA,pct.1=NA,pct.2=NA,p_val_adj=NA, group1=NA, group2=NA, gene=NA)
# for(j in 1:ncol(group_combinations)){
#   g1 = group_combinations[1, j]
#   g2 = group_combinations[2, j]
#   
#   de_genes1 = FindMarkers(s_obj_1, ident.1 = g1, ident.2 = g2, assay = "RNA",
#                           test.use = 'wilcox', logfc.threshold = log(0), 
#                           min.pct = 0.1, only.pos = FALSE,
#                           print.bar = FALSE)
#   de_genes1$group1=g1
#   de_genes1$group2=g2
#   de_genes1$gene=rownames(de_genes1)
#   res = rbind(res, de_genes1)
# }
# res=res[-1,]
# 
# 
# cell_table_de=subset(res, abs(avg_logFC)>0.75 & p_val_adj<0.05)
# good_genes = as.character(unique(cell_table_de$gene))
# cell_table_filt = subset(res, gene %in% good_genes)
# 
# if(nrow(cell_table_filt)>0){
#   de_mat = cast(cell_table_filt, gene ~ group1+group2, value='avg_logFC')
#   rownames(de_mat) = de_mat[,1]
#   de_mat=data.matrix(de_mat[,-1])
#   
#   de_mat[,1]=de_mat[,1]*-1
#   de_mat[,3]=de_mat[,3]*-1
#   colnames(de_mat) = c("Pouchitis vs. Normal Pouch", "UC Inflamed vs. Normal Pouch",
#                        "Pouchitis vs. UC Inflamed")
#   de_mat[is.na(de_mat)]<-0
#   if(nrow(de_mat)>30){row_size=4} else{row_size=8}
#   
#   namer <- paste0("OBJECTS/Myeloid_cells/clusters-MinorPopulations-clust5/diff-expression-orig.ident/IL1B_LYZ_FC_heatmap.pdf")
#   pdf(namer, height = 10, width = 5)
#   print(Heatmap(de_mat, show_row_names = T, show_column_names = T,
#                 #top_annotation = ha,
#                 heatmap_legend_param = list(title = "Fold Change"),
#                 cluster_rows = T, cluster_columns = F, #row_names_side = 'left',
#                 column_names_gp = gpar(fontsize=12),
#                 row_names_gp = gpar(fontsize=row_size),
#                 row_title_gp = gpar(fontsize = 10),
#                 #row_names_max_width = unit(10,'cm'),
#                 use_raster = T,
#                 #cluster_row_slices=F,
#                 #split = cluster2,
#                 #left_annotation = ha,
#                 col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
#   )
#   dev.off()
# }
##
#
###deseq on this pop?
#counts = data.matrix(s_obj@assays$RNA@counts[,rownames(plotter_int)])

#counts = counts[rowSums(counts)>10,]
#counts = counts[rowVars(counts)>0.1,]

#library("BiocParallel")
#param <- BatchtoolsParam(workers=10, cluster="slurm", template="SCRIPTS/slurm.tmpl")

#dds <- DESeqDataSetFromMatrix(counts, data.frame(pouch.status=plotter_int$pouch.status), ~ pouch.status)
#dds <- DESeq(dds, parallel=T, BPPARAM = param)
#saveRDS(dds, "OBJECTS/Myeloid_cells/Monomac2_DESeq.rds")
#res <- results(dds)

##umaps

markers = c("LYZ")
lyz_u = gene_umap(s_obj2, markers, pts=1)

markers = c("IL1B")
il1b_u = gene_umap(s_obj2, markers, pts=1)

uIL1B = arrangeGrob(lyz_u+theme_void(), il1b_u+theme_void(), ncol=1)
#
#

##
#### mono_mac3

plotter <- s_obj@meta.data
plotter$APOE <- s_obj@assays$integrated@data["APOE",]
plotter$C1QC <- s_obj@assays$integrated@data["C1QC",]

plotter <- subset(plotter, MajorPopulations == "mcells_monomac")

plotter$pouch.status2 <- plotter$pouch.status
plotter$pouch.status2[plotter$APOE < 0.25] <- "other"
plotter$pouch.status2[plotter$C1QC < 0.25] <- "other"

ggAPOE = ggplot(plotter, aes(APOE, C1QC, color=pouch.status2)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=colorer) +
  theme_bw() + facet_wrap(~pouch.status, nrow=1) +
  geom_vline(xintercept =0.25) + geom_hline(yintercept = 0.25)


plotter_int <- subset(plotter, APOE > 0.25 & C1QC > 0.25)
samples <- unique(plotter_int$orig.ident)
resser <- data.frame(sampleID = NA, pouchStatus = NA, perc = NA)
for(j in 1:length(samples)){
  curr <- subset(plotter_int, orig.ident == samples[j])
  perc = nrow(curr)/as.numeric(as.character(curr$cell_total[1]))
  adder = data.frame(sampleID = samples[j],
                     pouchStatus = curr$pouch.status[1],
                     perc = perc)
  resser <- rbind(resser, adder)
}
resser=resser[-1,]
foxresser = resser

##stats
ressy=pairwise.wilcox.test(x=foxresser$perc, g=foxresser$pouchStatus, p.adjust.method="fdr")$p.value
pvals=as.numeric(ressy)
comps=paste0(rep(rownames(ressy), ncol(ressy)), ";", rep(colnames(ressy), each=nrow(ressy)))
ressy3=data.frame(Comparison=comps, pvals=pvals)
ressy3$cluster="mono_mac3"

bbAPOE = ggplot(resser, aes(pouchStatus, perc*100, color=pouchStatus)) +
  geom_boxplot(outlier.shape=NA) +
  ggtitle("APOE+ C1QC+ Cell Percentages") +
  ylab("Cell Percentages") +
  scale_color_manual(values = colorer) +
  geom_jitter(width=0.2) + theme_bw() +
  theme(axis.text.x=element_blank())


#
#
#
celler3 = rownames(plotter_int)
s_obj_1 = s_obj2[,celler3]
Idents(s_obj_1) = s_obj_1@meta.data$pouch.status
group_combinations = combn(levels(s_obj_1), m = 2, simplify = TRUE)
s_obj_1@meta.data$Facs_pop = "APOE"
saveRDS(s_obj_1, "OBJECTS/Myeloid_cells/clusters-MinorPopulations-clust5/diff-expression-orig.ident/APOE/seurat_obj.rds")

# res = data.frame(p_val=NA, avg_logFC=NA,pct.1=NA,pct.2=NA,p_val_adj=NA, group1=NA, group2=NA, gene=NA)
# for(j in 1:ncol(group_combinations)){
#   g1 = group_combinations[1, j]
#   g2 = group_combinations[2, j]
#   
#   de_genes1 = FindMarkers(s_obj_1, ident.1 = g1, ident.2 = g2, assay = "RNA",
#                           test.use = 'wilcox', logfc.threshold = log(0), 
#                           min.pct = 0.1, only.pos = FALSE,
#                           print.bar = FALSE)
#   de_genes1$group1=g1
#   de_genes1$group2=g2
#   de_genes1$gene=rownames(de_genes1)
#   res = rbind(res, de_genes1)
# }
# res=res[-1,]
# 
# 
# cell_table_de=subset(res, abs(avg_logFC)>0.75 & p_val_adj<0.05)
# good_genes = as.character(unique(cell_table_de$gene))
# cell_table_filt = subset(res, gene %in% good_genes)
# 
# if(nrow(cell_table_filt)>0){
#   de_mat = cast(cell_table_filt, gene ~ group1+group2, value='avg_logFC')
#   rownames(de_mat) = de_mat[,1]
#   de_mat=data.matrix(de_mat[,-1])
#   
#   de_mat[,1]=de_mat[,1]*-1
#   de_mat[,3]=de_mat[,3]*-1
#   colnames(de_mat) = c("Pouchitis vs. Normal Pouch", "UC Inflamed vs. Normal Pouch",
#                        "Pouchitis vs. UC Inflamed")
#   de_mat[is.na(de_mat)]<-0
#   if(nrow(de_mat)>30){row_size=4} else{row_size=8}
#   
#   namer <- paste0("OBJECTS/Myeloid_cells/clusters-MinorPopulations-clust5/diff-expression-orig.ident/APOE_C1QC_FC_heatmap.pdf")
#   pdf(namer, height = 10, width = 5)
#   print(Heatmap(de_mat, show_row_names = T, show_column_names = T,
#                 #top_annotation = ha,
#                 heatmap_legend_param = list(title = "Fold Change"),
#                 cluster_rows = T, cluster_columns = F, #row_names_side = 'left',
#                 column_names_gp = gpar(fontsize=12),
#                 row_names_gp = gpar(fontsize=row_size),
#                 row_title_gp = gpar(fontsize = 10),
#                 #row_names_max_width = unit(10,'cm'),
#                 use_raster = T,
#                 #cluster_row_slices=F,
#                 #split = cluster2,
#                 #left_annotation = ha,
#                 col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
#   )
#   dev.off()
# }
#
#
##
###deseq on this pop?
#counts = data.matrix(s_obj@assays$RNA@counts[,rownames(plotter_int)])

#counts = counts[rowSums(counts)>10,]
#counts = counts[rowVars(counts)>0.1,]

#library("BiocParallel")
#param <- BatchtoolsParam(workers=10, cluster="slurm", template="SCRIPTS/slurm.tmpl")

#dds <- DESeqDataSetFromMatrix(counts, data.frame(pouch.status=plotter_int$pouch.status), ~ pouch.status)
#dds <- DESeq(dds, parallel=T, BPPARAM = param)
#saveRDS(dds, "OBJECTS/Myeloid_cells/Monomac3_DESeq.rds")
#res <- results(dds)


##umaps

markers = c("C1QC")
c1qc_u = gene_umap(s_obj2, markers, pts=1)

markers = c("APOE")
apoe_u = gene_umap(s_obj2, markers, pts=1)

uAPOE = arrangeGrob(c1qc_u+theme_void(), apoe_u+theme_void(), ncol=1)
#
#
#
#

pdf("FIGURES/FIGURE2/Figure2d_facs.pdf", height = 10, width = 14)
grid.arrange(uSox, ggSox, bbSox,
             uIL1B, ggIL1B, bbIL1B,
             uAPOE, ggAPOE, bbAPOE, nrow=3, widths = c(0.5,1.8,1))
dev.off()

png("FIGURES/FIGURE2/Figure2d_facs.png", height = 10, width = 14, units = "in", res=150)
grid.arrange(uSox, ggSox, bbSox,
             uIL1B, ggIL1B, bbIL1B,
             uAPOE, ggAPOE, bbAPOE, nrow=3, widths = c(0.5,1.8,1))
dev.off()

ressy = rbind(ressy1, ressy2, ressy3)
write.table(ressy, "FIGURES/FIGURE2/Figure2d_facs_stats.txt", sep='\t', quote=F, row.names=F)

#

#
##
#
###### LogFC heatmap between pouch conditions
library(reshape)

de_table = read.table("OBJECTS/Myeloid_cells/clusters-MinorPopulations2-clust6/diff-expression-orig.ident/de.clust6.orig.ident.wilcox.all.csv", T, sep=',')

cells = as.character(unique(de_table$cluster))
for(i in 1:length(cells)){
  cell_table=subset(de_table, cluster==cells[i])
  cell_table_de=subset(cell_table, abs(avg_logFC)>0.75 & p_val_adj<0.05)
  good_genes = as.character(unique(cell_table_de$gene))
  cell_table_filt = subset(cell_table, gene %in% good_genes)
  
  namer1<- paste0("OBJECTS/Myeloid_cells/clusters-MinorPopulations2-clust6/diff-expression-orig.ident/",
                         cells[i], "_FC.txt")
  write.table(cell_table_filt, namer1, row.names=F, sep='\t', quote=F)
  
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
  
  
  namer <- paste0("OBJECTS/Myeloid_cells/clusters-MinorPopulations2-clust6/diff-expression-orig.ident/",
                   cells[i], "_FC_heatmap.pdf")
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




### SPP1
genes = c("SPP1")
vln_plot = 
  VlnPlot(
    s_obj2, features = genes, 
    combine = TRUE, cols = colorer, ncol = 1
  ) 

gradient_cols = c("grey85", "purple4")

feat_plot = 
  FeaturePlot(
    s_obj2, features = genes, 
    combine = TRUE, cols = gradient_cols, ncol = 1
  ) 

plt = grid.arrange(vln_plot, feat_plot, nrow=1)

ggsave("OBJECTS/Myeloid_cells/SPP1_vln_plot.pdf", plot=plt, height = 4, width = 10, units='in')

