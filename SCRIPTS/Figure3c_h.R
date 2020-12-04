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
                 Th17_cells=colors_clusters[12],
                 cycling_b_cells=colors_clusters[13], 
                 follicular_cells=colors_clusters[14], 
                 gc_cells=colors_clusters[15], 
                 plasma_NFKBIA_hi=colors_clusters[16], 
                 plasma_NFKBIA_lo=colors_clusters[17],
                 mast_cells="navy", 
                 pdcs="pink", 
                 m1_HLA_DR_macs="yellow", 
                 m1_macs="purple", 
                 m2_macs="forestgreen")

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



####### population level boxplots first

s_obj <- readRDS("OBJECTS/T_cells/seurat_obj.rds")
s_obj@meta.data$MinorPopulations <- factor(s_obj@meta.data$MinorPopulations, 
                                           levels = level_fix)

cluster = s_obj@meta.data[,"MinorPopulations"]

Idents(s_obj) <- s_obj@meta.data$pouch.status

barplot = data.frame(PouchStatus=Idents(s_obj), Cluster = cluster)
barplot$Sample <- s_obj@meta.data$orig.ident
bb=ggplot(barplot, aes(PouchStatus, fill=cluster)) + geom_bar() + scale_fill_manual(values=colors_clusters) + theme_bw() + ylab("Cell count")
bb2=ggplot(barplot, aes(PouchStatus, fill=cluster)) + geom_bar(position="fill") + scale_fill_manual(values=colors_clusters) + theme_bw() + ylab("Proportion")

df = table(barplot$Cluster, barplot$Sample, barplot$PouchStatus)
dfm = melt(df)
dfm = subset(dfm, value > 0)
colnames(dfm) <- c("Cluster", "Sample", "PouchStatus", "value")
dfm$Cluster <- factor(dfm$Cluster, levels = level_fix)
dfm$PouchStatus <- factor(dfm$PouchStatus, levels = c("Normal_pouch", "Pouchitis", "UC_inflamed"))
melt_plot1 <- ggplot(dfm, aes(PouchStatus, value, color=PouchStatus)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(width = 0.2) + xlab("PouchStatus") +
  theme_bw() + scale_color_manual(values = colorer) + ylab("Cell Count") + 
  theme(axis.text.x=element_blank()) +
  facet_wrap(~Cluster, nrow=2, scales="free_y")

dfm$totals <- 0
samps <- as.character(unique(dfm$Sample))
for(i in 1:length(samps)){
  curr = subset(s_obj@meta.data, orig.ident == samps[i])
  totaler = as.numeric(curr$cell_total[1])
  dfm$totals[dfm$Sample == samps[i]] <- totaler
}

dfm$totals <- as.numeric(dfm$totals)
dfm$perc <- (dfm$value/dfm$totals)*100
write_to <- paste0("OBJECTS/T_cells/clusters-resolutions/clust_box_MinorPopulations_values.txt")
write.table(dfm, write_to, sep='\t', row.names=F, quote=F)

dfm$Cluster <- factor(dfm$Cluster, levels = level_fix)

melt_plot2 <- ggplot(dfm, aes(PouchStatus, perc, color=PouchStatus)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(width = 0.2) + xlab("PouchStatus") +
  theme_bw() + scale_color_manual(values = colorer) + ylab("Cell Percentage") + 
  theme(axis.text.x=element_blank()) +
  facet_wrap(~Cluster, nrow=2, scales="free_y")

bb3 = ggplot(barplot, aes(cluster, fill=PouchStatus)) + 
  geom_bar() + scale_fill_manual(values=colorer) + 
  theme_bw() + ylab("Cell Count") + theme(axis.text.x = element_text(angle = 90))
bb4 = ggplot(barplot, aes(cluster, fill=PouchStatus)) + 
  geom_bar(position="fill") + scale_fill_manual(values=colorer) + 
  theme_bw() + ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 90))
bb = arrangeGrob(bb,bb2,bb3,bb4, ncol=2)

#write_to <- list.files(args[1])[grep(gsub("\\.","", args[2]), list.files(args[1]))]
write_to <- paste0("OBJECTS/T_cells/clusters-resolutions/clust_bar_MinorPopulations_barplot.pdf")

ggsave(write_to, plot = bb, width = 10, height = 10, units = "in")

mm = arrangeGrob(melt_plot2, ncol=1)
write_to <- paste0("FIGURES/FIGURE3/Figure3b_MinorPopulations_boxplot.pdf")
ggsave(write_to, plot = mm, width = 14, height = 6, units = "in")

write_to <- paste0("FIGURES/SF2/FigureSF2B_MinorPopulations_boxplot.pdf")
ggsave(write_to, plot = mm, width = 14, height = 6, units = "in")

### grab foxp3, memory_cd8 and naive cd8 boxes for later

# dfm$PouchStatus = factor(dfm$PouchStatus, levels = c("Normal_pouch", "Pouchitis", "UC_inflamed"))
# foxp3 = ggplot(subset(dfm, Cluster == "Foxp3_cells"), aes(PouchStatus, perc, color=PouchStatus)) + 
#   geom_boxplot(outlier.shape=NA) + geom_jitter(width = 0.2) + xlab("PouchStatus") +
#   theme_bw() + scale_color_manual(values = colorer) + ylab("Cell Percentage") + 
#   ggtitle("Foxp3 Cells") +
#   theme(axis.text.x=element_blank())
# 
# cd8s = ggplot(subset(dfm, Cluster == "memory_CD8" | Cluster == "naive_CD8"), aes(PouchStatus, perc, color=PouchStatus)) + 
#   geom_boxplot(outlier.shape=NA) + geom_jitter(width = 0.2) + xlab("PouchStatus") +
#   theme_bw() + scale_color_manual(values = colorer) + ylab("Cell Percentage") + 
#   ggtitle("CD8+ Cells") +
#   theme(axis.text.x=element_blank())


### UMAP markers
source("SCRIPTS/gene_umap.R")

markers1 = c("CCL5", "IFNG", "TNF", "CD4", "IL7R", "FOXP3") 
markers2 = c("GIMAP7", "SPRY1", "CD8B", "GZMB", "TOX2", "CCL20")
uplots1 = gene_umap(s_obj, markers1, pts=1)
uplots2 = gene_umap(s_obj, markers2, pts=1)

uplots = arrangeGrob(uplots1+theme_void(), uplots2+theme_void())

ggsave("FIGURES/SF2/FigureSF2B.umap.png", plot = uplots, width = 22, height = 7, units = "in")
ggsave("FIGURES/FIGURE3/Figure3B.umap.png", plot = uplots, width = 22, height = 7, units = "in")


#
#
######





s_obj2 = s_obj
s_obj = readRDS("OBJECTS/seurat_obj.rds")
meta <- s_obj@meta.data

plotter <- s_obj@meta.data
plotter$FOXP3 <- s_obj@assays$integrated@data["FOXP3",]
plotter$BATF <- s_obj@assays$integrated@data["BATF",]

plotter <- subset(plotter, MajorPopulations == "tcells")

plotter$pouch.status2 <- plotter$pouch.status
plotter$pouch.status2[plotter$FOXP3 < 0.25] <- "other"
plotter$pouch.status2[plotter$BATF < 0.25] <- "other"

gg = ggplot(plotter, aes(FOXP3, BATF, color=pouch.status2)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=colorer) +
  theme_bw() + facet_wrap(~pouch.status, nrow=1) +
  geom_vline(xintercept =0.25) + geom_hline(yintercept = 0.25)


plotter_int <- subset(plotter, FOXP3 > 0.25 & BATF > 0.25)
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
ressy1$cluster="Foxp3"

bb = ggplot(resser, aes(pouchStatus, perc*100, color=pouchStatus)) +
  geom_boxplot(outlier.shape=NA) +
  ggtitle("FOXP3+ BATF+ Cell Percentages") +
  ylab("Cell Percentages") +
  scale_color_manual(values = colorer) +
  geom_jitter(width=0.2) + theme_bw() +
  theme(axis.text.x=element_blank())


#
#
#
celler = rownames(plotter_int)
s_obj_1 = s_obj2[,celler]
Idents(s_obj_1) = s_obj_1@meta.data$pouch.status
group_combinations = combn(levels(s_obj_1), m = 2, simplify = TRUE)
s_obj_1@meta.data$Facs_pop = "FOXP3"
saveRDS(s_obj_1, "OBJECTS/T_cells/clusters-MinorPopulations-clust12/diff-expression-orig.ident/FOXP3/seurat_obj.rds")

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
#   namer <- paste0("OBJECTS/T_cells/clusters-MinorPopulations-clust12/diff-expression-orig.ident/FOXP3_BATF_FC_heatmap.pdf")
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

###deseq on this pop?
#counts = data.matrix(s_obj@assays$RNA@counts[,rownames(plotter_int)])

#counts = counts[rowSums(counts)>10,]
#counts = counts[rowVars(counts)>0.1,]

#library("BiocParallel")
#param <- BatchtoolsParam(workers=10, cluster="slurm", template="SCRIPTS/slurm.tmpl")

#dds <- DESeqDataSetFromMatrix(counts, data.frame(pouch.status=plotter_int$pouch.status), ~ pouch.status)
#dds <- DESeq(dds, parallel=T, BPPARAM = param)
#saveRDS(dds, "OBJECTS/T_cells/Foxp3_DESeq.rds")
#res <- results(dds)

### umaps
source("SCRIPTS/gene_umap.R")

markers = c("BATF")
batf_u = gene_umap(s_obj2, markers, pts=0.3)

markers = c("FOXP3")
fox_u = gene_umap(s_obj2, markers, pts=0.3)

uFox = arrangeGrob(batf_u+theme_void(), fox_u+theme_void(), ncol=1)

##

pdf("FIGURES/FIGURE3/Figure3cde_facs.pdf", height = 4, width = 14)
grid.arrange(uFox, gg,bb, nrow=1, widths = c(0.5,1.8,1))
dev.off()

png("FIGURES/FIGURE3/Figure3cde_facs.png", height = 4, width = 14, units = "in", res=150)
grid.arrange(uFox, gg,bb, nrow=1, widths = c(0.5,1.8,1))
dev.off()





#
#
#
#
#####
#### CD8 t cells atf3, il2 and cd69


plotter <- s_obj@meta.data
plotter$CD8A <- s_obj@assays$integrated@data["CD8A",]
plotter$ATF3 <- s_obj@assays$integrated@data["ATF3",]
plotter$IL2 <- s_obj@assays$integrated@data["IL2",]
plotter$CD69 <- s_obj@assays$integrated@data["CD69",]
plotter$IFNG <- s_obj@assays$integrated@data["IFNG",]

plotter <- subset(plotter, MajorPopulations == "tcells")

plotter$pouch.status2 <- plotter$pouch.status
#plotter$pouch.status2[plotter$IFNG < 1.25] <- "other"
plotter$pouch.status2[plotter$IL2 < 0.25] <- "other"
plotter$pouch.status2[plotter$CD8A < 0.25] <- "other"

ggCD8 = ggplot(plotter, aes(IL2, CD8A, color=pouch.status2)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=colorer) +
  theme_bw() + facet_wrap(~pouch.status, nrow=1) +
  geom_vline(xintercept = 0.25) + geom_hline(yintercept = 0.25)


#plotter_int <- subset(plotter, IFNG > 1.25 & CD8A > 0.5)
plotter_int <- subset(plotter, IL2 > 0.25 & CD8A > 0.25)
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
ressy2=data.frame(Comparison=comps, pvals=pvals)
ressy2$cluster="CD8/IL2"

bb1 = ggplot(resser, aes(pouchStatus, perc*100, color=pouchStatus)) +
  geom_boxplot(outlier.shape=NA) +
  ggtitle("CD8a+ IL2+ Cell Percentages") +
  ylab("Cell Percentages") +
  scale_color_manual(values = colorer) +
  geom_jitter(width=0.2) + theme_bw() +
  theme(axis.text.x=element_blank())

#
#
#
celler2 = rownames(plotter_int)
s_obj_1 = s_obj2[,celler2]
Idents(s_obj_1) = s_obj_1@meta.data$pouch.status
group_combinations = combn(levels(s_obj_1), m = 2, simplify = TRUE)
s_obj_1@meta.data$Facs_pop = "CD8A"
saveRDS(s_obj_1, "OBJECTS/T_cells/clusters-MinorPopulations-clust12/diff-expression-orig.ident/CD8A/seurat_obj.rds")

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

res = res = read.table("OBJECTS/T_cells/clusters-MinorPopulations-clust12/diff-expression-orig.ident/CD8A/clusters-Facs_pop-clust1/diff-expression-orig.ident/de.clust1.orig.ident.wilcox.all.csv", sep=',', header=T)
cell_table_de=subset(res, abs(avg_logFC)>0.75 & p_val_adj<0.05)
good_genes = as.character(unique(cell_table_de$gene))
cell_table_filt = subset(res, gene %in% good_genes)

write.table(cell_table_filt,
            "OBJECTS/T_cells/clusters-MinorPopulations-clust12/diff-expression-orig.ident/CD8A/clusters-Facs_pop-clust1/diff-expression-orig.ident/CD8A_FC.txt", sep='\t', row.names=F, quote=F)

if(nrow(cell_table_filt)>0){
  de_mat = cast(cell_table_filt, gene ~ group1+group2, value='avg_logFC')
  rownames(de_mat) = de_mat[,1]
  de_mat=data.matrix(de_mat[,-1])

  de_mat[,1]=de_mat[,1]*-1
  de_mat[,3]=de_mat[,3]*-1
  colnames(de_mat) = c("Pouchitis vs. Normal Pouch", "UC Inflamed vs. Normal Pouch",
                       "Pouchitis vs. UC Inflamed")
  de_mat[is.na(de_mat)]<-0
  if(nrow(de_mat)>30){row_size=4} else{row_size=8}

  namer <- paste0("OBJECTS/T_cells/clusters-MinorPopulations-clust12/diff-expression-orig.ident/CD8A/clusters-Facs_pop-clust1/diff-expression-orig.ident/CD8A_FC_heatmap.pdf")
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


### umaps
source("SCRIPTS/gene_umap.R")

markers = c("CD8A")
cd8_u = gene_umap(s_obj2, markers, pts=0.3)

#markers = c("IFNG")
#ifng_u = gene_umap(s_obj2, markers, pts=0.5)

markers = c("IL2")
il2_u = gene_umap(s_obj2, markers, pts=0.3)

uCD8 = arrangeGrob(cd8_u+theme_void(), il2_u+theme_void(), ncol=1)


pdf("FIGURES/FIGURE3/Figure3fgh_facs.pdf", height = 4, width = 14)
grid.arrange(uCD8, ggCD8, bb1, nrow=1, widths = c(0.5,1.8,1))
dev.off()

png("FIGURES/FIGURE3/Figure3fgh_facs.png", height = 4, width = 14, units = "in", res=150)
grid.arrange(uCD8, ggCD8, bb1, nrow=1, widths = c(0.5,1.8,1))
dev.off()


ressy = rbind(ressy1, ressy2)
write.table(ressy, "FIGURES/FIGURE3/Figure3_facs_stats.txt", sep='\t', quote=F)


##
#
#
#

#
#
##
###### LogFC heatmap between pouch conditions
library(reshape)

de_table = read.table("OBJECTS/T_cells/clusters-MinorPopulations-clust12/diff-expression-orig.ident/de.clust12.orig.ident.wilcox.all.csv", T, sep=',')

cells = as.character(unique(de_table$cluster))
for(i in 1:length(cells)){
  cell_table=subset(de_table, cluster==cells[i])
  cell_table_de=subset(cell_table, abs(avg_logFC)>0.6 & p_val_adj<0.05)
  good_genes = as.character(unique(cell_table_de$gene))
  cell_table_filt = subset(cell_table, gene %in% good_genes)
  
  namer1<- paste0("OBJECTS/T_cells/clusters-MinorPopulations-clust12/diff-expression-orig.ident/",
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
    if(nrow(de_mat)>50){row_size=4} else{row_size=8}
    
    namer <- paste0("OBJECTS/T_cells/clusters-MinorPopulations-clust12/diff-expression-orig.ident/",
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



#
#
#
#
#
