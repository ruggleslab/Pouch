library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(gridExtra)

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
                 DC='darkseagreen',
                 mono_mac1="purple", 
                 mono_mac2="coral", 
                 mono_mac3="forestgreen")

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
              "Th17_cells",
              "cycling_b_cells", 
              "follicular_cells", 
              "gc_cells", 
              "plasma_NFKBIA_hi", 
              "plasma_NFKBIA_lo",
              "mast_cells", 
              "pdcs", 
              "DC",
              "mono_mac1", 
              "mono_mac2",
              "mono_mac3")

## read in a Marker DE table for each of 22 cell state (T,B and Myeloid)

T_ct <- read.table("OBJECTS/T_cells/clusters-MinorPopulations-clust12/markers-global/markers.clust12.wilcox.all.csv",
                   T, ",")

B_ct <- read.table("OBJECTS/B_cells/clusters-MinorPopulations-clust5/markers-global/markers.clust5.wilcox.all.csv",
                   T, ",")

M_ct <- read.table("OBJECTS/Myeloid_cells/clusters-MinorPopulations2-clust6/markers-global/markers.clust6.wilcox.all.csv",
                   T, ",")

#
#
##
#top 10%

ct_de <- rbind(T_ct, B_ct, M_ct)
ct_de$cluster = factor(ct_de$cluster, levels = level_fix)
ct_de=subset(ct_de, avg_logFC > 0)

lower=c("cycling_b_cells", "follicular_cells", 
        "plasma_NFKBIA_lo", "gc_cells", "plasma_NFKBIA_hi", "pdcs")
higher=c("act_CD8_NR4A2_hi", "activated_CD4", 'CD4_fos_lo', "memory_CD4", "memory_CD8", "naive_CD8")
higherh=c("CD4_fos_hi")

ct_final=ct_de[1,]
for(j in 1:length(level_fix)){
  curr=subset(ct_de, cluster==level_fix[j])
  rounder=as.numeric(round(quantile(1:nrow(curr),0.1)))
  if(level_fix[j] %in% lower){
    rounder=as.numeric(round(quantile(1:nrow(curr),0.04)))
  }
  if(level_fix[j] %in% higher){
    rounder=as.numeric(round(quantile(1:nrow(curr),0.25)))
  }
  if(level_fix[j] %in% higherh){
    rounder=as.numeric(round(quantile(1:nrow(curr),0.7)))
  }
  ct_final=rbind(ct_final,curr[1:rounder,])
}
ct_final=ct_final[-1,]


###
#
#


# ct_de <- rbind(T_ct, B_ct, M_ct)
# ct_de$cluster = factor(ct_de$cluster, levels = level_fix)
# ct_de = subset(ct_de, avg_logFC > 0.4 & p_val_adj < 0.01)
# 
# red = names(table(ct_de$cluster))[table(ct_de$cluster)>52]
# red_term = paste(red, collapse = "|")
# ct_add = subset(ct_de, grepl(red_term, ct_de$cluster))
# ct_add1 = subset(ct_add, avg_logFC > 1.2 & p_val_adj < 0.01)
# 
# suppl = names(table(ct_add1$cluster))[table(ct_add1$cluster)<15]
# supp_term = paste(suppl, collapse = "|")
# ct_add = subset(ct_de, grepl(supp_term, ct_de$cluster))
# ct_add2 = subset(ct_add, avg_logFC > 0.9 & p_val_adj < 0.01)
# 
# suppl = names(table(ct_add1$cluster))[table(ct_add2$cluster)<15 & table(ct_add1$cluster)<15]
# supp_term = paste(suppl, collapse = "|")
# ct_add = subset(ct_de, grepl(supp_term, ct_de$cluster))
# ct_add3 = subset(ct_add, avg_logFC > 0.7 & p_val_adj < 0.01)
# 
# suppl = names(table(ct_de$cluster))[table(ct_de$cluster)<52]
# supp_term = paste(suppl, collapse = "|")
# ct_add = subset(ct_de, grepl(supp_term, ct_de$cluster))
# #ct_add4 = subset(ct_add, avg_logFC > 0.7 & p_val_adj < 0.01)
# 
# ct_final = rbind(ct_add, ct_add1, ct_add2, ct_add3)
# #ct_final = ct_final[duplicated(ct_final)==FALSE,]
# #table(ct_final$cluster)


dups = as.character(ct_final$gene[duplicated(ct_final$gene)])
remover=c()
for(i in 1:length(dups)){ 
  checker = subset(ct_final, gene == dups[i])
  checker = checker[order(checker$avg_logFC, decreasing=T),]
  
  remover = c(remover, rownames(checker)[-1])
}
remover = unique(remover)
ct_final2 = ct_final[setdiff(rownames(ct_final), remover),]

ct_final2 = ct_final2[order(ct_final2$cluster),]
write.table(ct_final2, "OBJECTS/Signature/signatureRefCt.txt", row.names=F, sep='\t', quote=F)

####
s_obj <- readRDS("OBJECTS/seurat_obj.rds")
m_obj <- readRDS("OBJECTS/Myeloid_cells/seurat_obj.rds")
t_obj <- readRDS("OBJECTS/T_cells/seurat_obj.rds")
b_obj <- readRDS("OBJECTS/B_cells/seurat_obj.rds")

MinorPop=data.frame(
  cells=c(colnames(m_obj), colnames(b_obj), colnames(t_obj)),
  pops=c(m_obj@meta.data$MinorPopulations2, b_obj@meta.data$MinorPopulations, t_obj@meta.data$MinorPopulations)
)
rownames(MinorPop)=as.character(MinorPop$cells)
MinorPop=MinorPop[colnames(s_obj),]
s_obj@meta.data$MinorPopulations=MinorPop$pops

counts <- s_obj@assays$RNA@data[as.character(ct_final2$gene),]

count_agg <- aggregate(t(data.matrix(counts)), list(s_obj@meta.data$MinorPopulations), mean)
rownames(count_agg) <- count_agg[,1]
sigmat = t(count_agg[,-1])

sigmat = sigmat[,level_fix]



#
#
##### are any of these signature genes also GWAS associated?
gwas_genes <- as.character(read.table("additional_data/IBD_GWAS_genes.txt", F, '\t')[,1])
good_gwas = intersect(rownames(sigmat), gwas_genes)

gwasser = c()
new_namer = c()
for(i in 1:nrow(sigmat)){
  curr = rownames(sigmat)[i]
  adder = "No"
  adder2 = ""
  if(curr %in% good_gwas){
    adder = "Yes"
    adder2=curr
  }
  gwasser = c(gwasser, adder)
  new_namer = c(new_namer, adder2)
}

#####

ha = rowAnnotation(df=data.frame(cellState=ct_final2$cluster, GWAS=gwasser),
                   col = list(cellState=big_colorer,
                              GWAS=c("Yes"="black", "No"="white")))

heatty_small = grid.grabExpr(draw(
  Heatmap(t(scale(t(sigmat))), show_row_names = F, show_column_names = T,
          #top_annotation = ha2,
          heatmap_legend_param = list(title = "Scaled\nValue"),
          cluster_rows = F, cluster_columns = F, row_names_side = 'left',
          column_names_gp = gpar(fontsize=12),
          row_names_gp = gpar(fontsize=5),
          #column_names_rot = 60,
          row_title_gp = gpar(fontsize = 10),
          row_names_max_width = unit(10,'cm'),
          use_raster = T,
          #cluster_row_slices=F,
          #row_split = ct_de_sig$CellState,
          col = colorRamp2(c(-2,0,2), c("blue","white" ,"red")))+ha,
  padding = unit(c(2, 2, 3, 2), "cm")
))


namer <- "FIGURES/FIGURE5/Figure5a_sigMatrix_small.pdf"
pdf(namer, height = 10, width = 12)
grid.arrange(heatty_small)
dev.off()

sig_write = data.frame(GeneSymbol = rownames(sigmat), sigmat)
write.table(sig_write, "OBJECTS/Signature/PouchUCI_22_signatures.txt", sep='\t', row.names=F, quote=F)

#rownames(sigmat) = new_namer
heatty_large = grid.grabExpr(draw(
  Heatmap(t(scale(t(sigmat))), show_row_names = T, show_column_names = T,
          #top_annotation = ha2,
          heatmap_legend_param = list(title = "Scaled\nValue"),
          cluster_rows = F, cluster_columns = F, row_names_side = 'left',
          column_names_gp = gpar(fontsize=12),
          row_names_gp = gpar(fontsize=5),
          #column_names_rot = 60,
          row_title_gp = gpar(fontsize = 10),
          row_names_max_width = unit(10,'cm'),
          use_raster = T,
          #cluster_row_slices=F,
          #row_split = ct_de_sig$CellState,
          col = colorRamp2(c(-2,0,2), c("blue","white" ,"red")))+ha,
  padding = unit(c(2, 2, 3, 2), "cm")
))


namer <- "FIGURES/FIGURE5/Figure5a_sigMatrix_large.pdf"
pdf(namer, height = 45, width = 10)
grid.arrange(heatty_large)
dev.off()

