library(DESeq2)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)

colors_clusters = c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))

big_colorer <- c(act_CD8_NR4A2_hi=colors_clusters[1], 
                 act_CD8_NR4A2_lo=colors_clusters[2], 
                 activated_CD4=colors_clusters[3], 
                 CD4_fos_hi=colors_clusters[4], 
                 CD4_fos_lo=colors_clusters[5],
                 Foxp3_cells="darkseagreen1", 
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
                 mono_mac2="coral2", 
                 mono_mac1="plum2", 
                 mono_mac3="forestgreen")

### etrolizumab

counts = read.table("additional_data/Tew_et_al/etrolizumab_counts.txt", T, sep='\t')
rownames(counts) = counts[,1]
counts = counts[,-1]

meta = read.table("additional_data/Tew_et_al/etrolizumab_meta.txt", T, '\t')
rownames(meta) <- as.character(meta$Sample_geo_accession)
meta = meta[colnames(counts),]

meta = subset(meta, Remission != "N/A")
counts = counts[,rownames(meta)]

cds = DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~Remission)
cds
dds = DESeq(cds)

res = data.frame(results(dds, contrast = c("Remission", "Remitter", "Non-remitter")))
res$Gene = rownames(res)

ct_ref = read.table("OBJECTS/Signature/signatureRefCt.txt", T, "\t")
rownames(ct_ref) <- ct_ref$gene
ct_ref = ct_ref[rownames(res),]

res$cellState = ct_ref$cluster
res=subset(res, !is.na(cellState))

res = res[order(abs(res$log2FoldChange), decreasing=T),]

write.table(res, "FIGURES/SF8/Etro_DE_results.txt", row.names=F, sep='\t', quote=F)


res_plot = res[1:40,]
res_plot = res_plot[order(res_plot$log2FoldChange, decreasing=T),]

res_plot$Gene = factor(res_plot$Gene, levels = as.character(res_plot$Gene))
etro_p=ggplot(res_plot, aes(x=Gene, y = log2FoldChange, fill=cellState)) + 
  geom_col() + coord_flip() + ggtitle("Etrolizumab Responders (n=12) vs. Non-responders (n=58)") +
  scale_fill_manual(values = big_colorer) + theme_bw()

counts = counts(dds, normalized = TRUE)
mono_mac2_g = as.character(subset(ct_ref, cluster=="mono_mac2")$gene)
counts = counts[mono_mac2_g,]

ha = columnAnnotation(df=data.frame(Remission=meta$Remission),
                      col=list(Remission=c("Non-remitter"="red3", "Remitter"="steelblue2")))

heatty_e = grid.grabExpr(draw(
  Heatmap(t(scale(t(counts))), show_row_names = T, show_column_names = F,
          heatmap_legend_param = list(title = "Scaled\nValue"),
          cluster_rows = F, cluster_columns = T, row_names_side = 'left',
          column_names_gp = gpar(fontsize=10),
          column_title = "Etrolizumab Responders (n=12) vs. Non-responders (n=58)",
          top_annotation = ha,
          #row_names_gp = gpar(fontsize=12),
          column_names_rot = 60,
          #row_title_gp = gpar(fontsize = 10),
          #row_names_max_width = unit(10,'cm'),
          column_split = meta$Remission,
          use_raster = F,
          col = colorRamp2(c(-2,0,2), c("blue","white" ,"red"))),
  padding = unit(c(1, 1, 1, 1), "cm")
))







#
#

#
#
#
#
##
#microarray
# vedo
library(limma)
table = read.table("additional_data/Arijs_et_al/vedolizumab_counts.txt", T, '\t')
rownames(table) <- table[,1]
table <- table[,-1]

meta = read.table("additional_data/Arijs_et_al/old_vedolizumab_meta.txt", T, '\t')
rownames(meta) = meta$Sample_geo_accession
table = table[,rownames(meta)]

group = gsub("-", "_", meta$Remission)

design <- model.matrix(~0 +group)

fit <- lmFit(table, design)
contrast.matrix <- makeContrasts(groupResponder-groupNon_responder, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit2)
fit_filt <- subset(topTable(fit, n=Inf), adj.P.Val < 0.1)
#fit_filt <- fit_filt[order(fit_filt$logFC),]

fit_filt$Gene = rownames(fit_filt)

ct_ref = read.table("OBJECTS/Signature/signatureRefCt.txt", T, "\t")
rownames(ct_ref) <- ct_ref$gene
ct_ref = ct_ref[rownames(fit_filt),]

fit_filt$cellState = ct_ref$cluster
fit_filt=subset(fit_filt, !is.na(cellState))

fit_filt = fit_filt[order(abs(fit_filt$logFC), decreasing=T),]

write.table(fit_filt, "FIGURES/FIGURE6/Vedo_DE_results.txt", row.names=F, sep='\t', quote=F)

res_plot = fit_filt[1:40,]
res_plot = res_plot[order(res_plot$logFC, decreasing=T),]

res_plot$Gene = factor(res_plot$Gene, levels = as.character(res_plot$Gene))
vedo_p=ggplot(res_plot, aes(x=Gene, y = logFC, fill=cellState)) + 
  geom_col() + coord_flip() + ggtitle("Vedolizumab Responders (n=17) vs. Non-responders (n=47)") +
  scale_fill_manual(values = big_colorer) + theme_bw()

#subset the counts
meta2 = subset(meta, grepl(c("Responder|Non-responder"), meta$Remission))

mono_mac2_g = as.character(subset(ct_ref, cluster=="mono_mac2")$gene)
counts = table[mono_mac2_g,rownames(meta2)]

ha = columnAnnotation(df=data.frame(Remission=meta2$Remission),
                      col=list(Remission=c("Non-responder"="red3", "Responder"="steelblue2")))

heatty_v = grid.grabExpr(draw(
  Heatmap(t(scale(t(counts))), show_row_names = T, show_column_names = F,
          heatmap_legend_param = list(title = "Scaled\nValue"),
          cluster_rows = F, cluster_columns = T, row_names_side = 'left',
          column_names_gp = gpar(fontsize=10),
          top_annotation = ha,
          #row_names_gp = gpar(fontsize=12),
          column_title = "Vedolizumab Responders (n=17) vs. Non-responders (n=47)",
          column_names_rot = 60,
          #row_title_gp = gpar(fontsize = 10),
          #row_names_max_width = unit(10,'cm'),
          column_split = meta2$Remission,
          use_raster = F,
          col = colorRamp2(c(-2,0,2), c("blue","white" ,"red"))),
  padding = unit(c(1, 1, 1, 1), "cm")
))


##quick IL1B stats

il1b=as.numeric(counts["IL1B",])
il1b=data.frame(Remission=meta2$Remission, IL1B=il1b)

ggplot(il1b, aes(Remission, IL1B)) + geom_boxplot() + geom_jitter(width = 0.2)

#
#
#
#
# gol

table = read.table("additional_data/Li_et_al/golimumab_counts.txt", T, '\t')
rownames(table) <- table[,1]
table <- table[,-1]

meta = read.table("additional_data/Li_et_al/golimumab_meta.txt", T, '\t')
rownames(meta) = meta$Sample_geo_accession
table = table[,rownames(meta)]

group = meta$Remission

design <- model.matrix(~0 +group)

fit <- lmFit(table, design)
contrast.matrix <- makeContrasts(groupYes-groupNo, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit2)
fit_filt <- subset(topTable(fit, n=Inf))
#fit_filt <- fit_filt[order(fit_filt$logFC),]

fit_filt$Gene = rownames(fit_filt)

ct_ref = read.table("OBJECTS/Signature/signatureRefCt.txt", T, "\t")
rownames(ct_ref) <- ct_ref$gene
ct_ref = ct_ref[rownames(fit_filt),]

fit_filt$cellState = ct_ref$cluster
fit_filt=subset(fit_filt, !is.na(cellState))

fit_filt = fit_filt[order(abs(fit_filt$logFC), decreasing=T),]

write.table(fit_filt, "FIGURES/SF8/Gol_DE_results.txt", row.names=F, sep='\t', quote=F)

res_plot = fit_filt[1:40,]
res_plot = res_plot[order(res_plot$logFC, decreasing=T),]

res_plot$Gene = factor(res_plot$Gene, levels = as.character(res_plot$Gene))
gol_p=ggplot(res_plot, aes(x=Gene, y = logFC, fill=cellState)) + 
  geom_col() + coord_flip() + ggtitle("Golimumab Responders (n=29) vs. Non-responders (n=27)") +
  scale_fill_manual(values = big_colorer) + theme_bw()


#subset the counts
meta2 = subset(meta, grepl(c("Yes|No"), meta$Remission))

mono_mac2_g = as.character(subset(ct_ref, cluster=="mono_mac2")$gene)
counts = table[mono_mac2_g,rownames(meta2)]

ha = columnAnnotation(df=data.frame(Remission=meta2$Remission),
                      col=list(Remission=c("No"="red3", "Yes"="steelblue2")))


heatty_g = grid.grabExpr(draw(
  Heatmap(t(scale(t(counts))), show_row_names = T, show_column_names = F,
          heatmap_legend_param = list(title = "Scaled\nValue"),
          cluster_rows = F, cluster_columns = T, row_names_side = 'left',
          column_names_gp = gpar(fontsize=10),
          top_annotation = ha,
          column_title = c("Golimumab Responders (n=29) vs. Non-responders (n=27)"),
          #row_names_gp = gpar(fontsize=12),
          column_names_rot = 60,
          #row_title_gp = gpar(fontsize = 10),
          #row_names_max_width = unit(10,'cm'),
          column_split = meta2$Remission,
          use_raster = F,
          col = colorRamp2(c(-2,0,2), c("blue","white" ,"red"))),
  padding = unit(c(1, 1, 1, 1), "cm")
))



pdf("FIGURES/SF8/FigureSF8_clinicalTrials.pdf", height = 10, width = 20)
grid.arrange(etro_p, vedo_p, gol_p, nrow=1)
dev.off()

pdf("FIGURES/FIGURE6/Figure6C_clinicalTrials.pdf", height = 8, width = 20)
grid.arrange(heatty_e, heatty_v, heatty_g, nrow=1)
dev.off()


