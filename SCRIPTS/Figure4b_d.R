library(ggsci)
library(ggplot2)
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
                 mono_mac1="purple", 
                 mono_mac2="coral", 
                 mono_mac3="forestgreen")


library(Matrix)
mat = readMM("additional_data/Smillie_et_al/gene_sorted-Imm.matrix.mtx")

genes = as.character(read.table("additional_data/Smillie_et_al/Imm.genes.tsv", F, "\t")[,1])
barcodes = as.character(read.table("additional_data/Smillie_et_al/Imm.barcodes2.tsv", F, '\t')[,1])

colnames(mat) = barcodes
rownames(mat) = genes

sig = read.table("OBJECTS/Signature/signatureRefCt.txt", T, '\t')
overlapper = intersect(as.character(sig$gene), rownames(mat))
counts2 = as.matrix(mat[overlapper,])

meta = read.table("additional_data/Smillie_et_al/ramnik_meta.txt", T,'\t')
meta = meta[-1,]
rownames(meta) <- as.character(meta$NAME)
meta = meta[colnames(counts2),]


#
#
#
###
counts_agg <- aggregate(t(counts2), list(meta$Sample, meta$Cluster, meta$Health), mean)
count_summ = t(counts_agg[,4:ncol(counts_agg)])
colnames(count_summ) = paste0(counts_agg$Group.2, "_", counts_agg$Group.3, "_", counts_agg$Group.1)
meta_sum = counts_agg[,1:3]
colnames(meta_sum) <- c("SampleID", "cellState", "Inflammation")

library(ComplexHeatmap)
library(gridExtra)
library(circlize)


count_simp = aggregate(t(count_summ), by=list(meta_sum$cellState), mean)
rownames(count_simp) = count_simp[,1]
count_simp = t(count_simp[,-1])

counts_ours = read.table("OBJECTS/Signature/PouchUCI_22_signatures.txt", T, '\t')
rownames(counts_ours) = counts_ours[,1]
counts_ours=counts_ours[,-1]
counts_ours = counts_ours[rownames(count_simp),]


cc = cor(counts_ours, count_simp)

our_t  = c("act_CD8_NR4A2_hi", 
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

their_t = c("CD4+ Activated Fos-hi","CD4+ Activated Fos-lo","CD4+ Memory" ,"CD4+ PD1+",
            "CD8+ IELs","CD8+ IL17+","CD8+ LP", "Cycling T", "NKs", "ILCs", "Tregs", "MT-hi")

our_b = c("cycling_b_cells", 
          "follicular_cells", 
          "gc_cells", 
          "plasma_NFKBIA_hi", 
          "plasma_NFKBIA_lo")

their_b = c("Cycling B", "Plasma", "Follicular", "GC")

our_m = c("mono_mac1", "mono_mac2", "mono_mac3", "mast_cells", "pdcs")

their_m =c("CD69- Mast","CD69+ Mast", "Cycling Monocytes", "DC1", "DC2", "Macrophages", "Inflammatory Monocytes")



heatty_t = grid.grabExpr(draw(
  Heatmap(cc[our_t, their_t], show_row_names = T, show_column_names = T,
          heatmap_legend_param = list(title = "Correlation\nValue"),
          cluster_rows = T, cluster_columns = T, row_names_side = 'right',
          #column_names_gp = gpar(fontsize=10),
          #row_names_gp = gpar(fontsize=12),
          column_names_rot = 60,
          #row_title_gp = gpar(fontsize = 10),
          #row_names_max_width = unit(10,'cm'),
          #row_split = sig_imm$cluster,
          use_raster = F,
          col = colorRamp2(c(0.2,0.35,0.5,0.65), c("grey95" ,"steelblue1", "steelblue3", "purple4"))),
  padding = unit(c(2, 2, 3, 2), "cm")
))

heatty_b = grid.grabExpr(draw(
  Heatmap(cc[our_b, their_b], show_row_names = T, show_column_names = T,
          heatmap_legend_param = list(title = "Correlation\nValue"),
          cluster_rows = T, cluster_columns = T, row_names_side = 'right',
          #column_names_gp = gpar(fontsize=10),
          #row_names_gp = gpar(fontsize=12),
          column_names_rot = 60,
          #row_title_gp = gpar(fontsize = 10),
          #row_names_max_width = unit(10,'cm'),
          #row_split = sig_imm$cluster,
          use_raster = F,
          col = colorRamp2(c(0.2,0.35,0.5,0.65), c("grey95" ,"steelblue1", "steelblue3", "purple4"))),
  padding = unit(c(9, 2, 3, 2), "cm")
))

heatty_m = grid.grabExpr(draw(
  Heatmap(cc[our_m, their_m], show_row_names = T, show_column_names = T,
          heatmap_legend_param = list(title = "Correlation\nValue"),
          cluster_rows = T, cluster_columns = T, row_names_side = 'right',
          #column_names_gp = gpar(fontsize=10),
          #row_names_gp = gpar(fontsize=12),
          column_names_rot = 60,
          #row_title_gp = gpar(fontsize = 10),
          #row_names_max_width = unit(10,'cm'),
          #row_split = sig_imm$cluster,
          use_raster = F,
          col = colorRamp2(c(0.2,0.35,0.5,0.65), c("grey95" ,"steelblue1", "steelblue3", "purple4"))),
  padding = unit(c(7, 2, 3, 2), "cm")
))

namer <- "FIGURES/FIGURE4/Figure4b_smillie_heat.pdf"
pdf(namer, height = 8, width = 20)
grid.arrange(heatty_t, heatty_m, heatty_b, nrow=1, widths = c(1.5,1,1))
dev.off()

#
#
#
## clean up metadata
ram_totals = data.frame(table(as.character(meta$Sample)))
colnames(ram_totals) = c("Sample", "cellTotal")
samps = as.character(unique(ram_totals$Sample))
health = c()
for(i in 1:length(samps)){
  curr = subset(meta, Sample == samps[i])
  new_h = as.character(curr$Health[1])
  health = c(health, new_h)
}
ram_totals$Health = health

meta$cellTotals = NA
for( i in 1:length(samps)){
  CT = as.numeric(subset(ram_totals, Sample == samps[i])$cellTotal)
  meta$cellTotals[meta$Sample==samps[i]] <- CT
}
#

#
#
#
#
# look for SOX4/MAFb, Foxp3 and CD8 cell states in different conditions
# IL1B APOE, LYZ, C1QA,B,C
#
#
#

#
#
### 
##### sox4/mafa cells
counts3 = scale(counts2)

sox4 = counts3["SOX4",]
mafb = counts3["MAFB",]
il1B = counts3["IL1B",]
lyz = counts3["LYZ",]
apoe = counts3["APOE",]
c1qc = counts3["C1QC",]
foxp3 = counts3["FOXP3",]
batf = counts3["BATF",]
cd8a = counts3["CD8A",]
cd69 = counts3["CD69",]
il2 = counts3["IL2",]

colorer = c("Healthy"="dodgerblue2", "Inflamed"="red3", "Non-inflamed"="purple", "other"="grey85")


mac_val = data.frame(cells = colnames(counts3),
                     SOX4 = sox4, MAFB=mafb,
                     BATF=batf, FOXP3=foxp3,
                     CD8A=cd8a, IL2=il2,
                     IL1B=il1B,LYZ =lyz,
                     APOE = apoe,C1QC=c1qc,
                     Health = meta$Health,
                     Cluster = meta$Cluster,
                     cellTotal = meta$cellTotals,
                     sample = meta$Sample)

mac_val$Health = factor(mac_val$Health, levels = c("Healthy", "Non-inflamed", "Inflamed", "other"))


### mono_mac 1
plotter = subset(mac_val, SOX4 > 0 & MAFB > 0)
plotter = subset(mac_val, grepl("Monocyte|Macropha", mac_val$Cluster))

plotter$Inflammation2 <- plotter$Health
plotter$Inflammation2[plotter$SOX4 < 0] <- "other"
plotter$Inflammation2[plotter$MAFB < 0] <- "other"

ggSOX = ggplot(plotter, aes(SOX4, MAFB, color=Inflammation2)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=colorer) +
  theme_bw() + facet_wrap(~Health, nrow=1) +
  geom_vline(xintercept =0) + geom_hline(yintercept = 0)


plotter_int <- subset(plotter, SOX4 > 0 & MAFB > 0)
samples <- unique(plotter_int$sample)
resser <- data.frame(sampleID = NA, Inflammation = NA, perc = NA)
for(j in 1:length(samples)){
  curr <- subset(plotter_int, sample == samples[j])
  perc = nrow(curr)/as.numeric(as.character(curr$cellTotal[1]))
  adder = data.frame(sampleID = samples[j],
                     Inflammation = curr$Health[1],
                     perc = perc)
  resser <- rbind(resser, adder)
}
resser=resser[-1,]
soxresser = resser

ressy=pairwise.wilcox.test(x=soxresser$perc, g=soxresser$Inflammation, p.adjust.method="fdr")$p.value
pvals=as.numeric(ressy)
comps=paste0(rep(rownames(ressy), ncol(ressy)), ";", rep(colnames(ressy), each=nrow(ressy)))
ressy1=data.frame(Comparison=comps, pvals=pvals)
ressy1$cluster="SOX4"

resser$Inflammation = factor(resser$Inflammation, levels = c("Inflamed","Non-inflamed", "Healthy", "other"))
bbSOX = ggplot(resser, aes(Inflammation, perc*100, color=Inflammation)) +
  geom_boxplot(outlier.shape=NA) +
  ggtitle("SOX4+ MAFB+ Cell Percentages") +
  ylab("Cell Percentages") +
  scale_color_manual(values = colorer) +
  geom_jitter(width=0.2) + theme_bw()


#
#
#
#
####mono_mac2

plotter = subset(mac_val, IL1B > 0 & LYZ > 0)
plotter = subset(mac_val, grepl("Monocyte|Macropha", mac_val$Cluster))

plotter$Inflammation2 <- plotter$Health
plotter$Inflammation2[plotter$IL1B < 0.5] <- "other"
plotter$Inflammation2[plotter$LYZ < 0.5] <- "other"

ggIL1B = ggplot(plotter, aes(IL1B, LYZ, color=Inflammation2)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=colorer) +
  theme_bw() + facet_wrap(~Health, nrow=1) +
  geom_vline(xintercept =0.5) + geom_hline(yintercept = 0.5)


plotter_int <- subset(plotter, IL1B > 0.5 & LYZ > 0.5)
samples <- as.character(unique(plotter_int$sample))
resser <- data.frame(sampleID = NA, Inflammation = NA, perc = NA)
for(j in 1:length(samples)){
  curr <- subset(plotter_int, sample == samples[j])
  perc = nrow(curr)/as.numeric(as.character(curr$cellTotal[1]))
  adder = data.frame(sampleID = samples[j],
                     Inflammation = curr$Health[1],
                     perc = perc)
  resser <- rbind(resser, adder)
}
resser=resser[-1,]
il1bresser = resser

ressy=pairwise.wilcox.test(x=il1bresser$perc, g=il1bresser$Inflammation, p.adjust.method="fdr")$p.value
pvals=as.numeric(ressy)
comps=paste0(rep(rownames(ressy), ncol(ressy)), ";", rep(colnames(ressy), each=nrow(ressy)))
ressy2=data.frame(Comparison=comps, pvals=pvals)
ressy2$cluster="IL1B"

resser$Inflammation = factor(resser$Inflammation, levels = c("Inflamed","Non-inflamed", "Healthy", "other"))
bbIL1B = ggplot(resser, aes(Inflammation, perc*100, color=Inflammation)) +
  geom_boxplot(outlier.shape=NA) +
  ggtitle("IL1B+ LYZ+ Cell Percentages") +
  ylab("Cell Percentages") +
  scale_color_manual(values = colorer) +
  geom_jitter(width=0.2) + theme_bw()

#

#
#
#
#
#
####mono_mac3

plotter = subset(mac_val, APOE > 0 & C1QC > 0)
plotter = subset(mac_val, grepl("Monocyte|Macropha", mac_val$Cluster))

plotter$Inflammation2 <- plotter$Health
plotter$Inflammation2[plotter$APOE < 0] <- "other"
plotter$Inflammation2[plotter$C1QC < 0] <- "other"

ggAPOE = ggplot(plotter, aes(APOE, C1QC, color=Inflammation2)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=colorer) +
  theme_bw() + facet_wrap(~Health, nrow=1) +
  geom_vline(xintercept =0) + geom_hline(yintercept = 0)


plotter_int <- subset(plotter, APOE > 0 & C1QC > 0)
samples <- unique(plotter_int$sample)
resser <- data.frame(sampleID = NA, Inflammation = NA, perc = NA)
for(j in 1:length(samples)){
  curr <- subset(plotter_int, sample == samples[j])
  perc = nrow(curr)/as.numeric(as.character(curr$cellTotal[1]))
  adder = data.frame(sampleID = samples[j],
                     Inflammation = curr$Health[1],
                     perc = perc)
  resser <- rbind(resser, adder)
}
resser=resser[-1,]
apoeresser = resser

ressy=pairwise.wilcox.test(x=apoeresser$perc, g=apoeresser$Inflammation, p.adjust.method="fdr")$p.value
pvals=as.numeric(ressy)
comps=paste0(rep(rownames(ressy), ncol(ressy)), ";", rep(colnames(ressy), each=nrow(ressy)))
ressy3=data.frame(Comparison=comps, pvals=pvals)
ressy3$cluster="APOE"

resser$Inflammation = factor(resser$Inflammation, levels = c("Inflamed","Non-inflamed", "Healthy", "other"))
bbAPOE = ggplot(resser, aes(Inflammation, perc*100, color=Inflammation)) +
  geom_boxplot(outlier.shape=NA) +
  ggtitle("APOE+ C1QC+ Cell Percentages") +
  ylab("Cell Percentages") +
  scale_color_manual(values = colorer) +
  geom_jitter(width=0.2) + theme_bw()


#
#
#
#
#
#
#

plotter = subset(mac_val, FOXP3 > 0 & BATF > 0)
plotter = subset(mac_val, grepl("CD8|CD4|Treg|Cycling T|NKs", mac_val$Cluster))

plotter$Inflammation2 <- plotter$Health
plotter$Inflammation2[plotter$FOXP3 < 0] <- "other"
plotter$Inflammation2[plotter$BATF < 0] <- "other"

ggFOXP3 = ggplot(plotter, aes(FOXP3, BATF, color=Inflammation2)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=colorer) +
  theme_bw() + facet_wrap(~Health, nrow=1) +
  geom_vline(xintercept =0) + geom_hline(yintercept = 0)


plotter_int <- subset(plotter, FOXP3 > 0 & BATF > 0)
samples <- unique(plotter_int$sample)
resser <- data.frame(sampleID = NA, Inflammation = NA, perc = NA)
for(j in 1:length(samples)){
  curr <- subset(plotter_int, sample == samples[j])
  perc = nrow(curr)/as.numeric(as.character(curr$cellTotal[1]))
  adder = data.frame(sampleID = samples[j],
                     Inflammation = curr$Health[1],
                     perc = perc)
  resser <- rbind(resser, adder)
}
resser=resser[-1,]
foxresser = resser

ressy=pairwise.wilcox.test(x=foxresser$perc, g=foxresser$Inflammation, p.adjust.method="fdr")$p.value
pvals=as.numeric(ressy)
comps=paste0(rep(rownames(ressy), ncol(ressy)), ";", rep(colnames(ressy), each=nrow(ressy)))
ressy4=data.frame(Comparison=comps, pvals=pvals)
ressy4$cluster="FOXP3"

resser$Inflammation = factor(resser$Inflammation, levels = c("Inflamed","Non-inflamed", "Healthy", "other"))
bbFOXP3 = ggplot(resser, aes(Inflammation, perc*100, color=Inflammation)) +
  geom_boxplot(outlier.shape=NA) +
  ggtitle("FOXP3+ BATF+ Cell Percentages") +
  ylab("Cell Percentages") +
  scale_color_manual(values = colorer) +
  geom_jitter(width=0.2) + theme_bw()


#
#
#
###
plotter = subset(mac_val, CD8A > 0 & IL2 > 0)
plotter = subset(mac_val, grepl("CD8|CD4|Treg|Cycling T|NKs", mac_val$Cluster))

plotter$Inflammation2 <- plotter$Health
plotter$Inflammation2[plotter$CD8A < 0] <- "other"
plotter$Inflammation2[plotter$IL2 < 0] <- "other"

ggCD8A = ggplot(plotter, aes(CD8A, IL2, color=Inflammation2)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=colorer) +
  theme_bw() + facet_wrap(~Health, nrow=1) +
  geom_vline(xintercept =0) + geom_hline(yintercept = 0)


plotter_int <- subset(plotter, CD8A > 0 & IL2 > 0)
samples <- unique(plotter_int$sample)
resser <- data.frame(sampleID = NA, Inflammation = NA, perc = NA)
for(j in 1:length(samples)){
  curr <- subset(plotter_int, sample == samples[j])
  perc = nrow(curr)/as.numeric(as.character(curr$cellTotal[1]))
  adder = data.frame(sampleID = samples[j],
                     Inflammation = curr$Health[1],
                     perc = perc)
  resser <- rbind(resser, adder)
}
resser=resser[-1,]
cd8resser = resser

ressy=pairwise.wilcox.test(x=cd8resser$perc, g=cd8resser$Inflammation, p.adjust.method="fdr")$p.value
pvals=as.numeric(ressy)
comps=paste0(rep(rownames(ressy), ncol(ressy)), ";", rep(colnames(ressy), each=nrow(ressy)))
ressy5=data.frame(Comparison=comps, pvals=pvals)
ressy5$cluster="CD8/IL2"

resser$Inflammation = factor(resser$Inflammation, levels = c("Inflamed","Non-inflamed", "Healthy", "other"))
bbCD8A = ggplot(resser, aes(Inflammation, perc*100, color=Inflammation)) +
  geom_boxplot(outlier.shape=NA) +
  ggtitle("CD8A+ IL2+ Cell Percentages") +
  ylab("Cell Percentages") +
  scale_color_manual(values = colorer) +
  geom_jitter(width=0.2) + theme_bw()


namer <- "FIGURES/FIGURE4/Figure4cd_smillie_pops.pdf"
pdf(namer, height = 5, width = 11)
grid.arrange(ggIL1B, bbIL1B,
             #ggAPOE, bbAPOE,
             ggFOXP3, bbFOXP3,
             nrow=2, widths = c(2,1.2))
dev.off()

namer <- "FIGURES/FIGURE4/Figure4cd_smillie_pops.png"
png(namer, height = 5, width = 11, units = "in", res=200)
grid.arrange(ggIL1B, bbIL1B,
             #ggAPOE, bbAPOE,
             ggFOXP3, bbFOXP3,
             nrow=2, widths = c(2,1.2))
dev.off()

namer <- "FIGURES/FIGURE4/Figure4cd_smillie_pops_simp.pdf"
pdf(namer, height = 5, width = 4)
grid.arrange(bbIL1B,
             bbFOXP3,
             nrow=2)
dev.off()



####
namer <- "FIGURES/SF6/FigureSF6_smillie_othr_pops.pdf"
pdf(namer, height = 5, width = 11)
grid.arrange(ggSOX, bbSOX,
             ggAPOE, bbAPOE,
             nrow=2, widths = c(2,1.2))
dev.off()

namer <- "FIGURES/SF6/FigureSF6_smillie_othr_pops.png"
png(namer, height = 5, width = 11, units = "in", res=200)
grid.arrange(ggSOX, bbSOX,
             ggAPOE, bbAPOE,
             nrow=2, widths = c(2,1.2))
dev.off()

ressy = rbind(ressy1, ressy2, ressy3, ressy4, ressy5)
write.table(ressy, "FIGURES/FIGURE4/Figure4_facs_stats_smillie.txt", sep='\t', quote=F)
