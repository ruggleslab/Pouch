library(Seurat)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
library(RColorBrewer)
library(ggsci)
library(gridExtra)
library(reshape)

colors_samples = c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), pal_igv("default")(51))
colors_clusters = c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))

colorer <- c("Normal_pouch"="dodgerblue2","Pouchitis"="red3", "UC_colon"="forestgreen", "UC_inflamed"="darkorange1")

s_obj = readRDS("OBJECTS/seurat_obj.rds")
cluster = s_obj@meta.data[,"MajorPopulations"]

Idents(s_obj) <- s_obj@meta.data$pouch.status

barplot = data.frame(PouchStatus=Idents(s_obj), Cluster = cluster)
barplot$Sample <- s_obj@meta.data$orig.ident
bb=ggplot(barplot, aes(PouchStatus, fill=cluster)) + geom_bar() + scale_fill_manual(values=colors_clusters) + theme_bw() + ylab("Cell count")
bb2=ggplot(barplot, aes(PouchStatus, fill=cluster)) + geom_bar(position="fill") + scale_fill_manual(values=colors_clusters) + theme_bw() + ylab("Proportion")

df = table(barplot$Cluster, barplot$Sample, barplot$PouchStatus)
dfm = melt(df)
dfm = subset(dfm, value > 0)
colnames(dfm) <- c("Cluster", "Sample", "PouchStatus", "value")
dfm$Cluster <- factor(dfm$Cluster, levels = levels(s_obj@meta.data$MajorPopulations))
melt_plot1 <- ggplot(dfm, aes(PouchStatus, value, color=PouchStatus)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(width = 0.2) + xlab("PouchStatus") +
  theme_bw() + scale_color_manual(values = colorer) + ylab("Cell Count") + 
  theme(axis.text.x=element_blank()) +
  facet_wrap(~Cluster, nrow=1, scales="free_y")

dfm$totals <- 0
samps <- as.character(unique(dfm$Sample))
for(i in 1:length(samps)){
  curr = subset(s_obj@meta.data, orig.ident == samps[i])
  totaler = as.numeric(curr$cell_total[1])
  dfm$totals[dfm$Sample == samps[i]] <- totaler
}

dfm$totals <- as.numeric(dfm$totals)
dfm$perc <- (dfm$value/dfm$totals)*100
write_to <- paste0("OBJECTS/clusters-resolutions/clust_box_MajorPopulations_values.txt")
write.table(dfm, write_to, sep='\t', row.names=F, quote=F)

dfm$Cluster <- factor(dfm$Cluster, levels = levels(s_obj@meta.data$MajorPopulations))

melt_plot2 <- ggplot(dfm, aes(PouchStatus, perc, color=PouchStatus)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(width = 0.2) + xlab("PouchStatus") +
  theme_bw() + scale_color_manual(values = colorer) + ylab("Cell Percentage") + 
  theme(axis.text.x=element_blank()) +
  facet_wrap(~Cluster, nrow=1, scales="free_y")

bb3 = ggplot(barplot, aes(cluster, fill=PouchStatus)) + geom_bar() + scale_fill_manual(values=colorer) + theme_bw() + ylab("Cell Count")
bb4 = ggplot(barplot, aes(cluster, fill=PouchStatus)) + geom_bar(position="fill") + scale_fill_manual(values=colorer) + theme_bw() + ylab("Proportion")
bb = arrangeGrob(bb,bb2,bb3,bb4, ncol=2)

#write_to <- list.files(args[1])[grep(gsub("\\.","", args[2]), list.files(args[1]))]
write_to <- paste0("OBJECTS/clusters-resolutions/clust_bar_MajorPopulations_values.pdf")

ggsave(write_to, plot = bb, width = 10, height = 10, units = "in")

mm = arrangeGrob(melt_plot1, melt_plot2, ncol=1)
write_to <- paste0("FIGURES/FIGURE1/Figure1D_MajorPopulations_boxplot.pdf")
ggsave(write_to, plot = mm, width = 15, height = 6, units = "in")

#
#
#
##
#
#
#### specific UMAPs

source("SCRIPTS/gene_umap.R")

markers = c("CD3D", "BANK1", "STMN1", "MZB1", "IL1B", "KIT")
uplots = arrangeGrob(gene_umap(s_obj, markers)+theme_void())

ggsave("FIGURES/FIGURE1/Figure1E_UMAPS.png", plot = uplots, width = 35, height = 6, units = "in")


#
#
#
#
#
#
#
#
#


###### PCA and barplots
meta = s_obj@meta.data
Freq <- data.frame(table(meta$orig.ident, meta$MajorPopulations))
Freq <- cast(Freq, Var1 ~ Var2)
rownames(Freq) <- Freq[,1]
Freq=Freq[,-1]
Freq=data.matrix(Freq)
totals <- data.frame(samples=unique(meta$orig.ident), totals=as.numeric(unique(meta$cell_total)))
rownames(totals) <- totals$samples
totals=totals[rownames(Freq),]

for(i in 1:ncol(Freq)){
  Freq[,i] <- Freq[,i]/totals$totals
}
#

pc = prcomp(Freq)
why=pc$rotation
rownames(why) <- colnames(Freq)

PoV <- pc$sdev^2/sum(pc$sdev^2)

meta2 = data.frame(table(meta$orig.ident, meta$pouch.status))
meta2 = subset(meta2, Freq > 0)
colnames(meta2) <- c("Sample.ID", "pouch.status", "cells")
rownames(meta2) = meta2[,1]
meta2 = meta2[rownames(Freq),]

pc_plot <- data.frame(meta2, PC1=pc$x[,1], PC2=pc$x[,2])

pp1=ggplot(pc_plot, aes(PC1, PC2, color=pouch.status)) +
  geom_point(size=4) + scale_color_manual(values = colorer) +
  theme_bw() +
  xlab(paste0("PC1 (Explained Variance ", round(PoV[1],4)*100, "%)")) +
  ylab(paste0("PC2 (Explained Variance ", round(PoV[2],4)*100, "%)")) 



#
#
###
colors_clusters <- c("tcells"="forestgreen", 
             "bcells_gc_f"="darkred", "bcells_cycling"="orangered", "bcells_plasma"="lightcoral", 
             "mcells_monomac"="navy", "mcells_mast"="deepskyblue4", "mcells_pdcs"="dodgerblue2")

barplot = data.frame(PouchStatus=meta$pouch.status, Cluster = meta$MajorPopulations)
barplot$Sample <- s_obj@meta.data$orig.ident

Freq2 <- Freq[order(Freq[,1], Freq[,2], Freq[,4], Freq[,5], Freq[,3], Freq[,6], decreasing=T),]
barplot$Sample <- factor(barplot$Sample, levels = rownames(Freq2))

barplot$Cluster <- factor(barplot$Cluster, levels = levels(s_obj@meta.data$MajorPopulations))


uci=ggplot(subset(barplot, PouchStatus=="UC_inflamed"), aes(Sample, fill=Cluster)) + 
  geom_bar(position="fill") + ggtitle("UC Inflamed") +
  scale_fill_manual(values=colors_clusters) + theme_bw() + ylab("Cell Proportion") +
  theme(legend.position='none',
        axis.text.x=element_blank())
npouch=ggplot(subset(barplot, PouchStatus=="Normal_pouch"), aes(Sample, fill=Cluster)) + 
  geom_bar(position="fill") + ggtitle("Normal Pouch") +
  scale_fill_manual(values=colors_clusters) + theme_bw() + ylab("Cell Proportion") +
  theme(legend.position='none',
        axis.text.x=element_blank())
pouchitis=ggplot(subset(barplot, PouchStatus=="Pouchitis"), aes(Sample, fill=Cluster)) + 
  geom_bar(position="fill") + ggtitle("Pouchitis") +
  scale_fill_manual(values=colors_clusters) + theme_bw() + ylab("Cell Proportion") +
  theme(axis.text.x=element_blank())

bb = arrangeGrob(uci,npouch,pouchitis, nrow=1, widths = c(1,1,1.4))


write_to <- paste0("FIGURES/SF2/FigureSF1C.pdf")
pdf(write_to, height = 5, width = 17)
grid.arrange(pp1, bb, nrow=1, widths = c(1,1.6))
dev.off()
