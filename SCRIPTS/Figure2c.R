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

s_obj <- readRDS("OBJECTS/Myeloid_cells/seurat_obj.rds")
s_obj@meta.data$MinorPopulations2 <- factor(s_obj@meta.data$MinorPopulations2, 
                                           levels = c("mono_mac1", "mono_mac2", "mono_mac3", 
                                                      "mast_cells", "DC", "pdcs"))

cluster = s_obj@meta.data[,"MinorPopulations2"]

Idents(s_obj) <- s_obj@meta.data$pouch.status

barplot = data.frame(PouchStatus=Idents(s_obj), Cluster = cluster)
barplot$Sample <- s_obj@meta.data$orig.ident
bb=ggplot(barplot, aes(PouchStatus, fill=cluster)) + geom_bar() + scale_fill_manual(values=colors_clusters) + theme_bw() + ylab("Cell count")
bb2=ggplot(barplot, aes(PouchStatus, fill=cluster)) + geom_bar(position="fill") + scale_fill_manual(values=colors_clusters) + theme_bw() + ylab("Proportion")

df = table(barplot$Cluster, barplot$Sample, barplot$PouchStatus)
dfm = melt(df)
dfm = subset(dfm, value > 0)
colnames(dfm) <- c("Cluster", "Sample", "PouchStatus", "value")
dfm$Cluster <- factor(dfm$Cluster, levels = c("mono_mac1", "mono_mac2", "mono_mac3", "mast_cells", "DC", "pdcs"))
dfm$PouchStatus <- factor(dfm$PouchStatus, levels = c("Normal_pouch", "Pouchitis", "UC_inflamed"))
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
write_to <- paste0("OBJECTS/Myeloid_cells/clusters-resolutions/clust_box_MinorPopulations2_values.txt")
write.table(dfm, write_to, sep='\t', row.names=F, quote=F)

dfm$Cluster <- factor(dfm$Cluster, levels = c("mono_mac1", "mono_mac2", "mono_mac3", "mast_cells","DC", "pdcs"))

melt_plot2 <- ggplot(dfm, aes(PouchStatus, perc, color=PouchStatus)) + 
geom_boxplot(outlier.shape=NA) + geom_jitter(width = 0.2) + xlab("PouchStatus") +
theme_bw() + scale_color_manual(values = colorer) + ylab("Cell Percentage") + 
  theme(axis.text.x=element_blank()) +
  facet_wrap(~Cluster, nrow=1, scales="free_y")

bb3 = ggplot(barplot, aes(cluster, fill=PouchStatus)) + geom_bar() + scale_fill_manual(values=colorer) + theme_bw() + ylab("Cell Count")
bb4 = ggplot(barplot, aes(cluster, fill=PouchStatus)) + geom_bar(position="fill") + scale_fill_manual(values=colorer) + theme_bw() + ylab("Proportion")
bb = arrangeGrob(bb,bb2,bb3,bb4, ncol=2)

#write_to <- list.files(args[1])[grep(gsub("\\.","", args[2]), list.files(args[1]))]
write_to <- paste0("OBJECTS/Myeloid_cells/clusters-resolutions/clust_bar_MinorPopulations2_values.pdf")

ggsave(write_to, plot = bb, width = 10, height = 10, units = "in")

mm = arrangeGrob(melt_plot1, melt_plot2, ncol=1)
write_to <- paste0("FIGURES/FIGURE2/Figure2C_MinorPopulations2_boxplot.pdf")
ggsave(write_to, plot = mm, width = 12, height = 6, units = "in")

