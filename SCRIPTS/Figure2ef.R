library(slingshot)
library(Seurat)
library(destiny)
library(gam)
library(ggplot2)
library(SingleCellExperiment)
library(RColorBrewer)
library(cowplot)
library(gridExtra)

colorer_pouch <- c("Normal_pouch"="dodgerblue2","Pouchitis"="red3", "UC_colon"="forestgreen", "UC_inflamed"="darkorange1")

colorer_cells <- c("mast_cells"="navy", 
                   "pdcs"="pink", "DC"="darkorange1",
                   "mono_mac1"="red3", "mono_mac2"="purple", "mono_mac3"="forestgreen")

s_obj <- readRDS("OBJECTS/Myeloid_cells/seurat_obj.rds")
s_obj@meta.data$MinorPopulations2 <- factor(s_obj@meta.data$MinorPopulations2, 
                                           levels = c("mono_mac1", "mono_mac2", "mono_mac3", "mast_cells","DC", "pdcs"))
Idents(s_obj) = s_obj@meta.data$MinorPopulations2
s_obj <- subset(s_obj, MinorPopulations2 %in% c("mono_mac1", "mono_mac2", "mono_mac3"))


counts <- counts <- s_obj@assays$RNA@counts
sim <- SingleCellExperiment(assays=List(counts=counts))

meta <- s_obj@meta.data

#filter genes

geneFilter <- apply(assays(sim)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sim <- sim[geneFilter, ]


##
## quant norm

FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sim)$norm <- FQnorm(assays(sim)$counts)


#
#

### diffusion map
dm <- DiffusionMap(t(log1p(assays(sim)$norm)))
rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)

#add reductions
reducedDims(sim) <- SimpleList(DiffMap = rd2)


## add clusters
colData(sim)$cell_states <- s_obj@meta.data$MinorPopulations2
colData(sim)$orig.ident <- s_obj@meta.data$pouch.status

##

###perform slingshot
### diff map


sim_diff <- slingshot(sim, clusterLabels = 'cell_states', reducedDim = 'DiffMap')

plotter <- data.frame(sim_diff@reducedDims$DiffMap, 
                      cell_states=colData(sim_diff)$cell_states, 
                      pouch.status=colData(sim_diff)$orig.ident,
                      psuedotime=sim_diff$slingPseudotime_1)
pc = SlingshotDataSet(sim_diff)@curves$curve1
fittedLine <- data.frame(pc$s[pc$ord, ])



p1 = ggplot(sample(plotter), aes(DC1, DC2, color=cell_states)) + geom_point() +
  geom_path(data = fittedLine, col = 'black') +
  scale_color_manual(values=colorer_cells) + theme_bw()

p2 = ggplot(sample(plotter), aes(DC1, DC2, color=pouch.status)) + geom_point() +
  geom_path(data = fittedLine, col = 'black') +
  scale_color_manual(values=colorer_pouch) + theme_bw()

ggsave(paste0("FIGURES/FIGURE2/Figure2e_psuedotime_cells.png"), plot = p1, width = 7, height = 5, units = "in")
ggsave(paste0("OBJECTS/Myeloid_cells/slingshot/diffMap_pouch.png"), plot = p2, width = 7, height = 5, units = "in")

#
#
#### genes change with pseudotime

t <- sim_diff$slingPseudotime_1
Y <- log1p(assays(sim_diff)$norm)

var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:round(quantile(1:nrow(Y),0.1),0)]
Y <- Y[var100,]

# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})

#
m1 = (rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 == "mono_mac1"]])-rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 != "mono_mac1"]]))/rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 == "mono_mac1"]])
m2 = (rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 == "mono_mac2"]])-rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 != "mono_mac2"]]))/rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 == "mono_mac2"]])
m3 = (rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 == "mono_mac3"]])-rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 != "mono_mac3"]]))/rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 == "mono_mac3"]])
m2m3 = (rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 != "mono_mac1"]])-rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 == "mono_mac1"]]))/rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 != "mono_mac1"]])

m1 = rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 == "mono_mac1"]])
m2 = rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 == "mono_mac2"]])
m3 = rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 == "mono_mac3"]])
m2m3 = rowMeans(Y[,rownames(meta)[meta$MinorPopulations2 != "mono_mac1"]])

topgenes <- data.frame(gene=rownames(Y), m1=m1,m2=m2,m3=m3,m2m3=m2m3, gam.pval=gam.pval)
topgenes <- subset(topgenes, gam.pval < 0.05 & !grepl("MT-|HSP", rownames(topgenes)))
topgenes23 <- topgenes[order(topgenes$m2m3, decreasing=T),]
topgenes1 <- topgenes[order(topgenes$m1, decreasing=T),]
topgenes2 <- topgenes[order(topgenes$m2, decreasing=T),]

write.table(rbind(topgenes1), "OBJECTS/Myeloid_cells/slingshot/topgenes.txt", sep='\t', row.names=F, quote=F)

best_genes <- c(as.character(topgenes1$gene[1:20]), 
                as.character(topgenes2$gene[1:20]),
                as.character(topgenes23$gene[1:20]))

best_genes = unique(best_genes)[c(1:16, 18:21)]

plist=list()
pouchlist=list()
cellslist=list()
for(i in 1:length(best_genes)){
  plotter$gene_color <- Y[best_genes[i],]
  
  p2 = ggplot(sample(plotter), aes(DC1, DC2, color=gene_color)) + geom_point() +
    geom_path(data = fittedLine, col = 'black') +
    ggtitle(best_genes[i]) +
    scale_colour_gradient(low="grey85", high="purple4") + theme_bw()
  
  plist[[i]] <- arrangeGrob(p2)
  
  ppouch = ggplot(sample(plotter), aes(psuedotime, y=gene_color, color=pouch.status)) + geom_point() +
    geom_smooth(col='black') +
    ggtitle(best_genes[i]) +
    scale_colour_manual(values=colorer_pouch) + theme_bw()
  
  pcells = ggplot(sample(plotter), aes(psuedotime, y=gene_color, color=cell_states)) + geom_point() +
    geom_smooth(col='black') +
    ggtitle(best_genes[i]) +
    scale_colour_manual(values=colorer_cells) + theme_bw()
  
  plist[[i]] <- arrangeGrob(p2)
  pouchlist[[i]] <- arrangeGrob(ppouch)
  cellslist[[i]] <- arrangeGrob(pcells)
  
}

png("OBJECTS/Myeloid_cells/slingshot/DiffMap_topGenes_psuedotime2.png", height =20, width = 30, units='in', res=100)
plot_grid(plotlist=plist, ncol = 5)
dev.off()

png("OBJECTS/Myeloid_cells/slingshot/DiffMap_topGenes_psuedotime_pouch2.png", height =20, width = 30, units='in', res=100)
plot_grid(plotlist=pouchlist, ncol = 5)
dev.off()

png("OBJECTS/Myeloid_cells/slingshot/DiffMap_topGenes_psuedotime_cells2.png", height =20, width = 30, units='in', res=100)
plot_grid(plotlist=cellslist, ncol = 5)
dev.off()

saveRDS(sim_diff, "OBJECTS/Myeloid_cells/slingshot/slingshot_obj.rds")

### Figure 2f
# plotter$gene_color <- Y["VAMP8",]
# vamp8 = ggplot(sample(plotter), aes(DC1, DC2, color=gene_color)) + geom_point() +
#   geom_path(data = fittedLine, col = 'black') +
#   ggtitle("VAMP8") +
#   scale_colour_gradient(low="grey85", high="purple4") + theme_bw()
# 
# plotter$gene_color <- Y["CLEC10A",]
# clec10a = ggplot(sample(plotter), aes(DC1, DC2, color=gene_color)) + geom_point() +
#   geom_path(data = fittedLine, col = 'black') +
#   ggtitle("CLEC10A") +
#   scale_colour_gradient(low="grey85", high="purple4") + theme_bw()
# 
# plotter$gene_color <- Y["CXCL10",]
# sox4 = ggplot(sample(plotter), aes(DC1, DC2, color=gene_color)) + geom_point() +
#   geom_path(data = fittedLine, col = 'black') +
#   ggtitle("CXCL10") +
#   scale_colour_gradient(low="grey85", high="purple4") + theme_bw()
# 
# pdf("FIGURES/FIGURE2/Figure2f_psuedotime_genes.pdf", height = 4, width = 15)
# grid.arrange(vamp8, clec10a, sox4, nrow=1)
# dev.off()



##
#
#
### Figure 2f
plotter$gene_color <- Y["FOS",]
vamp8 = ggplot(sample(plotter), aes(DC1, DC2, color=gene_color)) + geom_point() +
  geom_path(data = fittedLine, col = 'black') +
  ggtitle("FOS") +
  scale_colour_gradient(low="grey85", high="purple4") + theme_bw()

plotter$gene_color <- Y["JUN",]
clec10a = ggplot(sample(plotter), aes(DC1, DC2, color=gene_color)) + geom_point() +
  geom_path(data = fittedLine, col = 'black') +
  ggtitle("JUN") +
  scale_colour_gradient(low="grey85", high="purple4") + theme_bw()

plotter$gene_color <- Y["LYZ",]
sox4 = ggplot(sample(plotter), aes(DC1, DC2, color=gene_color)) + geom_point() +
  geom_path(data = fittedLine, col = 'black') +
  ggtitle("LYZ") +
  scale_colour_gradient(low="grey85", high="purple4") + theme_bw()

pdf("FIGURES/FIGURE2/Figure2f_psuedotime_genes2.pdf", height = 4, width = 15)
grid.arrange(vamp8, clec10a, sox4, nrow=1)
dev.off()

