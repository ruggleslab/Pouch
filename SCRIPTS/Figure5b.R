library(limma)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(broom)
library(matrixStats)

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

disease_cols=c("NP"="grey85", "FAP"="purple", "PI"="red3")

table = read.table("additional_data/Morgan_et_al/ipaa_counts.txt", T, '\t')
rownames(table) <- table[,1]
table <- table[,-1]

meta = read.table("additional_data/Morgan_et_al/ipaa_meta.txt", T, '\t')
rownames(meta) = meta$Sample_geo_accession
table = table[,rownames(meta)]

ct_ref = read.table("OBJECTS/Signature/signatureRefCt.txt", T, "\t")
rownames(ct_ref) <- ct_ref$gene
ct_ref = ct_ref[rownames(table),]

group = gsub("-", "_", meta$Outcome)

group2 = gsub("CDL", "PI", group)
group2 = gsub("AP", "PI", group2)
group2 = gsub("CP", "PI", group2)
group2 = gsub("FPI", "FAP", group2)

group2 = factor(group2, levels = c("NP", "FAP", "PI"))

pp = prcomp(t(table))

PoV <- pp$sdev^2/sum(pp$sdev^2)

meta$PC1 = pp$x[,1]
meta$PC2 = pp$x[,2]
g1=ggplot(meta, aes(PC1, PC2, color=group2))+geom_point() + theme_bw() +
  scale_color_manual(values=disease_cols) +
  xlab(paste0("PC1 (Explained Variance ", round(PoV[1],4)*100, "%)")) +
  ylab(paste0("PC2 (Explained Variance ", round(PoV[2],4)*100, "%)")) 
g2=ggplot(meta, aes(PC1, PC2, color=ISCORE))+geom_point() + theme_bw() +
  scale_color_gradient(low='grey85', high="red3") +
  xlab(paste0("PC1 (Explained Variance ", round(PoV[1],4)*100, "%)")) +
  ylab(paste0("PC2 (Explained Variance ", round(PoV[2],4)*100, "%)")) 

pdf("FIGURES/FIGURE5/Figure5b_ipaa_pca.pdf", height =4, width = 10)
grid.arrange(g1,g2, nrow=1)
dev.off()

#
#
#
#
#
#####
###
#

ggplot(meta, aes(group2,ISCORE, color=group2))+
  #geom_boxplot() +
  geom_jitter(width = 0.2) + theme_bw() +
  scale_color_manual(values=disease_cols)
##


##
#
#
#
# #mono_mac2_g = as.character(ct_ref$gene[ct_ref$cluster=="mono_mac2"])
# #counts = table[mono_mac2_g,]
# counts = table
# 
# ha = rowAnnotation(df=data.frame(cellState=ct_ref$cluster),
#                    col = list(cellState=big_colorer))
# 
# ha2 = columnAnnotation(df=data.frame(ISCORE=meta$ISCORE, Disease=group2),
#                        col = list(Disease=c("NP"="grey85", "FAP"="purple", "PI"="red3")))
# 
# ct_ref$cluster = factor(ct_ref$cluster, levels = unique(ct_ref$cluster))
# 
# heatty = grid.grabExpr(draw(
#   Heatmap(t(scale(t(counts))), show_row_names = T, show_column_names = F,
#           heatmap_legend_param = list(title = "Scaled\nValue"),
#           top_annotation = ha2,
#           cluster_rows = F, cluster_columns = T, row_names_side = 'left',
#           column_names_gp = gpar(fontsize=10),
#           #row_names_gp = gpar(fontsize=12),
#           column_names_rot = 60,
#           #row_title_gp = gpar(fontsize = 10),
#           #row_names_max_width = unit(10,'cm'),
#           #column_split = meta$Remission,
#           row_split = ct_ref$cluster,
#           cluster_row_slices = F,
#           use_raster = F,
#           col = colorRamp2(c(-2,0,2), c("blue","white" ,"red")))+ha,
#   padding = unit(c(2, 2, 3, 2), "cm")
# ))
# 
# namer <- "FIGURES/SF6/FigureSF6a_ipaa_details_large.pdf"
# pdf(namer, height = 40, width = 15)
# grid.arrange(heatty)
# #grid.arrange(ipaa_p, heatty, nrow=2, heights = c(1.3,1))
# dev.off()
# 
# 
# ##
# heatty = grid.grabExpr(draw(
#   Heatmap(t(scale(t(counts))), show_row_names = F, show_column_names = F,
#           heatmap_legend_param = list(title = "Scaled\nValue"),
#           top_annotation = ha2,
#           cluster_rows = F, cluster_columns = T, row_names_side = 'left',
#           column_names_gp = gpar(fontsize=10),
#           #row_names_gp = gpar(fontsize=12),
#           column_names_rot = 60,
#           #row_title_gp = gpar(fontsize = 10),
#           #row_names_max_width = unit(10,'cm'),
#           #column_split = meta$Remission,
#           row_split = ct_ref$cluster,
#           cluster_row_slices = F,
#           use_raster = F,
#           col = colorRamp2(c(-2,0,2), c("blue","white" ,"red")))+ha,
#   padding = unit(c(2, 2, 3, 2), "cm")
# ))
# 
# namer <- "FIGURES/SF6/FigureSF6a_ipaa_details_small.pdf"
# pdf(namer, height = 8, width = 10)
# grid.arrange(heatty)
# #grid.arrange(ipaa_p, heatty, nrow=2, heights = c(1.3,1))
# dev.off()



#
#
#
#
#
#
#
#
#
#
#### repeat limma!!!
group3 = rep("Uninflamed", nrow(meta))
group3[meta$ISCORE>2] <- "Inflamed"

fitter = data.frame(logFC=NA,AveExpr=NA, t=NA, P.Value=NA, adj.P.Val=NA, B=NA, Gene=NA, comp=NA)
  design <- model.matrix(~0 +group3)
  fit <- lmFit(table, design)
  cc = paste0("group3Inflamed", "-group3Uninflamed")
  contrast.matrix <- makeContrasts(cc, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit2)
  fit_filt <- subset(topTable(fit, n=Inf))
  #fit_filt <- fit_filt[order(fit_filt$logFC),]
  fit_filt$Gene = rownames(fit_filt)
  fit_filt$comp = "Inflamed"
  fit_filt=fit_filt[rownames(table),]
  fitter = rbind(fitter, fit_filt)
fitter = fitter[-1,]


ct_ref = read.table("OBJECTS/Signature/signatureRefCt.txt", T, "\t")
rownames(ct_ref) <- ct_ref$gene
ct_ref = ct_ref[rownames(table),]

fitter$cellState = ct_ref$cluster

fitter = fitter[order(abs(fitter$logFC), decreasing=T),]

res_plot = fitter[1:30,]
res_plot = res_plot[order(res_plot$logFC, decreasing=T),]

res_plot$Gene = factor(res_plot$Gene, levels = as.character(unique(res_plot$Gene)))
ipaa_p=ggplot(res_plot, aes(x=Gene, y = logFC, fill=cellState)) + 
  geom_col() + coord_flip() + ggtitle("Inflamed vs. Uninflamed") +
  scale_fill_manual(values = big_colorer) + theme_bw() + #facet_wrap(~comp)
  theme(axis.text = element_text(color='black'))

##
il1b_iscore <- data.frame(il1b=as.numeric(table["IL1B",]), 
                          ISCORE=as.numeric(meta$ISCORE),
                          Gender=meta$Gender, Drug=meta$Antibiotics)
pval = round(glance(lm(il1b~ISCORE+Gender+Drug, il1b_iscore))$p.value,10)

## add covariates
lm_fit <- lm(il1b~ISCORE+Gender+Drug, il1b_iscore)
summary(lm_fit)

predicted_df <- data.frame(mpg_pred = predict(lm_fit, il1b_iscore), ISCORE=il1b_iscore$ISCORE)
  
###

comp_p = ggplot(meta, aes(as.numeric(table["IL1B",]), ISCORE, color=group3)) + 
  geom_point() + theme_bw() + xlab("IL1B Expression") +
  geom_smooth(data = predicted_df, aes(x=mpg_pred, y=ISCORE),
              method = lm, se = FALSE, color="black", linetype="dashed") +
  annotate("text", x = 11.2, y = 10, label = paste0("p < 2.2e-16"), color="black") +
  scale_color_manual(values = c("Inflamed"="tomato", "Uninflamed"="steelblue1")) +
  theme(axis.text = element_text(color='black'),
        legend.position = c(0.15, 0.87))

pdf("FIGURES/FIGURE5/Figure5c_ipaa_il1b_comp.pdf", height =4, width = 10)
grid.arrange(ipaa_p,comp_p, nrow=1, widths = c(0.55, 0.45))
dev.off()

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#### if multiple comps
multi_group = make.names(paste0(group2, ":", group3))

contrasts = c("PI.Inflamed", "PI.Uninflamed", "FAP.Uninflamed")
fitter = data.frame(logFC=NA,AveExpr=NA, t=NA, P.Value=NA, adj.P.Val=NA, B=NA, Gene=NA, comp=NA)
for(i in 1:length(contrasts)){
  design <- model.matrix(~0 +multi_group)
  fit <- lmFit(table, design)
  cc = paste0("multi_groupNP.Uninflamed", "-multi_group", contrasts[i])
  contrast.matrix <- makeContrasts(cc, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit2)
  fit_filt <- subset(topTable(fit, n=Inf))
  #fit_filt <- fit_filt[order(fit_filt$logFC),]
  fit_filt$Gene = rownames(fit_filt)
  fit_filt$comp = contrasts[i]
  fit_filt=fit_filt[rownames(table),]
  fitter = rbind(fitter, fit_filt)
}
fitter = fitter[-1,]

fitter$logFC=fitter$logFC*-1

ct_ref = read.table("OBJECTS/Signature/signatureRefCt.txt", T, "\t")
rownames(ct_ref) <- ct_ref$gene
ct_ref = ct_ref[rownames(table),]

fitter$cellState = rep(ct_ref$cluster,3)

fitter = fitter[order(abs(fitter$logFC), decreasing=T),]

#good_genes = unique(fitter[1:30,]$Gene)
good_genes = c(
  subset(fitter, comp == "PI.Inflamed")$Gene[1:10],
  subset(fitter, comp == "PI.Uninflamed")$Gene[1:10],
  subset(fitter, comp == "FAP.Uninflamed")$Gene[1:10]
)

res_plot=subset(fitter, Gene %in% good_genes)
res_plot = res_plot[order(res_plot$logFC, decreasing=T),]

res_plot$comp = factor(res_plot$comp, levels = c("PI.Inflamed", "PI.Uninflamed", "FAP.Uninflamed"))

res_plot$Gene = factor(res_plot$Gene, levels = as.character(unique(res_plot$Gene)))
ipaa_p=ggplot(res_plot, aes(x=Gene, y = logFC, fill=cellState)) + 
  geom_col() + coord_flip() + ggtitle("FAP and Pouchitis Immune Cell Differences") +
  scale_fill_manual(values = big_colorer) + theme_bw() + #facet_wrap(~comp)
  theme(axis.text = element_text(color='black'))+facet_wrap(~comp, nrow=1)

pdf("FIGURES/FIGURE5/Figure5c_ipaa_FAP_comp.pdf", height =4, width = 10)
ipaa_p
dev.off()
#
#

#
#
#
##

#### association with microial taxa?

library(mixOmics)

##read in and organize counts, otus and meta
otus = read.table("additional_data/Morgan_et_al/ipaa_otu_table.txt", T, '\t')
counts = read.table("additional_data/Morgan_et_al/ipaa_counts.txt", T, '\t')
meta = read.table("additional_data/Morgan_et_al/ipaa_meta.txt", T, '\t')
meta$group3 = group3

otu_names = otus[,c(1,ncol(otus))]
rownames(otu_names) <- otu_names[,1]
rownames(otus) = otus[,1]
otus = otus[,-1]

write.table(otu_names, "additional_data/Morgan_et_al/ipaa_otu_name_ref.txt", row.names=F, sep='\t', quote=F)

overlapper_MID = intersect(meta$MID, colnames(otus))
otus = otus[,overlapper_MID]

rownames(meta) = meta$MID
meta_otu = meta[overlapper_MID,]

colnames(otus) = as.character(meta_otu$Sample_geo_accession)
rownames(meta_otu) = as.character(meta_otu$Sample_geo_accession)

counts_otus = counts[,rownames(meta_otu)]
rownames(counts_otus) = counts$GeneSymbol
otus = otus[,rownames(meta_otu)]

ct_ref = read.table("OBJECTS/Signature/signatureRefCt.txt", T, "\t")
rownames(ct_ref) <- ct_ref$gene
ct_ref = ct_ref[rownames(counts_otus),]
###norm

otus.raw =t(otus) + 1

low.count.removal = function(
  data, # OTU count data frame of size n (sample) x p (OTU)
  percent=0.01 # cutoff chosen
){
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

result.filter = low.count.removal(otus.raw, percent=0.01)
otus_norm = result.filter$data.filter

TSS.divide = function(x){
  x/sum(x)
}

norm.TSS = t(apply(otus_norm, 1, TSS.divide))
otu_names.norm = otu_names[colnames(norm.TSS),]

Y  =  as.factor(meta_otu$group3)
sample = as.factor(meta_otu$Sample_geo_accession)

#
#
#
###### sPLS-DA
## Tuning sPLS-DA

#The tuning of sPLS-DA is being performed one component at a time inside the function and the optimal number of variables to select is automatically retrieved for each component. We set ncomp = 5 and we used 10-fold cross validation (folds = 5 repeated 10 times).

# set the tuning grid for some possible values of keepX
list.keepX<- c(5:10, seq(15, 50, 5), seq(60, 100, 10))
ncomp = 4 
splsda.tune = tune.splsda(norm.TSS, Y, 
                          ncomp = 5, # (ncomp + 1)
                          test.keepX = list.keepX,
                          logratio = 'CLR',
                          validation = 'Mfold',
                          folds = 5, dist = "centroids.dist", nrepeat = 10) #recomended between 10-50
# may show some convergence issues for some of the cases, it is ok for tuning

####
plot(splsda.tune, optimal = TRUE, sd = TRUE)

### good one based on tune

select.keepX = splsda.tune$choice.keepX

res.splsda = splsda(X = norm.TSS,
                    Y = Y,
                    ncomp = 5,
                    keepX = select.keepX,
                    logratio= "CLR")

saveRDS(res.splsda, "additional_data/Morgan_et_al/ipaa_res.splda.rds")
###how did we do after tune

#on the first two components
plotIndiv(res.splsda,
          comp =c(1,2),
          col.per.group = color.mixo(1:2), # the number of groups 
          ellipse = TRUE,
          legend = TRUE,
          title = 'IPAA, sPLSDA comp 1 - 2')


plotter_main = data.frame(comp1=res.splsda$variates$X[,1], 
                     comp2=res.splsda$variates$X[,2],
                     Inflammation=meta_otu$group3)

g_main = ggplot(plotter_main, aes(comp1, comp2, color = Inflammation)) + 
  geom_point(size=2) + scale_color_manual(values=c("Inflamed"="tomato", "Uninflamed"="steelblue1")) +
  theme_bw() + theme_void() + theme(legend.position = 'none')


plotter_1 = selectVar(res.splsda, comp = 1)$value
plotter_1$otu.id = rownames(plotter_1)
plotter_1$group = "Inflamed"
plotter_1$group[plotter_1$value.var > 0] <- "Uninflamed"

plotter_1 = plotter_1[order(abs(plotter_1$value.var), decreasing=T),]
plotter_1$otu.id = factor(plotter_1$otu.id, levels = rev(as.character(plotter_1$otu.id)))

plotter1.1 = plotter_1[1:15,]
plotter1.1$names = make.unique(as.character(otu_names[as.character(plotter1.1$otu.id),]$Consensus.Lineage))
plotter1.1$names = factor(plotter1.1$names, levels = as.character(plotter1.1$names))

g1 = ggplot(plotter1.1, aes(value.var, names, fill = group)) + 
  geom_col() + scale_fill_manual(values=c("Inflamed"="tomato", "Uninflamed"="steelblue1")) +
  theme_bw() + theme(legend.position = 'none') +
  xlab("Importance") + ylab("Taxa")


plotter_2 = selectVar(res.splsda, comp = 2)$value
plotter_2$otu.id = rownames(plotter_2)
plotter_2$group = "Inflamed"
plotter_2$group[plotter_2$value.var > 0] <- "Uninflamed"

plotter_2 = plotter_2[order(abs(plotter_2$value.var), decreasing=T),]
plotter_2$otu.id = factor(plotter_2$otu.id, levels = as.character(plotter_2$otu.id))

plotter_2$names = make.unique(as.character(otu_names[as.character(plotter_2$otu.id),]$Consensus.Lineage))
plotter_2$names = factor(plotter_2$names, levels = as.character(plotter_2$names))

g2 = ggplot(plotter_2, aes(value.var, names, fill = group)) + 
  geom_col() + scale_fill_manual(values=c("Inflamed"="tomato", "Uninflamed"="steelblue1")) +
  theme_bw() + coord_flip() + 
  theme(legend.position = 'none',
        axis.text.x=element_text(angle=90, hjust=1)) +
  xlab("Importance") + ylab("Taxa")

g_blank = ggplot() + theme_void()

pdf("FIGURES/FIGURE5/Figure5d1_ipaa_16s_comp.pdf", height =5, width = 5)
grid.arrange(g_main)
dev.off()

pdf("FIGURES/FIGURE5/Figure5d2_ipaa_16s_comp.pdf", height =5, width = 10)
grid.arrange(g1)
dev.off()

pdf("FIGURES/FIGURE5/Figure5d3_ipaa_16s_comp.pdf", height =10, width = 5)
grid.arrange(g2)
dev.off()




#
#
#
#
### or just full on sPLS
X = t(counts_otus)
Y = norm.TSS

ncomp = 10
result.spls <- spls(X, Y, ncomp = ncomp, keepX = c(rep(10, ncomp)), mode = 'canonical')

cord.X = cor(result.spls$X, result.spls$variates$X[, 1:ncomp], use = "pairwise")
cord.Y = cor(result.spls$Y, result.spls$variates$Y[, 1:ncomp], use = "pairwise")
true_sim = cord.X %*% t(cord.Y)


selector = subset(ct_ref, cluster == "mono_mac2")
mono_mac2_sim = true_sim[rownames(selector),]
mono_mac2_sim[is.na(mono_mac2_sim)] <- 0
keeper_c = which(colVars(mono_mac2_sim)>0.02)

mono_mac2_sim = mono_mac2_sim[,keeper_c]
sim_names = otu_names[colnames(mono_mac2_sim),]
#colnames(mono_mac2_sim) <- gsub(".*;", "", sim_names$Consensus.Lineage)
colnames(mono_mac2_sim) <- sim_names$Consensus.Lineage
dim(mono_mac2_sim)


temp=t(mono_mac2_sim)
temp_agg=aggregate(temp, by=list(sim_names$Consensus.Lineage), 'mean')
rownames(temp_agg)=temp_agg[,1]
temp_agg=temp_agg[,-1]
mono_mac_sim3=t(temp_agg)

mono_mac_sim3=mono_mac_sim3[which(rowVars(mono_mac_sim3)>0),]

heatty = grid.grabExpr(draw(
  Heatmap(mono_mac_sim3, show_row_names = T, show_column_names = T,
          heatmap_legend_param = list(title = "Covariate\nValue"),
          cluster_rows = T, cluster_columns = T, row_names_side = 'left',
          column_names_gp = gpar(fontsize=10),
          row_names_gp = gpar(fontsize=9),
          column_names_rot = 60, 
          show_row_dend=F,
          #row_title_gp = gpar(fontsize = 10),
          #row_names_max_width = unit(10,'cm'),
          use_raster = F,
          col = colorRamp2(c(-0.5,0,0.5), c("blue","white" ,"red"))),
  padding = unit(c(2, 2, 3, 2), "cm")
))

namer <- "FIGURES/FIGURE5/Figure5d_otu_array_spls.pdf"
pdf(namer, height = 9, width = 9)
grid.arrange(heatty)
#grid.arrange(ipaa_p, heatty, nrow=2, heights = c(1.3,1))
dev.off()


##### null distr


true_flat <- c(cord.X %*% t(cord.Y))
true_flat[is.na(true_flat)] <- 0

num_sim <- 999
null_mat <- matrix(NA, 1, ncol(X)*ncol(Y))
for(j in 1:num_sim){
  null.spls <- spls(X, Y[sample(rownames(Y)),], ncomp = ncomp)
  cord.X = cor(null.spls$X, null.spls$variates$X[, 1:ncomp], use = "pairwise")
  cord.Y = cor(null.spls$Y, null.spls$variates$Y[, 1:ncomp], use = "pairwise")
  null_sim = cord.X %*% t(cord.Y)
  
  null_flat <- c(null_sim)
  null_mat <- rbind(null_mat, null_flat)
}
null_mat <- null_mat[-1,]
null_mat[is.na(null_mat)] <- 0

p.value = NULL
for (k in 1:ncol(null_mat)){
  p=mean(abs(null_mat[,k])>=abs(true_flat[k]))
  p.value = c(p.value,p)
}
p.adj = p.adjust(p.value, "BH")
p.adj = p.value
p_mat <- matrix(p.adj, ncol(X), ncol(Y))
rownames(p_mat) <- colnames(X)
colnames(p_mat) <- colnames(Y)


#
#
#
true_sim2=true_sim

for(jj in 1:ncol(true_sim2)){
  heyer=true_sim2[,jj]
  heyer[p_mat[,jj]>0.025]<-0
  true_sim2[,jj]=heyer
}

selector = subset(ct_ref, cluster == "mono_mac2")
mono_mac2_sim = true_sim2[rownames(selector),]
mono_mac2_sim[is.na(mono_mac2_sim)] <- 0
keeper_c = which(colSums(mono_mac2_sim)>0)

mono_mac2_sim = mono_mac2_sim[rownames(selector),keeper_c]
mono_mac2_sim[is.na(mono_mac2_sim)] <- 0
sim_names = otu_names[colnames(mono_mac2_sim),]

temp=t(mono_mac2_sim)
temp_agg=aggregate(temp, by=list(sim_names$Consensus.Lineage), 'mean')
rownames(temp_agg)=temp_agg[,1]
temp_agg=temp_agg[,-1]
mono_mac_sim3=t(temp_agg)

mono_mac_sim3=mono_mac_sim3[which(rowSums(mono_mac_sim3)>0),]

#colnames(mono_mac_sim3) <- gsub("eae;g__", "eae", colnames(mono_mac_sim3))
#colnames(mono_mac_sim3) <- gsub(".*;", "", colnames(mono_mac_sim3))
dim(mono_mac_sim3)

heatty = grid.grabExpr(draw(
  Heatmap(mono_mac_sim3, show_row_names = T, show_column_names = T,
          heatmap_legend_param = list(title = "Covariate\nValue"),
          cluster_rows = T, cluster_columns = T, row_names_side = 'left',
          column_names_gp = gpar(fontsize=10),
          show_column_dend = F, show_row_dend = F,
          row_names_gp = gpar(fontsize=9),
          column_names_rot = 60,
          #row_title_gp = gpar(fontsize = 10),
          #row_names_max_width = unit(10,'cm'),
          use_raster = F,
          col = colorRamp2(c(-0.5,0,0.5), c("blue","white" ,"red"))),
  padding = unit(c(2, 2, 3, 2), "cm")
))

namer <- "FIGURES/FIGURE5/Figure5d_otu_array_spls.pdf"
pdf(namer, height = 8, width = 9)
grid.arrange(heatty)
#grid.arrange(ipaa_p, heatty, nrow=2, heights = c(1.3,1))
dev.off()
