library(ggplot2)

sig_means = read.table("OBJECTS/Interactions/significant_means.txt", T, '\t')

library(ComplexHeatmap)
library(circlize)
library(grid)
library(scales)
library(reshape)

meta = sig_means[,1:12]
values = sig_means[13:ncol(sig_means)]

cells = c("act_CD8_NR4A2_hi", 
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
          "pdcs", "DC",
          "mono_mac1", 
          "mono_mac2",
          "mono_mac3")

#####
library(igraph)
library(reshape)

geneInts = make.unique(paste0(meta$gene_a, "__", meta$gene_b))

name_df <- melt(values)
name_df$interaction = rep(geneInts, 529)
name_df$ct1 = gsub("\\..*", "", name_df$variable)
name_df$ct2 = gsub(".*\\.", "", name_df$variable)
name_df$value[is.na(name_df$value)] <- 0
cell_interactions = subset(name_df, value > 0)
cell_interactions = data.frame(cell_interactions[,c(4,5,3,2,1)])

keeper = c()
for(i in 1:nrow(cell_interactions)){
  ct1 = cell_interactions$ct1[i]
  ct2 = cell_interactions$ct2[i]
  if(ct1 != ct2){
    keeper = c(keeper, rownames(cell_interactions)[i])
  }
  
}
cell_interactions = cell_interactions[keeper,]
cell_interactions = subset(cell_interactions, value > 2)

cell_types = data.frame(cellTypes = unique(name_df$ct1))

g <- graph_from_data_frame(cell_interactions, directed=FALSE, vertices=cell_types)
E(g)$weight <- cell_interactions$value

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

E(g)$color[E(g)$weight > 0] <- add.alpha('red', alpha=0.25)

#edge thickness

edge_df <- data.frame(as_edgelist(g), weight = E(g)$weight)
edge_df$thickness <- rescale(abs(edge_df$weight), to=c(0,8))

# set labels
cell_types$dataType = "T_cells"
cell_types$dataType[grepl("follic|b_ce|gc|plasma", cell_types$cellTypes)] <- "B_cells"
cell_types$dataType[grepl("mono|mast|pdc|DC", cell_types$cellTypes)] <- "Myeloid_cells"

coul <- c(Myeloid_cells="dodgerblue2", T_cells="darkseagreen1", B_cells="coral")
my_color=coul[cell_types$dataType]

# vertex adjust
rownames(cell_types) = cell_types$cellTypes
cells = as.character(cell_types$cellTypes)
totaler = c()
for(j in 1:length(cells)){
  ct1_tot = nrow(subset(cell_interactions, ct1 == cells[j]))
  ct2_tot = nrow(subset(cell_interactions, ct2 == cells[j]))
  tot = as.numeric(ct1_tot)+as.numeric(ct2_tot)
  totaler = c(totaler, tot)
}
cell_types$edge_num = totaler
cell_types$edge_adj = (round(cell_types$edge_num/20,0)+1)*2.5
#



bound_x_hi = c(60,80,100)
bound_x_lo = c(0,-10,75)
bound_y_hi = c(60,100,60)
bound_y_lo = c(0,50,20)

new_axes=0
new_axes=data.frame(x=rep(0, nrow(cell_types)), y=rep(0, nrow(cell_types)))
for(i in 1:nrow(cell_types)){
  ct = cell_types$dataType[i]
  if(ct == "T_cells"){
    new_axes[,1][i] <- sample(bound_x_lo[1]:bound_x_hi[1],1)
    new_axes[,2][i] <- sample(bound_y_lo[1]:bound_y_hi[1],1)
  } else if(ct == "B_cells"){
    new_axes[,1][i] <- sample(bound_x_lo[3]:bound_x_hi[3],1)
    new_axes[,2][i] <- sample(bound_y_lo[3]:bound_y_hi[3],1)
  } else if (ct=="Myeloid_cells"){
    new_axes[,1][i] <- sample(bound_x_lo[2]:bound_x_hi[2],1)
    new_axes[,2][i] <- sample(bound_y_lo[2]:bound_y_hi[2],1)
  }
}
new_axes2=as.matrix(new_axes)
#cell_types$x = new_axes[,1]
#cell_types$y = new_axes[,2]

pdf("OBJECTS/Interactions/interaction_network2.pdf", height = 25, width = 30)
par(bg="white", mar=c(0,0,0,0))
set.seed(4)
plot(g, #layout=as.matrix(better_coords[,5:6]),
     #layout=as.matrix(p_plot[,5:6]),
     #vertex.size=5,
     layout=new_axes2,
     edge.width = as.matrix(edge_df$thickness),
     vertex.size=as.matrix(cell_types$edge_adj),
     vertex.shape='circle',
     vertex.color=my_color, 
     vertex.label.cex=1.2,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     vertex.frame.color="transparent"
)
dev.off()



#
#
#
#

#
#
#

#### and a dotplot for figure 6
library(reshape2)

all_means = read.table("OBJECTS/Interactions/means.txt", T, '\t')
all_ps = read.table("OBJECTS/Interactions/pvalues.txt", T, '\t')

meta = all_means[,1:11]
values = all_means[12:ncol(all_means)]
pvalues = all_ps[12:ncol(all_ps)]

geneInts = make.unique(paste0(meta$gene_a, "__", meta$gene_b))

name_df <- melt(values)
p_df = melt(pvalues)
name_df$pvalue = p_df$value
name_df$interaction = rep(geneInts, ncol(values))
name_df$ct1 = gsub("\\..*", "", name_df$variable)
name_df$ct2 = gsub(".*\\.", "", name_df$variable)
name_df$value[is.na(name_df$value)] <- 0
cell_interactions = subset(name_df, value > 0)

cell_interactions2 = subset(cell_interactions, value > 0 & pvalue < 0.05)
cell_interactions2 = cell_interactions2[order(cell_interactions2$ct1, cell_interactions2$ct2),]
cell_interactions2$variable = factor(cell_interactions2$variable, levels = unique(cell_interactions2$variable))

cell_interactions2$int = log2(cell_interactions2$value+1)
cell_interactions2 = subset(cell_interactions2, grepl("mono|DC|Th17|Foxp3", cell_interactions2$ct1))
cell_interactions2 = subset(cell_interactions2, grepl("mono|DC|Th17|Foxp3", cell_interactions2$ct2))

int_pairs = c("ALOX5__ALOX5AP",
              "CCR6__CCL20",
              "CCL5__CCR5",
              "CCL5__CCR4",
              "CCL5__CCR1",
              "CCL4__CCR5",
              "CCL3__CCR5",
              "CCL3__CCR1",
              "CD40__CD40LG",
              "CD74__APP",
              "CSF1R__CSF1",
              "CXCR6__CXCL16",
              "CTLA4__CD86",
              "CTLA4__CD80",
              "CSF2__CSF1R",
              "ICOSLG__ICOS",
              "ICAM1__AREG",
              "LTA__TNFRSF1A",
              "PDCD1__FAM3C",
              "PDCD1__CD274",
              "NRP2__VEGFA")

int_pairs = unique(c("LGALS9__HAVCR2",
  "LGALS9__CD44",
  "AXL__GAS6",
  "TNFRSF10B__FASLG",
  "TNF__TNFRSF1B",
  "TNF__TNFRSF1A",
  "TNFRSF10B__TNFSF10",
  "CCR6__CCL20",
  "CCL5__CCR5",
  "CCL5__CCR4",
  "CCL5__CCR1",
  "CCL4__CCR5",
  "CCL3__CCR5",
  "CCL3__CCR1",
  "CD40__CD40LG",
  "CD74__APP",
  "CSF1R__CSF1",
  "CXCR6__CXCL16",
  "CTLA4__CD86",
  "CTLA4__CD80",
  "CSF2__CSF1R",
  "ICOSLG__ICOS",
  "LTA__TNFRSF1A",
  "PDCD1__CD274",
  "NRP2__VEGFA", #all next is DCs
  "HLA-DPB1__TNFSF13B",
  "HLA-DPA1__TNFSF9",
  "IL1B__ADRB2", "CD74__APP",
  "CCL3__CCR1"))

cell_interactions2 = subset(cell_interactions2, interaction %in% int_pairs)


#### get order
cell_interactions3 = subset(cell_interactions, interaction %in% int_pairs)
cell_interactions3 = subset(cell_interactions3, grepl("mono|DC|Th17|Foxp3", cell_interactions3$ct1))
cell_interactions3 = subset(cell_interactions3, grepl("mono|DC|Th17|Foxp3", cell_interactions3$ct2))
#cell_interactions3=subset(cell_interactions3, ct1!=ct2)

new_mat = cast(cell_interactions3, ct1+ct2~interaction, value="value")
rownames(new_mat)=paste0(new_mat[,1], ":", new_mat[,2])
new_mat=new_mat[,3:ncol(new_mat)]
new_mat2=t(new_mat)
colnames(new_mat2)=rownames(new_mat)
rownames(new_mat2)=colnames(new_mat)

new_order = hclust(dist(new_mat2))$labels[hclust(dist(new_mat2))$order]

###

m4 = length(unique(subset(cell_interactions2, ct1 == "DC")$variable))+0.5
ffs = length(unique(subset(cell_interactions2, ct1 == "Foxp3_cells")$variable))+m4
m1 = length(unique(subset(cell_interactions2, ct1 == "mono_mac1")$variable))+ffs
m2 = length(unique(subset(cell_interactions2, ct1 == "mono_mac2")$variable))+m1
m3 = length(unique(subset(cell_interactions2, ct1 == "mono_mac3")$variable))+m2
ths = length(unique(subset(cell_interactions2, ct1 == "Th17")$variable))+m3

cell_interactions2$interaction = factor(cell_interactions2$interaction, levels = new_order)
gg_dots = ggplot(cell_interactions2, aes(variable, interaction, color = int)) + 
  geom_point() + theme_classic() +
  scale_color_gradient(low="navajowhite", high = "red3") + 
  geom_vline(xintercept = ffs, color = "black", size=0.6, linetype="dashed") +
  geom_vline(xintercept = m1, color = "black", size=0.6, linetype="dashed") +
  geom_vline(xintercept = m2, color = "black", size=0.6, linetype="dashed") +
  geom_vline(xintercept = m3, color = "black", size=0.6, linetype='dashed') +
  geom_vline(xintercept = m4, color = "black", size=0.6, linetype="dashed") +
  geom_vline(xintercept = ths, color = "black", size=0.6, linetype="dashed") +
  theme(axis.text.x = element_text(angle = 60, color="black", hjust=1),
        axis.text.y = element_text(color="black"))

pdf("FIGURES/FIGURE6/Figure6b_dotplot_rev.pdf", height = 5, width = 7)
gg_dots
dev.off()




