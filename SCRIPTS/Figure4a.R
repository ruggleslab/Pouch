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
library(broom)

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

s_obj = readRDS("OBJECTS/seurat_obj.rds")
meta <- s_obj@meta.data

#
#
#
#
#
#
### mono_mac1
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
soxresser = resser

#

#
#
#
#

##
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
il1bresser = resser


#
#
#
#
#

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
apoeresser = resser
#
#
#
#
# foxp3
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



#
#
#
#### any relationship between foxp3 and mono_macs?

fox_sox <- data.frame(foxresser[,1:2], foxp3=foxresser[,3], sox4=soxresser[,3])
pval = round(glance(lm(foxp3~sox4 + pouchStatus, fox_sox))$p.value,4)

r2 = round(summary(lm(foxp3~sox4 + pouchStatus, fox_sox))$adj.r.squared,3)

fox_sox_p = 
  ggplot(fox_sox, aes(foxp3*100, sox4*100, color=pouchStatus)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, color="black", linetype="dashed") +
  ggtitle("Foxp3 versus mono_mac1") +
  annotate("text", x = 8, y = 90*max(fox_sox$sox4), label = paste0("p = ", pval, ", R^2 = ", r2), color="red3") +
  xlab("FOXP3+/BATF+ Cell percentage") + ylab("SOX4+/MAFA+ Cell Percentage") +
  #geom_label_repel() + 
  theme_bw() +
  scale_color_manual(values=colorer)

##

fox_il1b <- data.frame(foxresser[,1:2], foxp3=foxresser[,3], il1b=il1bresser[,3])
pval = round(glance(lm(foxp3~il1b + pouchStatus, fox_il1b))$p.value,4)

r2 = round(summary(lm(foxp3~il1b + pouchStatus, fox_il1b))$adj.r.squared,3)

fox_il1b_p = 
  ggplot(fox_il1b, aes(foxp3*100, il1b*100, color=pouchStatus)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, color="black", linetype="dashed") +
  ggtitle("Foxp3 versus mono_mac2") +
  annotate("text", x = 8, y = 90*max(fox_il1b$il1b), label = paste0("p = ", pval, ", R^2 = ", r2), color="red3") +
  xlab("FOXP3+/BATF+ Cell percentage") + ylab("IL1B+/LYZ+ Cell Percentage") +
  #geom_label_repel() + 
  theme_bw() +
  scale_color_manual(values=colorer)

##
fox_apoe <- data.frame(foxresser[,1:2], foxp3=foxresser[,3], apoe=apoeresser[,3])
pval = round(glance(lm(foxp3~apoe + pouchStatus, fox_apoe))$p.value,4)

r2 = round(summary(lm(foxp3~apoe + pouchStatus, fox_apoe))$adj.r.squared,3)

fox_apoe_p = 
  ggplot(fox_apoe, aes(foxp3*100, apoe*100, color=pouchStatus)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE, color="black", linetype="dashed") +
  ggtitle("Foxp3 versus mono_mac3") +
  annotate("text", x = 8, y = 90*max(fox_apoe$apoe), label = paste0("p = ", pval, ", R^2 = ", r2), color="red3") +
  xlab("FOXP3+/BATF+ Cell percentage") + ylab("APOE+/C1QC+ Cell Percentage") +
  #geom_label_repel() + 
  theme_bw() +
  scale_color_manual(values=colorer)


pdf("FIGURES/FIGURE4/Figure4a_lm.pdf", height = 3.5, width = 14)
grid.arrange(fox_sox_p,
             fox_il1b_p,
             fox_apoe_p, nrow=1)
dev.off()