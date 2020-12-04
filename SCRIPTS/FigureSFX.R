library(Seurat)
library(ggplot2)
library(gridExtra)

colorer <- c("Normal_pouch"="dodgerblue2","Pouchitis"="red3", "UC_colon"="forestgreen", "UC_inflamed"="darkorange1")

s_obj = readRDS("OBJECTS/seurat_obj.rds")


df = data.frame(table(s_obj@meta.data$orig.ident, s_obj@meta.data$pouch.status))
df =subset(df, Freq>0)
df=df[order(df$Freq, decreasing=T),]

df$Var1=factor(df$Var1, levels=unique(df$Var1))
colnames(df)[2] = "DiseaseCondition"

gg1=ggplot(df, aes(Var1, Freq, fill=DiseaseCondition)) + 
  geom_col() + theme_bw() +
  ylab("Number of cells") + xlab("Sample") +
  scale_fill_manual(values=colorer) +
  theme(axis.text.x=element_blank())

pairwise.wilcox.test(df$Freq, df$DiseaseCondition)

gg2=ggplot(df, aes(DiseaseCondition, Freq, color=DiseaseCondition)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(width=0.2, size=0.5, col='black')+
  theme_bw() +
  ylab("Number of cells") + xlab("") +
  scale_color_manual(values=colorer) +
  theme(axis.text.x=element_blank())

gg=arrangeGrob(gg1,gg2, nrow=1)

ggsave("FIGURES/FIGURE1/cell_hist.pdf", gg, height=2, width = 8, units='in')