library(Seurat)

s_obj <- readRDS("OBJECTS/seurat_obj.rds")

counts = s_obj@assays$RNA@counts[rownames(s_obj@assays$integrated@data),]
counts = data.frame(Gene=rownames(counts), counts)
write.table(counts, "OBJECTS/Interactions/all_counts.txt", row.names=F, sep='\t', quote=F)

meta = data.frame(Cell=colnames(s_obj), cell_type=s_obj@meta.data$MinorPopulations)
meta$Cell = make.names(meta$Cell)
write.table(meta, "OBJECTS/Interactions/all_meta.txt", sep='\t', row.names=F, quote=F)