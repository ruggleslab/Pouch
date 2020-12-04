library(Seurat)
library(ggplot2)

s_obj <- readRDS("OBJECTS/seurat_obj.rds")

tcells <- subset(s_obj, subset = snn_res.0.4 %in% c("C02", "C03", "C05", "C06", "C07", "C08"))
bcells <- subset(s_obj, subset = snn_res.0.4 %in% c("C01", "C10", "C11", "C04", "C13"))
mcells <- subset(s_obj, subset = snn_res.0.4 %in% c("C09", "C12", "C14"))

#dir.create("OBJECTS/T_cells")
#dir.create("OBJECTS/B_cells")
#dir.create("OBJECTS/Myeloid_cells")

saveRDS(tcells, "OBJECTS/T_cells/seurat_obj.rds")
saveRDS(bcells, "OBJECTS/B_cells/seurat_obj.rds")
saveRDS(mcells, "OBJECTS/Myeloid_cells/seurat_obj.rds")
