library(Seurat)
library(ggplot2)
library(future)

plan(strategy = "multicore", workers = 10)
options(future.globals.maxSize = 50 * 1024 ^ 3)

##### Smillie

library(Matrix)
# counts_matrix = readMM("/gpfs/data/ruggleslab/home/devlij03/IBD/MANUSCRIPT_FILES/additional_data/Smillie_et_al/matrix.mtx")
# 
# genes = as.character(read.table("/gpfs/data/ruggleslab/home/devlij03/IBD/MANUSCRIPT_FILES/additional_data/Smillie_et_al/features.tsv", F, "\t")[,1])
# barcodes = as.character(read.table("/gpfs/data/ruggleslab/home/devlij03/IBD/MANUSCRIPT_FILES/additional_data/Smillie_et_al/barcodes.tsv", F, '\t')[,1])
# 
# colnames(counts_matrix) = barcodes
# rownames(counts_matrix) = genes
# 
#
# print("data read")
# 
# smillie_obj = CreateSeuratObject(counts = counts_matrix, project = "Smillie",
#                            names.field = 1, names.delim = ":")
#
#
# meta = read.table("additional_data/Smillie_et_al/ramnik_meta.txt", T,'\t')
# meta = meta[-1,]
# rownames(meta) <- as.character(meta$NAME)
# meta = meta[colnames(smillie_obj),]
#
# smillie_obj@meta.data = cbind(smillie_obj@meta.data, meta)
#
# 
# print("s_obj created")
# 
# smillie_obj = SCTransform(smillie_obj, ncells = 50000)
# 
# print("SCT complete")
# 
# saveRDS(smillie_obj, "OBJECTS/Smillie/seurat_obj.rds")

smillie_obj = readRDS("OBJECTS/Smillie/seurat_obj.rds")

##### Martins

load("additional_data/Martin_et_al/model_and_samples_ileum.rd")
counts_matrix = ileum_ldm$dataset$umitab

meta=read.table("additional_data/Martin_et_al/GIMAT_full_meta.txt", T, sep='\t')

martin_obj = CreateSeuratObject(counts = counts_matrix, project = "Martins",
                           names.field = 1, names.delim = ":")

print("s_obj created")

martin_obj@meta.data = cbind(martin_obj@meta.data, meta)
  
martin_obj = SCTransform(martin_obj, ncells = 30000)

saveRDS(martin_obj, "OBJECTS/Martin/seurat_obj.rds")

#####

##### pouch

pouch_obj = readRDS("OBJECTS/seurat_obj.rds")
pouch_obj = SCTransform(pouch_obj, ncells=25000)
pouch_obj@meta.data$orig.ident2 = pouch_obj@meta.data$orig.ident
pouch_obj@meta.data$orig.ident = "Pouch"

####

####

####

hca.list=list(pouch_obj, smillie_obj, martin_obj)

features <- SelectIntegrationFeatures(object.list = hca.list, nfeatures = 3000)
hca.list <- PrepSCTIntegration(object.list = hca.list, anchor.features = features)
hca.list <- lapply(X = hca.list, FUN = RunPCA, verbose = FALSE, features = features)
anchors <- FindIntegrationAnchors(object.list = hca.list, anchor.features = features, normalization.method = "SCT", reduction = "rpca")
hca.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

hca.integrated <- RunPCA(hca.integrated, verbose = FALSE)
hca.integrated <- RunUMAP(hca.integrated, dims = 1:30)
feat_map = DimPlot(hca.integrated, group.by = "orig.ident")

ggsave("OBJECTS/Atlas/atlas_umap.png", feat_map, height=10, width = 12, units='in')

saveRDS(hca.integrated, "OBJECTS/Atlas/seurat_obj.rds")
