library(Seurat)
library(ggplot2)
library(future)

plan(strategy = "multicore", workers = 10)
options(future.globals.maxSize = 50 * 1024 ^ 3)



##### Smillie epi

library(Matrix)
counts_matrix = readMM("/gpfs/data/ruggleslab/home/devlij03/IBD/MANUSCRIPT_FILES/additional_data/Smillie_et_al/gene_sorted-Fib.matrix.mtx")

genes = as.character(read.table("/gpfs/data/ruggleslab/home/devlij03/IBD/MANUSCRIPT_FILES/additional_data/Smillie_et_al/Fib.genes.tsv", F, "\t")[,1])
barcodes = as.character(read.table("/gpfs/data/ruggleslab/home/devlij03/IBD/MANUSCRIPT_FILES/additional_data/Smillie_et_al/Fib.barcodes2.tsv", F, '\t')[,1])

colnames(counts_matrix) = barcodes
rownames(counts_matrix) = genes


print("data read")

smillie_obj = CreateSeuratObject(counts = counts_matrix, project = "Smillie",
                                 names.field = 1, names.delim = ":")


meta = read.table("additional_data/Smillie_et_al/ramnik_meta.txt", T,'\t')
meta = meta[-1,]
rownames(meta) <- as.character(meta$NAME)
meta = meta[colnames(smillie_obj),]

smillie_obj@meta.data = cbind(smillie_obj@meta.data, meta)


print("s_obj created")

smillie_obj = SCTransform(smillie_obj, ncells = 50000)

print("SCT complete")

saveRDS(smillie_obj, "OBJECTS/Smillie/seurat_fib_obj.rds")







##### Smillie epi

library(Matrix)
counts_matrix = readMM("/gpfs/data/ruggleslab/home/devlij03/IBD/MANUSCRIPT_FILES/additional_data/Smillie_et_al/gene_sorted-Epi.matrix.mtx")

genes = as.character(read.table("/gpfs/data/ruggleslab/home/devlij03/IBD/MANUSCRIPT_FILES/additional_data/Smillie_et_al/Epi.genes.tsv", F, "\t")[,1])
barcodes = as.character(read.table("/gpfs/data/ruggleslab/home/devlij03/IBD/MANUSCRIPT_FILES/additional_data/Smillie_et_al/Epi.barcodes2.tsv", F, '\t')[,1])

colnames(counts_matrix) = barcodes
rownames(counts_matrix) = genes


print("data read")

smillie_obj = CreateSeuratObject(counts = counts_matrix, project = "Smillie",
                           names.field = 1, names.delim = ":")


meta = read.table("additional_data/Smillie_et_al/ramnik_meta.txt", T,'\t')
meta = meta[-1,]
rownames(meta) <- as.character(meta$NAME)
meta = meta[colnames(smillie_obj),]

smillie_obj@meta.data = cbind(smillie_obj@meta.data, meta)


print("s_obj created")

smillie_obj = SCTransform(smillie_obj, ncells = 50000)

print("SCT complete")

saveRDS(smillie_obj, "OBJECTS/Smillie/seurat_epi_obj.rds")

