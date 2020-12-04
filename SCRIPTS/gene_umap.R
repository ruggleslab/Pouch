library(Seurat)
library(ggplot2)

gradient_colors = c("gray85", "purple4")

gene_umap = function(s_obj, markers, pts = 2, cols=NA){
  
  if(is.na(cols)){
    
    plot_umap =
      FeaturePlot(s_obj, features=markers, 
                  reduction = "umap", 
                  cells = sample(colnames(s_obj)), pt.size=pts, 
                  cols = gradient_colors, ncol = length(markers), 
                  min.cutoff=0)
  } else {
    plot_umap =
      FeaturePlot(s_obj, features=markers, 
                  reduction = "umap", 
                  cells = sample(colnames(s_obj)), pt.size=pts,
                  cols = gradient_colors, ncol = cols, 
                  min.cutoff=0)
  }
  return(plot_umap)
}

