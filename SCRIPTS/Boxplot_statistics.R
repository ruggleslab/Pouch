library(stats)

stats_files = c("OBJECTS/clusters-resolutions/clust_box_MajorPopulations_values.txt",
                "OBJECTS/T_cells/clusters-resolutions/clust_box_MinorPopulations_values.txt",
                "OBJECTS/B_cells/clusters-resolutions/clust_box_MinorPopulations_values.txt",
                "OBJECTS/Myeloid_cells/clusters-resolutions/clust_box_MinorPopulations2_values.txt")

paths = c("OBJECTS/clusters-resolutions/clust_box_MajorPopulations_stats.txt",
          "OBJECTS/T_cells/clusters-resolutions/clust_box_MinorPopulations_stats.txt",
          "OBJECTS/B_cells/clusters-resolutions/clust_box_MinorPopulations_stats.txt",
          "OBJECTS/Myeloid_cells/clusters-resolutions/clust_box_MinorPopulations2_stats.txt")

for(i in 1:length(stats_files)){
  dfm=read.table(stats_files[i], header=T, sep='\t')

  cells = as.character(unique(dfm$Cluster))
  groups="PouchStatus"
  
  res_df = data.frame(cluster=NA, comps=NA, pval=NA)
  for(j in 1:length(cells)){
    curr=subset(dfm, Cluster==cells[j])
    x=curr[,"perc"]
    g=curr[,groups]
    ressy=pairwise.wilcox.test(x,g)$p.value
    pvals=as.numeric(ressy)
    comps_add=paste0(rep(rownames(ressy), ncol(ressy)), ";", rep(colnames(ressy), each=nrow(ressy)))
    cells_add=rep(cells[j], length(pvals))
    adder=data.frame(cluster=cells_add, comps=comps_add, pval=pvals)
    adder=subset(adder, !is.na(pval))
    res_df=rbind(res_df, adder)
  }
  res_df=res_df[-1,]
  write.table(res_df, paths[i], quote=F, row.names=F, sep='\t')
}