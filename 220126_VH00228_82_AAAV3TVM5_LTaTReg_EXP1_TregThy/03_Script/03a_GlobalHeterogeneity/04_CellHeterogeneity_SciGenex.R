# #########################################
# This script reads and filters sc10x  data
# #########################################




# Study heterogeneity with SciGeneX
# .................................

## @knitr cell_heterogeneity_scigenex

sc10x_rna_df <- as.data.frame(sc10x@assays$RNA@data)

# Run SciGenex
res_scigenex <- DBFMCL( sc10x_rna_df, k = 80, inflation = 8)

# Plot simulated distances and observed distances for the k nearest neighbors of each gene
viz_dist( res_scigenex)

# Select top 20 genes and use them for the heatmap visualization
res_scigenex <- top_genes(res_scigenex, top = 20)
plot_heatmap(res_scigenex, use_top_genes = TRUE)

sc10x_rna_df[sc10x_rna_df > 0] <- 1

for( cluster_name in unique( res_scigenex@cluster)){
  genes = names(res_scigenex@cluster)[ which( res_scigenex@cluster == cluster_name)]
  sc10x = AddModuleScore( sc10x, list( genes), name = paste0( "scigenex.", cluster_name))
  FeaturePlot( sc10x, features = paste0( "scigenex.", cluster_name, "1"))
}
