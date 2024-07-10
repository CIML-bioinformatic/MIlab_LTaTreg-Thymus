# ########################################################################################
# This script aims to perform specific analysis of chosen KEGG Pathways
# ########################################################################################

## @knitr kegg_pathway_analysis

# Get the KEGG reference
KEGG_pathway_gene_df = getGeneKEGGLinks( species="mmu")
KEGG_pathway_gene_df$Symbol = mapIds( org.Mm.eg.db, KEGG_pathway_gene_df$GeneID, column="SYMBOL", keytype="ENTREZID")

# For each KEGG pathway to analyse perform a AddModuleScore analysis with plits and statistics
for( kegg_pathway in names( MONITORED_KEGG_PATHWAY_LIST)){
  
  cat(" \n \n"); # Required for '.tabset'
  cat("### Pathway", kegg_pathway)
  cat(" \n \n"); # Required for '.tabset'
  
  # Get the KEGG ID of the pathway
  kegg_pathway_id = MONITORED_KEGG_PATHWAY_LIST[[ kegg_pathway]]
  
  # Get the genes symbols for the gene EntrezID
  gene_symbols = na.exclude( KEGG_pathway_gene_df[ which( KEGG_pathway_gene_df$PathwayID == paste0( "path:", kegg_pathway_id)), "Symbol"])
  gene_symbols = intersect( gene_symbols, rownames( sc10x.seurat))

  # Build the matrix of means of expression by clusters and conditions
  heatmap_matrix_WT = data.frame()
  heatmap_matrix_KO = data.frame()
  Idents( sc10x.seurat) = "seurat_clusters"
  for( current_cluster in levels( Idents( sc10x.seurat))){
    WT_cells = which( sc10x.seurat@meta.data$HTO_classification == SAMPLE_WT & sc10x.seurat@meta.data$seurat_clusters == current_cluster)
    pathway_data_matrix_WT = data.frame( GetAssayData(object = sc10x.seurat, slot = "data")[ gene_symbols, WT_cells])
    heatmap_matrix_WT = rbind( heatmap_matrix_WT, rowMeans( pathway_data_matrix_WT))
    
    KO_cells = which( sc10x.seurat@meta.data$HTO_classification == SAMPLE_KO & sc10x.seurat@meta.data$seurat_clusters == current_cluster)
    pathway_data_matrix_KO = data.frame( GetAssayData(object = sc10x.seurat, slot = "data")[ gene_symbols, KO_cells])
    heatmap_matrix_KO = rbind( heatmap_matrix_KO, rowMeans( pathway_data_matrix_KO))
  }
  colnames( heatmap_matrix_WT) = gene_symbols
  rownames( heatmap_matrix_WT) = levels( Idents( sc10x.seurat))
  colnames( heatmap_matrix_KO) = gene_symbols
  rownames( heatmap_matrix_KO) = levels( Idents( sc10x.seurat))
  
  # Print the heatmap of the WT matrix
  print(
  Heatmap( as.matrix(  heatmap_matrix_WT), 
           cluster_rows = FALSE, 
           cluster_columns = TRUE,
           height = unit( length( levels( Idents( sc10x.seurat))) * 1, "cm"),
           column_title = paste0( "WT sample : Mean expression by cluster of\n", kegg_pathway, " pathway genes"),
           column_names_gp = gpar(fontsize = 6))
  )
  
  # Print the heatmap of the KO matrix
  print(
  Heatmap( as.matrix( heatmap_matrix_KO), 
           cluster_rows = FALSE, 
           cluster_columns = TRUE,
           height = unit( length( levels( Idents( sc10x.seurat))) * 1, "cm"),
           column_title = paste0( "KO sample : Mean expression by cluster of\n", kegg_pathway, " pathway genes"),
           column_names_gp = gpar(fontsize = 6))
  )
  
  # Print the heatmap of the KO - WT matrix and export it as PNG and CSV files
  heatmap_ko_wt = Heatmap( as.matrix( heatmap_matrix_KO - heatmap_matrix_WT), 
           cluster_rows = TRUE, 
           cluster_columns = TRUE,
           height = unit( length( levels( Idents( sc10x.seurat))) * 1, "cm"), 
           column_title = paste0( "Difference of mean expression (KO - WT) by cluster of\n", kegg_pathway, " pathway genes"),
           column_names_gp = gpar(fontsize = 12))
  
  print(
    heatmap_ko_wt
  )
  
  write.table( x= heatmap_matrix_KO - heatmap_matrix_WT, 
               file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "KEGG_", kegg_pathway, "_KO-WT_Heatmap.csv")),
               sep =",", col.names = NA, row.names = TRUE)
  
  png(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( "KEGG_", kegg_pathway, "_KO-WT_Heatmap.png")), 
      width = 19 * length( gene_symbols), 
      height = 600, units="px")
    draw( heatmap_ko_wt)
  dev.off()
  
  # Perform the AddModuleScore on the list of pathway genes
  sc10x.seurat = AddModuleScore( sc10x.seurat, features = list( gene_symbols), name = kegg_pathway)
  pathway_score_name = paste0( kegg_pathway, "1")

  # Plot the modules scores by cluster and condition
  pathway_df = sc10x.seurat@meta.data[ , c( "seurat_clusters", "HTO_classification", pathway_score_name)]
  pathway_df$HTO_classification = factor( pathway_df$HTO_classification, levels = c( SAMPLE_WT, SAMPLE_KO))
  print(
    ggplot( pathway_df, 
            aes_string( x="seurat_clusters", y=pathway_score_name, fill = "HTO_classification", color = "HTO_classification")) +
      geom_violin( position = "dodge") + 
      stat_summary( fun = mean, geom = "point", color = "black", position = position_dodge( width = 0.9)) +
      stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", size = 0.5, width = 0.3, position = position_dodge( width = 0.9)) +
      scale_fill_manual( values = SAMPLE_COLOR[ c( SAMPLE_WT, SAMPLE_KO)]) +
      scale_color_manual( values = SAMPLE_COLOR[ c( SAMPLE_WT, SAMPLE_KO)]) +
      geom_hline( yintercept=0, linetype="dashed", color = "black", size = 1) + 
      theme_minimal() + ggtitle( paste( "Module scores of", kegg_pathway, "KEGG pathway (", length( gene_symbols), "genes kept)\nby cluster and condition with mean and standard deviation"))
  )
  
  # Compute the comparison statistics of module score by cluster between conditions with wilcoxon test
  
  data_cluster = sc10x.seurat@meta.data[ , c( "seurat_clusters", "HTO_classification", pathway_score_name)]
  data_cluster$HTO_classification_x_seurat_clusters = paste( data_cluster$HTO_classification, data_cluster$seurat_clusters, sep = "-")
  k_test = kruskal.test( formula = as.formula( paste( pathway_score_name, " ~ HTO_classification_x_seurat_clusters")), 
                data = data_cluster )
  
  cat("<BR>Kruskal-Wallis chi-squared =", k_test$statistic, ", df =", k_test$parameter[[ "df"]] , ", p-value = ", k_test$p.value, "<BR>")
  
  pathway_stats_df = data.frame()
  for( clusterid in levels( sc10x.seurat@meta.data$seurat_clusters)){
    data_cluster = sc10x.seurat@meta.data[ which( sc10x.seurat@meta.data$seurat_clusters == clusterid), c( "seurat_clusters", "HTO_classification", pathway_score_name)]
    stats_df = wilcox_test( data = data_cluster, formula = as.formula( paste( pathway_score_name, "~ HTO_classification")))
    stats_df$seurat_clusters = clusterid
    pathway_stats_df = rbind( pathway_stats_df,
                             stats_df
    )
  }
  
  # Adjust the p-value for multitesting with BH method to get FDR
  pathway_stats_df$FDR = p.adjust( pathway_stats_df$p, method = "BH")
  pathway_stats_df$signif = sapply( pathway_stats_df$FDR, function( fdr){ if( fdr < 0.05) "*" else "" })
  pathway_stats_df = pathway_stats_df[ , c( "seurat_clusters", "group1", "group2", "n1", "n2", "p", "FDR", "signif")]
  
  # Display the statistics result table
  print(
    pathway_stats_df %>% 
      kable( caption = paste( kegg_pathway, " : Comparison between conditions by clusters with Wilcoxon test")) %>% 
      kable_styling( full_width = FALSE)
  )
  
  cat(" \n \n"); # Required for '.tabset'
}

