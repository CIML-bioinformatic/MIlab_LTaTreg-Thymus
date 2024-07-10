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
  cat("<H3> Pathway", kegg_pathway, "</H3>")
  cat(" \n \n"); # Required for '.tabset'
  
  # Get the KEGG ID of the pathway
  kegg_pathway_id = MONITORED_KEGG_PATHWAY_LIST[[ kegg_pathway]]
  
  # Get the genes symboles for the gene EntrezID
  gene_symbols = KEGG_pathway_gene_df[ which( KEGG_pathway_gene_df$PathwayID == paste0( "path:", kegg_pathway_id)), "Symbol"]
  
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
      theme_minimal() + ggtitle( paste( "Module scores of", kegg_pathway, "KEGG pathway\nby cluster and condition with mean and standard deviation"))
  )
  
  # Compute the comparison statistics of module score by cluster between conditions with wilcoxon test
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

