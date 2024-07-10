# ##################################################
# This script reads data from previous analysis
# ##################################################




# Compare expression of gene in cluster between samples
# -------------------------------------------------------

## @knitr compare_cluster_expression


# Detect the DEG between samples in each cluster
# ...............................................
sample_cluster_markergenes_df = data.frame()

Idents( sc10x.seurat) = "seurat_clusters"
for( clusterid in levels( Idents( sc10x.seurat))){
  current_sample_cluster_markergenes_df = FindMarkers( sc10x.seurat, 
                                                               ident.1 = SAMPLE_KO, 
                                                               group.by = 'HTO_classification', 
                                                               subset.ident = clusterid,
                                                               test.use = DEG_TEST_USE)
  
  current_sample_cluster_markergenes_df = current_sample_cluster_markergenes_df[ which( current_sample_cluster_markergenes_df$p_val_adj <= MARKER_GENES_ALPHA_THRESHOLD), ]
    
  if( nrow( current_sample_cluster_markergenes_df) > 0){
    current_sample_cluster_markergenes_df$gene = row.names( current_sample_cluster_markergenes_df)
    current_sample_cluster_markergenes_df$cluster = clusterid
    row.names( current_sample_cluster_markergenes_df) = NULL
    
    sample_cluster_markergenes_df = rbind( sample_cluster_markergenes_df, current_sample_cluster_markergenes_df)
  }else{
    cat("<BR><b>No marker genes found for cluster", clusterid, " at alpha error :", MARKER_GENES_ALPHA_THRESHOLD,"</b>")
  }
}

all_sample_cluster_markergenes = unique( sample_cluster_markergenes_df$gene)

write.csv( sample_cluster_markergenes_df,
           file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "clusterDEG_KOvsWT_", DEG_TEST_USE, ".csv")),
           row.names = sample_cluster_markergenes_df$gene)


# ...................................................................
# Look at the enrichment of Chosen group of genes in DEG of clusters
# ...................................................................
cat("<H3>Enrichment of marker genes of clusters in chosen groups on genes</H3>")

universe_gene = rownames( sc10x.seurat@assays$RNA@data)
global_result = data.frame()
for( group_name in names( GROUP_GENES)){
  
  group_genes_universe = intersect( GROUP_GENES[[ group_name]], universe_gene)
  for( clusterid in unique( sample_cluster_markergenes_df$cluster)){
    markers = sample_cluster_markergenes_df[ which( sample_cluster_markergenes_df$cluster == clusterid), "gene"]
    intersec = intersect( group_genes_universe, markers)
    phyper_stats = phyper( max( 0, length( intersect) - 1), 
                          length( group_genes_universe), 
                          length( universe_gene) - length( group_genes_universe),
                          length( markers),
                          lower.tail = FALSE)
    
    global_result = rbind( global_result, data.frame( group = group_name,
                                                      cluster = clusterid,
                                                      universe.size = length( universe_gene),
                                                      group.size = length( group_genes_universe),
                                                      markers.size = length( markers),
                                                      markers.in.group = length( intersec),
                                                      hypergeom.pvalue = signif( phyper_stats, 3)
                                                      ))
  }
}

global_result$hypergeom.FDR = signif( p.adjust( global_result$hypergeom.pvalue, method = "BH"), 3)

datatable( global_result, 
           colnames = c( "Group name", "Cluster ID", "# all genes", "# group", "# cluster markers", "# markers in group", "p.value", "FDR"),
           caption = "Enrichment of Markers genes of clusters in chosen groups of genes")

# ................................................................
# Visualize the difference of expression of DEG between WT and KO
# ................................................................

cat("<H3>Difference of expression of DEG between WT and KO clusters</H3>")

# Get the mean expression of cells in cluster for WT sample
all_sample_cluster_markergenes_expression_WT_df = as.data.frame( t( GetAssayData( subset( sc10x.seurat,
                                                                                       cells = names( cells_sample[ cells_sample == SAMPLE_WT]),
                                                                                       features= all_sample_cluster_markergenes), 
                                                                             assay = "RNA", slot = "data")))
all_sample_cluster_markergenes_expression_WT_df$cluster = cells_clusterid[ row.names( all_sample_cluster_markergenes_expression_WT_df)]
all_sample_cluster_markergenes_expression_clustermean_WT_df = as.data.frame( all_sample_cluster_markergenes_expression_WT_df %>% group_by( cluster) %>% summarise( across( everything(), list( mean))))
row.names( all_sample_cluster_markergenes_expression_clustermean_WT_df) = all_sample_cluster_markergenes_expression_clustermean_WT_df$cluster
all_sample_cluster_markergenes_expression_clustermean_WT_df$cluster = NULL

Heatmap( t( all_sample_cluster_markergenes_expression_clustermean_WT_df), 
         cluster_columns = FALSE,
         column_title = paste0( "Mean expression of DEG in WT (", DEG_TEST_USE," method)"),
         row_names_gp = gpar(fontsize = 6)
         )

write.csv( t( all_sample_cluster_markergenes_expression_clustermean_WT_df),
           file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "clusterDEG_KOvsWT_", DEG_TEST_USE, "_WTMeanExpression.csv")))

# Get the mean expression of cells in cluster for KO sample
all_sample_cluster_markergenes_expression_KO_df = as.data.frame( t( GetAssayData( subset( sc10x.seurat,
                                                                                       cells = names( cells_sample[ cells_sample == SAMPLE_KO]),
                                                                                       features= all_sample_cluster_markergenes), 
                                                                               assay = "RNA", slot = "data")))
all_sample_cluster_markergenes_expression_KO_df$cluster = cells_clusterid[ row.names( all_sample_cluster_markergenes_expression_KO_df)]
all_sample_cluster_markergenes_expression_clustermean_KO_df = as.data.frame( all_sample_cluster_markergenes_expression_KO_df %>%group_by( cluster) %>% summarise( across( everything(), list( mean))))
row.names( all_sample_cluster_markergenes_expression_clustermean_KO_df) = all_sample_cluster_markergenes_expression_clustermean_KO_df$cluster
all_sample_cluster_markergenes_expression_clustermean_KO_df$cluster = NULL

Heatmap( t( all_sample_cluster_markergenes_expression_clustermean_KO_df), 
         cluster_columns = FALSE,
         column_title = paste0( "Mean expression of DEG in KO (", DEG_TEST_USE," method)"),
         row_names_gp = gpar(fontsize = 6)
         )

write.csv( t( all_sample_cluster_markergenes_expression_clustermean_KO_df),
           file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "clusterDEG_KOvsWT_", DEG_TEST_USE, "_KOMeanExpression.csv")))

# Compute the difference of the mean of expression between the samples
diff_df = all_sample_cluster_markergenes_expression_clustermean_KO_df - all_sample_cluster_markergenes_expression_clustermean_WT_df

Heatmap( t( diff_df), 
         cluster_columns = FALSE,
         column_title = paste0( "Difference of mean expression of DEG (KO - WT)"),
         row_names_gp = gpar(fontsize = 6)
         )

write.csv( t( diff_df),
           file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "clusterDEG_KOvsWT_", DEG_TEST_USE, "_DiffMeanExpression.csv")))


# ....................................................................................
# Visualize the difference of expression of CHOSEN MARKER GENES between WT and KO
# ....................................................................................

if( exists( "CHOSEN_MARKER_GENES")){

  cat("<H3>Difference of expression of Chosen marker genes between WT and KO clusters</H3>")
  
  # Get the mean expression of cells in cluster for WT sample
  all_sample_cluster_chosenmarkergenes_expression_WT_df = as.data.frame( t( GetAssayData( subset( sc10x.seurat,
                                                                                            cells = names( cells_sample[ cells_sample == SAMPLE_WT]),
                                                                                            features= CHOSEN_MARKER_GENES), 
                                                                                    assay = "RNA", slot = "data")))
  all_sample_cluster_chosenmarkergenes_expression_WT_df$cluster = cells_clusterid[ row.names( all_sample_cluster_chosenmarkergenes_expression_WT_df)]
  all_sample_cluster_chosenmarkergenes_expression_clustermean_WT_df = as.data.frame( all_sample_cluster_chosenmarkergenes_expression_WT_df %>% group_by( cluster) %>% summarise( across( everything(), list( mean))))
  row.names( all_sample_cluster_chosenmarkergenes_expression_clustermean_WT_df) = all_sample_cluster_chosenmarkergenes_expression_clustermean_WT_df$cluster
  all_sample_cluster_chosenmarkergenes_expression_clustermean_WT_df$cluster = NULL
  
  Heatmap( t( all_sample_cluster_chosenmarkergenes_expression_clustermean_WT_df), 
           cluster_columns = FALSE,
           column_title = paste0( "Mean expression of Chosen Marker Genes in WT (", DEG_TEST_USE," method)"),
           row_names_gp = gpar(fontsize = 6)
  )
  
  
  # Get the mean expression of cells in cluster for KO sample
  all_sample_cluster_chosenmarkergenes_expression_KO_df = as.data.frame( t( GetAssayData( subset( sc10x.seurat,
                                                                                            cells = names( cells_sample[ cells_sample == SAMPLE_KO]),
                                                                                            features= CHOSEN_MARKER_GENES), 
                                                                                    assay = "RNA", slot = "data")))
  all_sample_cluster_chosenmarkergenes_expression_KO_df$cluster = cells_clusterid[ row.names( all_sample_cluster_chosenmarkergenes_expression_KO_df)]
  all_sample_cluster_chosenmarkergenes_expression_clustermean_KO_df = as.data.frame( all_sample_cluster_chosenmarkergenes_expression_KO_df %>%group_by( cluster) %>% summarise( across( everything(), list( mean))))
  row.names( all_sample_cluster_chosenmarkergenes_expression_clustermean_KO_df) = all_sample_cluster_chosenmarkergenes_expression_clustermean_KO_df$cluster
  all_sample_cluster_chosenmarkergenes_expression_clustermean_KO_df$cluster = NULL
  
  Heatmap( t( all_sample_cluster_chosenmarkergenes_expression_clustermean_KO_df), 
           cluster_columns = FALSE,
           column_title = paste0( "Mean expression of chosen marker genes in KO (", DEG_TEST_USE," method)"),
           row_names_gp = gpar(fontsize = 6)
  )
  
  # Compute the difference of the mean of expression between the samples
  diff_df = all_sample_cluster_chosenmarkergenes_expression_clustermean_KO_df - all_sample_cluster_chosenmarkergenes_expression_clustermean_WT_df
  
  Heatmap( t( diff_df), 
           cluster_columns = FALSE,
           column_title = paste0( "Difference of mean expression of chosen marker genes (KO - WT)"),
           row_names_gp = gpar(fontsize = 6)
  )
}
