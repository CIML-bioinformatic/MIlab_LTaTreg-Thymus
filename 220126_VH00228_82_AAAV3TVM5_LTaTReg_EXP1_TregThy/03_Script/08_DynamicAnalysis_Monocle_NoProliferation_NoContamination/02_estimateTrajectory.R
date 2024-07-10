# ####################################################
# This script aim to use Monocle3 pipeline
# to analyse the cell trajectory
# ####################################################

## @knitr estimate_trajectory


# Trajectory analysis

# A subtab for the clustering
# ...............................................................................
cat( "\n \n")
cat( "\n### Clusters")
cat( "\n \n")

# -- Cluster the cells
cds_sc10X_seurat = cluster_cells( cds_sc10X_seurat, resolution = CLUSTER_RESOLUTION)


# -- Plot the clusters on UMAP
p1 = plot_cells( cds_sc10X_seurat, 
            color_cells_by = "cluster", 
            group_label_size = 10,
            show_trajectory_graph = FALSE, 
            cell_size = 1) + ggtitle( paste( "Cluster distribution\n for resolution", CLUSTER_RESOLUTION))

print( p1)
cat( "<BR>")


# A subtab for the trajectory
# ...............................................................................
cat( "\n \n")
cat( "\n### Trajectory Graph")
cat( "\n \n")

# -- Learn the trajectory from clusters and expression data
cds_sc10X_seurat = learn_graph( cds_sc10X_seurat, use_partition = FALSE, verbose = FALSE)

# -- Plot the trajectory on UMAP 
print( plot_cells( cds_sc10X_seurat,
            color_cells_by = "cluster",
            label_groups_by_cluster=FALSE,
            label_leaves=FALSE,
            label_branch_points=FALSE, 
            label_principal_points = TRUE,
            graph_label_size = 4,
            cell_size = 0.5)  + ggtitle( paste( "Trajectory for\n cluster resolution", CLUSTER_RESOLUTION)))

cat( "<BR>")

modulated_genes <- graph_test( cds_sc10X_seurat, neighbor_graph = "principal_graph", cores = 4)
modulated_genes = modulated_genes[ order( modulated_genes$q_value, modulated_genes$morans_I, decreasing = c( FALSE, TRUE)),]
path_DEG_genes <- row.names( subset( modulated_genes, q_value == 0 & morans_I > 0.25))

if( length( path_DEG_genes) > 50){
  path_DEG_genes = path_DEG_genes[ 1:50]  
}


pt.matrix <- exprs( cds_sc10X_seurat)[ match( path_DEG_genes, rownames( rowData( cds_sc10X_seurat))), order( pseudotime( cds_sc10X_seurat))]
pt.matrix <- t( apply( pt.matrix, 1, function( x){ smooth.spline( x, df=3)$y}))
pt.matrix <- t( apply( pt.matrix,1, function(x){ (x-mean( x)) / sd( x)}))
rownames(pt.matrix) <- path_DEG_genes;

# Plot heatmap with K-mean for clustering
Heatmap( pt.matrix,
          name                         = "z-score",
          col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
          show_row_names               = TRUE,
          show_column_names            = FALSE,
          row_names_gp                 = gpar(fontsize = 6),
          km = 3,
          row_title_rot                = 0,
          cluster_rows                 = TRUE,
          cluster_row_slices           = FALSE,
          cluster_columns              = FALSE,
          column_title = "Expression Z-score of DEG along pseudotime clustered by K-mean"
  )

# Plot heatmap with Ward.D2 Hierarchical Clustering
Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows       = "ward.D2",
  clustering_method_columns    = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  column_title = "Expression Z-score of DEG along pseudotime clustered by Ward method"
  )

# ...............................................................................
# ...............................................................................
# For some resolution values, get some extra details
# ...............................................................................
# ...............................................................................

if( CLUSTER_RESOLUTION == 0.003){
  
  
  # A subtab for the differentially expressed genes along trajectory
  # ...............................................................................
  cat( "\n \n")
  cat( "\n### Trajectories detail  {.tabset}")
  cat( "\n \n")
  
  # Identify the nereast cell of each principal vertex
  closest_vertex = cds_sc10X_seurat@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
  closest_to = list()
  for( vertex in sort( unique( closest_vertex[ ,1]))){
    closest_to[[ paste0( "Y_", vertex)]] = row.names(closest_vertex)[ which( closest_vertex == vertex)]
  }  
  
  # Identify the vertex of degree 1 in the principal graph (source or destination)
  vertex_degrees <- degree( cds_sc10X_seurat@principal_graph$UMAP)
  vertex_degree_1 = names( vertex_degrees[ vertex_degrees == 1])
  pairs_vertex_degree_1 = combn( vertex_degree_1, 2)
  
  # Look at the genes to monitor
  filtered_monitored_genes = intersect( MONITORED_GENES, row.names( sc10x.seurat@assays$RNA@data))
  
  # Look at the cells for each condition
  Idents( sc10x.seurat) = "HTO_classification"
  WT_cells = names( Idents( sc10x.seurat))[ which( Idents( sc10x.seurat) == SAMPLE_WT)]
  KO_cells = names( Idents( sc10x.seurat))[ which( Idents( sc10x.seurat) == SAMPLE_KO)]
  
  # For each pair of source/destination node, analyse the expression of gene along the
  # shortest path between the nodes
  pairs_vertex_degree_1 = cbind( pairs_vertex_degree_1, TRAJECTORY_CUSTOM_VERTEX_DF)
  
  for( pair_index in 1:ncol( pairs_vertex_degree_1)){
    
    origin = pairs_vertex_degree_1[ 1, pair_index]
    destination = pairs_vertex_degree_1[ 2, pair_index]
    
    cat( "\n \n")
    cat( "\n#### Trajectory:", origin, "to", destination, " {.tabset}")
    cat( "\n \n")
    
    # Identify the shortest path between two principal vertex
    principal_graph = cds_sc10X_seurat@principal_graph@listData$UMAP
    shortest_path = shortest_paths( principal_graph, from = origin, to = destination, mode="all", output = "both", predecessors = TRUE)
    
    # For each principal vertex in path, retrieve the nearest cells
    path_cells = data.frame()
    for( vertex in shortest_path$vpath[[1]]){
      nearest_cells = closest_to[[ vertex]]
      path_cells = rbind( path_cells, data.frame( principal.vertex.id = vertex, cell.id = nearest_cells))
    }
    path_cells$principal.vertex.id = factor( path_cells$principal.vertex.id, levels = shortest_path$vpath[[1]])
    
    # Compute the pseudotime of cells from the nearest cells of the origin principal vertex
    cds_sc10X_seurat <- order_cells( cds_sc10X_seurat, root_pr_nodes = origin)
    path_cells$pseudotime = cds_sc10X_seurat@principal_graph_aux@listData$UMAP$pseudotime[ path_cells$cell.id]
    
    # Provide the expression of monitored genes
    path_cells_monitored_expression = cbind( path_cells, t(cbind( apply( path_cells, 1, function( row){
         cell_id = as.character( row[[ "cell.id"]])
         return( sc10x.seurat@assays$RNA@data[ filtered_monitored_genes, cell_id])
    }))))
    
    path_cells_monitored_expression = path_cells_monitored_expression[ order( path_cells_monitored_expression$pseudotime, decreasing = FALSE), ]
    
    # Provide the expression of Path DEG genes
    path_cells_DEG_expression = cbind( path_cells, t(cbind( apply( path_cells, 1, function( row){
      cell_id = as.character( row[[ "cell.id"]])
      return( sc10x.seurat@assays$RNA@data[ path_DEG_genes, cell_id])
    }))))
    
    path_cells_DEG_expression = path_cells_DEG_expression[ order( path_cells_DEG_expression$pseudotime, decreasing = FALSE), ]
    
    # 
    # sc10x.seurat = AddModuleScore( sc10x.seurat, features = GROUP_GENES, names = names( GROUP_GENES))
    
    # Provide the score of the group of genes
    group_score = sc10x.seurat@meta.data[ path_cells$cell.id, paste0( names( GROUP_GENES), "1")]
    group_score$cell.id = row.names( group_score)
    path_group_score = merge( path_cells, group_score, by = "cell.id")
    path_group_score = path_group_score[ order( path_group_score$pseudotime, decreasing = FALSE), ]
    
    # Inject the nearest principal vertex information as metadata of the Seurat object
    row.names( path_cells) = path_cells$cell.id
    meta_principal_vertex = rep(NA, length( Cells( sc10x.seurat)))
    names( meta_principal_vertex) = Cells( sc10x.seurat)
    meta_principal_vertex[ row.names( path_cells)] = paste0( "Y_", path_cells$principal.vertex.id)
    
    sc10x.seurat = AddMetaData( sc10x.seurat, metadata = meta_principal_vertex, col.name = "principal.vertex.id")
    Idents( sc10x.seurat) = "principal.vertex.id"
    print( DimPlot( sc10x.seurat) + theme( legend.position = "None")) + ggtitle( "Nearest cells to principal vertices")
    
    # Inject the pseudotime information as metadata of the Seurat object
    meta_pseudotime = rep( NA, length( Cells( sc10x.seurat)))
    names( meta_pseudotime) = Cells( sc10x.seurat)
    meta_pseudotime[ row.names( path_cells)] = as.numeric( path_cells$pseudotime)
    
    sc10x.seurat = AddMetaData( sc10x.seurat, metadata = meta_pseudotime, col.name = "pseudotime")
    print( FeaturePlot( sc10x.seurat, features = c( pseudotime = "pseudotime")) +
             scale_colour_viridis_c( option = "plasma") +
             ggtitle( "Cells pseudotime")
           )
    cat("<BR>")
    
    cat( "\n \n")
    cat( "\n##### Monitored gene expression  {.tabset}")
    cat( "\n \n")
    
    # Change the "-" character to "." in gene names because of issues with "-" in ggplot
    names( path_cells_monitored_expression) = gsub( "-", ".", names( path_cells_monitored_expression))
    
    # For each monitored gene, look at its expression along the trajectory
    for( monitored_gene in filtered_monitored_genes){
      
      cat( "\n \n")
      cat( "\n######", monitored_gene)
      cat( "\n \n")
      
      # Convert eventual "-" to "." in the gene name
      monitored_gene = gsub( "-", ".", monitored_gene)
      
      # Plot the computed statistics and expression of gene along path
      print( ggplot( data = path_cells_monitored_expression, aes_string(  x="pseudotime", y = monitored_gene)) + 
               geom_point( aes( color = pseudotime)) + geom_smooth( formula = 'y ~ x', method = "loess") +
               scale_colour_viridis_c( option = "plasma") +
               theme_minimal() + ggtitle( paste( "Expression of", monitored_gene, "along pseudotime for all cells"))
      )

      print( ggplot(  data = NULL, aes_string(  x="pseudotime", y = monitored_gene)) + 
               geom_point( data = path_cells_monitored_expression[ which( path_cells_monitored_expression$cell.id %in% WT_cells), ], color = SAMPLE_COLOR[[ SAMPLE_WT]]) +
               geom_smooth( data = path_cells_monitored_expression[ which( path_cells_monitored_expression$cell.id %in% WT_cells), ], formula = 'y ~ x', method = "loess", color = SAMPLE_COLOR[[ SAMPLE_WT]]) + 
               geom_point( data = path_cells_monitored_expression[ which( path_cells_monitored_expression$cell.id %in% KO_cells), ], color = SAMPLE_COLOR[[ SAMPLE_KO]]) +
               geom_smooth( data = path_cells_monitored_expression[ which( path_cells_monitored_expression$cell.id %in% KO_cells), ], formula = 'y ~ x', method = "loess", color = SAMPLE_COLOR[[ SAMPLE_KO]]) + 
               theme_minimal() + ggtitle( paste( "Expression of", monitored_gene, "along pseudotime for WT and KO cells"))
      )
      
      cat("<BR>")
   
    }
    
    
    cat( "\n \n")
    cat( "\n##### DEG gene along path expression  {.tabset}")
    cat( "\n \n")
    
    # Change the "-" character to "." in gene names because of issues with "-" in ggplot
    names( path_cells_expression) = gsub( "-", ".", names( path_cells_expression))
    
    # For each monitored gene, look at its expression along the trajectory
    for( monitored_gene in filtered_monitored_genes){
      
      cat( "\n \n")
      cat( "\n######", monitored_gene)
      cat( "\n \n")
      
      # Convert eventual "-" to "." in the gene name
      monitored_gene = gsub( "-", ".", monitored_gene)
      
      # Plot the computed statistics and expression of gene along path
      print( ggplot( data = path_cells_expression, aes_string(  x="pseudotime", y = monitored_gene)) + 
               geom_point( aes( color = pseudotime)) + geom_smooth( formula = 'y ~ x', method = "loess") +
               scale_colour_viridis_c( option = "plasma") +
               theme_minimal() + ggtitle( paste( "Expression of", monitored_gene, "along pseudotime for all cells"))
      )
      
      print( ggplot(  data = NULL, aes_string(  x="pseudotime", y = monitored_gene)) + 
               geom_point( data = path_cells_expression[ which( path_cells_expression$cell.id %in% WT_cells), ], color = SAMPLE_COLOR[[ SAMPLE_WT]]) +
               geom_smooth( data = path_cells_expression[ which( path_cells_expression$cell.id %in% WT_cells), ], formula = 'y ~ x', method = "loess", color = SAMPLE_COLOR[[ SAMPLE_WT]]) + 
               geom_point( data = path_cells_expression[ which( path_cells_expression$cell.id %in% KO_cells), ], color = SAMPLE_COLOR[[ SAMPLE_KO]]) +
               geom_smooth( data = path_cells_expression[ which( path_cells_expression$cell.id %in% KO_cells), ], formula = 'y ~ x', method = "loess", color = SAMPLE_COLOR[[ SAMPLE_KO]]) + 
               theme_minimal() + ggtitle( paste( "Expression of", monitored_gene, "along pseudotime for WT and KO cells"))
      )
      
      cat("<BR>")
      
    }
    
    
    
    
    cat( "\n \n")
    cat( "\n##### Chosen groups of genes score {.tabset}")
    cat( "\n \n")
    
    # For each group of genes, look at its score along the trajectory
    for( groupid in names( GROUP_GENES)){
      
      cat( "\n \n")
      cat( "\n######", groupid)
      cat( "\n \n")
      
      # Plot the scores of group of gene along path
      print( ggplot( data = path_group_score, aes_string(  x="pseudotime", y = paste0( groupid, "1"))) + 
               geom_point( aes( color = pseudotime)) + geom_smooth( formula = 'y ~ x', method = "loess") +
               scale_colour_viridis_c( option = "plasma") +
               theme_minimal() + ggtitle( paste( "Score of", groupid, "along pseudotime for all cells"))
      )
      
      print( ggplot(  data = NULL, aes_string(  x="pseudotime", y = paste0( groupid, "1"))) + 
               geom_point( data = path_group_score[ which( path_group_score$cell.id %in% WT_cells), ], color = SAMPLE_COLOR[[ SAMPLE_WT]]) +
               geom_smooth( data = path_group_score[ which( path_group_score$cell.id %in% WT_cells), ], formula = 'y ~ x', method = "loess", color = SAMPLE_COLOR[[ SAMPLE_WT]]) + 
               geom_point( data = path_group_score[ which( path_group_score$cell.id %in% KO_cells), ], color = SAMPLE_COLOR[[ SAMPLE_KO]]) +
               geom_smooth( data = path_group_score[ which( path_group_score$cell.id %in% KO_cells), ], formula = 'y ~ x', method = "loess", color = SAMPLE_COLOR[[ SAMPLE_KO]]) + 
               theme_minimal() + ggtitle( paste( "Score of", groupid, "along pseudotime for WT and KO cells"))
      )
      
      cat("<BR>")
      
    }
  }
}
