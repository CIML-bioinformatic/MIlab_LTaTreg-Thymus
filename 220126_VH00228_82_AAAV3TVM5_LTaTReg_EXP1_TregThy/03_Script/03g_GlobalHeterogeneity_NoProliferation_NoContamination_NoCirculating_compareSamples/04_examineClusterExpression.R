# ##################################################
# This script reads data from previous analysis
# ##################################################




# Examine expression of specific genes in clusters
# -------------------------------------------------------

## @knitr examine_cluster_expression

all_cluster_marker_genes = vector()

Idents( sc10x.seurat) = "seurat_clusters"
for( module_name in names( MODULES_GENES)){
  
  cat("<H3>Genes in module", module_name, "</H3>")
  gene_set = MODULES_GENES[[ module_name]]
  
  print( 
    DoHeatmap( sc10x.seurat, features = gene_set, slot = "data") + 
      scale_fill_gradientn(colors = c( "blue", "red")) +
      ggtitle ( paste( "Module", module_name))
  )
  
  if( grepl( "Marker_Cluster_", module_name)){
    all_cluster_marker_genes = c( all_cluster_marker_genes, gene_set)
  }
}

if( length( all_cluster_marker_genes) > 0){

  
  cat("<H3>All Clusters Markers</H3>")
  
  print( 
    DoHeatmap( sc10x.seurat, features = all_cluster_marker_genes, slot = "data") + 
      scale_fill_gradientn(colors = c( "blue", "red")) +
      ggtitle (  "All Clusters Markers")
  )
  
  # Plot the heatmap with ComplexHeatmap
  # .........................................
  
  ## Select the monitored genes that are on the expression matrix
  kept_genes = intersect( all_cluster_marker_genes, row.names( sc10x.seurat[["RNA"]]@data))
  
  ## Get the expression data
  mat <- sc10x.seurat[["RNA"]]@data[ kept_genes, ] %>% as.matrix()
  
  ## Scale the rows
  mat <- t(scale(t(mat)))
  
  ## Get the clusters
  Idents( sc10x.seurat) = "seurat_clusters"
  cluster_anno<- Idents( sc10x.seurat)
  
  ## Get the min and max values for expression range, and define the colors to use
  quantil = quantile(mat, c(0.1, 0.95))
  col_fun = circlize::colorRamp2(c(quantil[ 1], quantil[2]), c("blue", "red"))
  
  # Plot the Heatmap without row clustering
  Heatmap(mat, name = "Expression (No row clustering)",  
          column_split = factor(cluster_anno),
          cluster_columns = FALSE,
          show_column_dend = FALSE,
          cluster_column_slices = TRUE,
          column_title_gp = gpar(fontsize = 8),
          column_gap = unit(0.5, "mm"),
          cluster_rows = FALSE,
          show_row_dend = FALSE,
          col = col_fun,
          row_names_gp = gpar(fontsize = 6),
          column_title_rot = 90,
          top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
          show_column_names = FALSE,
          use_raster = FALSE)
  
  # Plot the Heatmap with row clustering
  Heatmap(mat, name = "Expression (with row clustering)",  
          column_split = factor(cluster_anno),
          cluster_columns = FALSE,
          show_column_dend = FALSE,
          cluster_column_slices = TRUE,
          column_title_gp = gpar(fontsize = 8),
          column_gap = unit(0.5, "mm"),
          cluster_rows = TRUE,
          show_row_dend = TRUE,
          col = col_fun,
          row_names_gp = gpar(fontsize = 6),
          column_title_rot = 90,
          top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
          show_column_names = FALSE,
          use_raster = FALSE)
}
