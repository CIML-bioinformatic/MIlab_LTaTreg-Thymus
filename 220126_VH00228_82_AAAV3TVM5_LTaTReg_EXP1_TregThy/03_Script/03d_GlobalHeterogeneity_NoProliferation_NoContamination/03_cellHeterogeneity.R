# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################




# PCA
# ---

## @knitr heterogeneity_pca

# Compute PCA on selected variable genes
nbPC=PCA_NPC
if(PCA_NPC>length(Cells(sc10x)))
{
  warning( paste0( "Number of cells in object (", length(Cells(sc10x)), ") smaller than requested number of PCs (", PCA_NPC,"), setting lower PC number..." ))
  nbPC = length(Cells(sc10x))
}           
sc10x <- RunPCA( object   = sc10x,
                 features = VariableFeatures( sc10x),
                 npcs     = nbPC,
                 verbose  = .VERBOSE);

# Identify clusters of cells by graph approach
nbPC_findclusters=FINDCLUSTERS_USE_PCA_NBDIMS
if(FINDCLUSTERS_USE_PCA_NBDIMS>nbPC)
{
  warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'findclusters' (", FINDCLUSTERS_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
  nbPC_findclusters = nbPC
}  
sc10x <- FindNeighbors(object    = sc10x,
                       k.param   = FINDNEIGHBORS_K, 
                       reduction = "pca",
                       dims      = 1:nbPC_findclusters,
                       verbose   = .VERBOSE);

# Plot PCA, highlighting seurat clusters for combination of dimensions
invisible( apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims)
{
  print( DimPlot( object=sc10x, reduction="pca", dims = dims, group.by = "orig.ident") +
           theme( legend.position = "none"));                               # Remove legend in this case...
}));




## @knitr heterogeneity_pca_umisCounts

# Same PCA plots but highlight UMI counts
invisible( apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims)
{
  print( FeaturePlot( sc10x, feature = "nCount_RNA", reduction = "pca", dims = dims) +
           theme( legend.position = "none",
                  plot.margin = margin(0, 0, 0, 0, "cm"),
                  plot.title = element_blank()));
}));




## @knitr heterogeneity_pca_genesCounts

# Same PCA plots but highlight feature counts
invisible( apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims)
{
  print( FeaturePlot( sc10x, feature = "nFeature_RNA", reduction = "pca", dims = dims) +
           theme( legend.position = "none",
                  plot.margin = margin(0, 0, 0, 0, "cm"),
                  plot.title = element_blank()));
}));




## @knitr heterogeneity_pca_correlations

# Isolate the PCA location of cells in first dimensions, together with UMI and genes counts for correlation analysis (makes use of cbind recycling to repeat values for each stacked PC)
relationToPC = suppressWarnings( cbind( stack( as.data.frame( Embeddings( sc10x, reduction = "pca")[ rownames( sc10x[[]]), paste0( "PC_", 1:PCA_PLOTS_NBDIMS) ])),
                                        sc10x[[ c( "nCount_RNA", "nFeature_RNA") ]],
                                        Cluster = Idents( sc10x)));

# Plot relationship of UMIs and genes counts with PCs (combine plots using '/' from 'patchwork' lib)
print( (ggplot( data = relationToPC, aes( x = values, y = nCount_RNA)) +
          facet_wrap( ~ind) +
          stat_binhex( bins = 60) +
          #geom_point( aes(col = Cluster), alpha = 0.5) +
          geom_smooth( method = 'lm') +
          stat_cor( method = "spearman") +
          ylab( "# UMIs") +
          theme( axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 plot.margin = unit( c( 1, 1, -0.5, 0.5), "lines")))
       /
         (ggplot( data = relationToPC, aes( x = values, y = nFeature_RNA)) +
            facet_wrap( ~ind) +
            stat_binhex( bins = 60) +
            #geom_point( aes(col = Cluster), alpha = 0.5) +
            geom_smooth( method = 'lm') +
            stat_cor( method = "spearman") +
            xlab( "PC values") +
            ylab( "# Genes")));




## @knitr heterogeneity_pca_loadings

# Plot PCA loadings
invisible( apply( combn( 1:PCA_PLOTS_NBDIMS, 2), 2, function(dims)
{
  namesPC=paste0( "PC_", dims);
  # Get the loading values for concerned PCs
  loadingsMatrix = Loadings( sc10x, reduction = "pca")[ , namesPC ];
  # Sort features by average absolute value and convert as DF with features names as column
  loadingsMatrix = head( loadingsMatrix[ order( apply( loadingsMatrix, 1, function(x){ mean( abs( x)) }), decreasing = TRUE), ], PCA_PLOTS_NBFEATURES);
  loadingsDF = data.frame( loadingsMatrix, features = rownames( loadingsMatrix));
  
  # Define symmetric and consistent axes for group of plots
  axesLimit = max( abs( loadingsMatrix));
  
  # Plot arrows and features name
  print( ggplot( data = loadingsDF, aes_( x = as.name( namesPC[1]), y = as.name( namesPC[2]))) +
           coord_cartesian( xlim = c( -axesLimit, axesLimit), ylim = c( -axesLimit, axesLimit)) +
           geom_text_repel( aes( label = features), max.iter = 10000) +
           geom_segment( x = 0 , y = 0, aes_( xend = as.name(namesPC[1]), yend = as.name(namesPC[2])), col = "#00000044", arrow = arrow( length = unit( 2, "mm"))));
}));




# DIMENSIONAL REDUCTION (TSNE/UMAP)
###################################

## @knitr heterogeneity_dimReduc
nbPC_dimreduc=DIMREDUC_USE_PCA_NBDIMS
if(DIMREDUC_USE_PCA_NBDIMS>nbPC)
{
  warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'dimreduc' (", DIMREDUC_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
  nbPC_dimreduc = nbPC
}

sc10x = RunUMAP( sc10x, dims = 1:nbPC_dimreduc);
sc10x = RunTSNE( sc10x, dims = 1:nbPC_dimreduc);

# Save resulting coordinates for all cells as 'tsv' files
write.table( Embeddings(sc10x, reduction = "umap"), 
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "cellsCoordinates_umap.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

write.table( Embeddings(sc10x, reduction = "tsne"), 
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "cellsCoordinates_tsne.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");



# IDENTIFY CLUSTERS
###################

## @knitr heterogeneity_identifyClusters

sc10x <- FindClusters(object             = sc10x,
                      resolution         = FINDCLUSTERS_RESOLUTION,
                      algorithm          = FINDCLUSTERS_ALGORITHM,
                      temp.file.location = "/tmp/",
                      verbose            = .VERBOSE);

#DimPlot(sc10x)

# Rename the clusters if required
if( exists( "RENAME_CLUSTERS")){
  if( SAMPLE_NAME %in% names( RENAME_CLUSTERS)){
    Idents( sc10x) = "seurat_clusters"
    sc10x = AddMetaData( sc10x, metadata = Idents( sc10x), col.name = "original_seurat_clusters")
    new_names = Idents( sc10x)
    for( cluster_id in levels( Idents( sc10x))){
      new_names[ which( Idents( sc10x) == cluster_id)] = RENAME_CLUSTERS[[ SAMPLE_NAME]][[ cluster_id]]
    }
    sc10x = AddMetaData( sc10x, metadata = new_names, col.name = "seurat_clusters")
    
  }
}
Idents( sc10x) = "seurat_clusters"
Idents( sc10x.all) = "seurat_clusters"

# Add the new clusters to the global Seurat Cluster object
sc10x.all = AddMetaData( sc10x.all, metadata = Idents( sc10x)[ Cells( sc10x.all)], col.name = paste0( CLUSTER_GROUP_ID, "_seurat_clusters"))


# Show number of cells in each cluster
clustersCount = as.data.frame( table( Cluster = sc10x[[ "seurat_clusters" ]]), responseName = "CellCount");

# Define a set of colors for clusters (based on ggplot default)
clustersColor = hue_pal()( nlevels( Idents( sc10x)));
names( clustersColor) = levels( Idents( sc10x));

clustersColorAll = hue_pal()( nlevels( Idents( sc10x.all)));
names( clustersColorAll) = levels( Idents( sc10x.all));

# Save cells cluster identity as determined with 'FindClusters'
write.table( data.frame(sc10x[["numID"]], identity = Idents(sc10x)), 
             file = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "cellsClusterIdentity.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

# Also save cluster color attribution for reference
# Save cells cluster identity as determined with 'FindClusters'
write.table( clustersColor, 
             file = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "clustersColor.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = FALSE,
             sep="\t");


# Create datatable
datatable( clustersCount,
           class = "compact",
           rownames = FALSE,
           colnames = c("Cluster", "Nb. Cells"),
           options = list(dom = "<'row'rt>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          columnDefs = list( # Center all columns
                                            list( targets = 0:(ncol(clustersCount)-1),
                                            className = 'dt-center')),
                          orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                          paging = FALSE, # Disable pagination (show all)
                          processing = TRUE,
                          scrollCollapse = TRUE,
                          scroller = TRUE,  # Only load visible data
                          scrollX = TRUE,
                          scrollY = "525px",
                          stateSave = TRUE)) %>%
  # Add color from cluster
  formatStyle( columns = "Cluster",
               color = styleEqual( names(clustersColor), clustersColor),
               fontWeight = 'bold');


# Compute a matrix of average expression value by cluster (for each gene)
geneExpByCluster = do.call( rbind, 
                            apply( as.matrix( GetAssayData( sc10x)), # expression values
                                   1,                                # by rows
                                   tapply,                           # apply by group
                                   INDEX = Idents( sc10x),           # clusters IDs
                                   mean,                             # summary function
                                   simplify = FALSE));               # do not auto coerce

# Save it as 'tsv' file
write.table( geneExpByCluster, 
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "normExpressionByCluster.tsv")), 
             quote = FALSE, 
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");





## @knitr heterogeneity_dimReduc_UMAPPlot

# Plot UMAP of new clusters
DimPlot( sc10x, reduction = "umap", label=TRUE, label.size = 10) +
  scale_color_manual( name = "Cluster", values = clustersColor) +  # Fix legend title + set colors
  ggtitle( 'UMAP of new cluster annotation') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         legend.position = "none",
         plot.margin = margin( 0, 0, 0, 0, "cm"),
         plot.title = element_text( face = "bold",
                                    size = rel( 16/14),
                                    hjust = 0.5,
                                    vjust = 1,
                                    margin = margin( b = 7)));

# Plot UMAP of original cluster annotation
Idents( sc10x) = "original_seurat_clusters"
DimPlot( sc10x, reduction = "umap", label=TRUE, label.size = 10) +
  scale_color_manual( name = "Cluster", values = clustersColorAll) +  # Fix legend title + set colors
  ggtitle( 'UMAP of original cluster annotation') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         legend.position = "none",
         plot.margin = margin( 0, 0, 0, 0, "cm"),
         plot.title = element_text( face = "bold",
                                    size = rel( 16/14),
                                    hjust = 0.5,
                                    vjust = 1,
                                    margin = margin( b = 7)));

# Plot UMAP of new cluster annotation in original Seurat Object
Idents( sc10x.all) = paste0( CLUSTER_GROUP_ID, "_seurat_clusters")
DimPlot( sc10x.all, reduction = "umap", label=TRUE, label.size = 10) +
  scale_color_manual( name = "Cluster", values = clustersColor) +  # Fix legend title + set colors
  ggtitle( 'UMAP of new cluster group annotation\nin original embedding') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         legend.position = "none",
         plot.margin = margin( 0, 0, 0, 0, "cm"),
         plot.title = element_text( face = "bold",
                                    size = rel( 16/14),
                                    hjust = 0.5,
                                    vjust = 1,
                                    margin = margin( b = 7)));


# Identify the cells of condition WT
Idents( sc10x) = "HTO_classification"
WT_cells_names = names( Idents( sc10x))[ which( Idents( sc10x) == SAMPLE_WT)]
Idents( sc10x) = "seurat_clusters"

# Plot UMAP of WT cells
DimPlot( subset( sc10x, cells = WT_cells_names), reduction = "umap", label=TRUE, label.size = 10) +
  scale_color_manual( name = "Cluster", values = clustersColor) +  # Fix legend title + set colors
  ggtitle( 'UMAP of WT Cells') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         legend.position = "none",
         plot.margin = margin( 0, 0, 0, 0, "cm"),
         plot.title = element_text( face = "bold",
                                    size = rel( 16/14),
                                    hjust = 0.5,
                                    vjust = 1,
                                    margin = margin( b = 7)));

# Identify the cells of condition KO 
Idents( sc10x) = "HTO_classification"
KO_cells_names = names( Idents( sc10x))[ which( Idents( sc10x) == SAMPLE_KO)]
Idents( sc10x) = "seurat_clusters"

# Plot UMAP of KO cells
DimPlot( subset( sc10x, cells = KO_cells_names), reduction = "umap", label=TRUE, label.size = 10) +
  scale_color_manual( name = "Cluster", values = clustersColor) +  # Fix legend title + set colors
  ggtitle( 'UMAP of KO Cells') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         legend.position = "none",
         plot.margin = margin( 0, 0, 0, 0, "cm"),
         plot.title = element_text( face = "bold",
                                    size = rel( 16/14),
                                    hjust = 0.5,
                                    vjust = 1,
                                    margin = margin( b = 7)));
  
# Look at the distribution of cells along clusters
cluster_vs_sample_df = as.data.frame.matrix( table( sc10x@meta.data$seurat_clusters, sc10x@meta.data$HTO_classification))
cluster_vs_sample_df %>% kable( caption = "Cluster versus Sample distribution") %>% kable_styling()
signif( chisq.test( cluster_vs_sample_df)$residuals, 3) %>% kable( caption = "Cluster versus Sample chi2 Pearson residuals") %>% kable_styling()



## @knitr heterogeneity_dimReduc_thumbnail
# Non-interactive with large labels for generating thumbnails when analyzing monitored and marker genes expression
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)
# If contains several values, loop on them so we don't have to call chunk separately (which creates a new <p>aragraph in html)

# Replot the projection with colored clusters and large labels (add title if several plotted in the loop)
reductionVector = if(exists("useReduction")) useReduction else "umap";
for(currentReduction in reductionVector)
{
  ggFigure = suppressWarnings( DimPlot( sc10x, reduction = currentReduction, label=TRUE, label.size = 10)) +
                                 scale_color_manual( name = "Cluster", values = clustersColor) +  # Fix legend title + set colors
                                 theme( axis.title.x = element_blank(),
                                        axis.title.y = element_blank(),
                                        legend.position = "none",
                                        plot.margin = margin( 0, 0, 0, 0, "cm"),
                                        plot.title = element_text( face = "bold",
                                                                   size = rel( 16/14),
                                                                   hjust = 0.5,
                                                                   vjust = 1,
                                                                   margin = margin( b = 7)));

  if( length( reductionVector) > 1) ggFigure = ggFigure + ggtitle( label = currentReduction);
  print( ggFigure);
}




# MARKER GENES
##############

## @knitr heterogeneity_markerGenes

# Identify marker genes
markers = FindAllMarkers( object          = sc10x,
                          test.use        = FINDMARKERS_METHOD,
                          only.pos        = FINDMARKERS_ONLYPOS,
                          min.pct         = FINDMARKERS_MINPCT,
                          logfc.threshold = FINDMARKERS_LOGFC_THR,
                          return.thresh   = FINDMARKERS_PVAL_THR,
                          random.seed     = SEED,
                          verbose         = .VERBOSE);

# Save markers list as 'tsv' table
write.table( markers,
             file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "MarkerGenes.tsv")),
             quote = FALSE,
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

# Filter markers by cluster (TODO: check if downstream code works when no markers found)
topMarkers = by( markers, markers[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_log2FC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( FINDMARKERS_SHOWTOP)) x else head( x, n = FINDMARKERS_SHOWTOP));
});

# Merge marker genes in a single data.frame and render it as datatable
topMarkersDF = do.call( rbind, topMarkers);
# Select and order columns to be shown in datatable
topMarkersDT = topMarkersDF[c("gene", "cluster", "avg_log2FC", "p_val_adj")]




## @knitr heterogeneity_markerGenes_table

# Create datatable
datatable( topMarkersDT,
           class = "compact",
           filter="top",
           rownames = FALSE,
           colnames = c("Gene", "Cluster", "Avg. LogFC", "Adj. Pvalue"),
           caption = paste(ifelse( is.null( FINDMARKERS_SHOWTOP), "All", paste("Top", FINDMARKERS_SHOWTOP)), "marker genes for each cluster"),
           extensions = c('Buttons', 'Select'),
           options = list(dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          buttons = exportButtonsListDT,
                          columnDefs = list(
                            list( # Center all columns except first one
                              targets = 1:(ncol( topMarkersDT)-1),
                              className = 'dt-center'),
                            list( # Set renderer function for 'float' type columns (LogFC)
                              targets = ncol( topMarkersDT)-2,
                              render = htmlwidgets::JS("function ( data, type, row ) {return type === 'export' ? data : data.toFixed(4);}")),
                            list( # Set renderer function for 'scientific' type columns (PValue)
                              targets = ncol( topMarkersDT)-1,
                              render = htmlwidgets::JS( "function ( data, type, row ) {return type === 'export' ? data : data.toExponential(4);}"))),
                          #fixedColumns = TRUE, # Does not work well with filter on this column
                          #fixedHeader = TRUE, # Does not work well with 'scrollX'
                          lengthMenu = list(c( 10, 50, 100, -1),
                                            c( 10, 50, 100, "All")),
                          orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                          processing = TRUE,
                          #search.regex= TRUE, # Does not work well with 'search.smart'
                          search.smart = TRUE,
                          select = TRUE, # Enable ability to select rows
                          scrollCollapse = TRUE,
                          scroller = TRUE,  # Only load visible data
                          scrollX = TRUE,
                          scrollY = "525px",
                          stateSave = TRUE)) %>%
  # Add bar relative to logFC
  formatStyle( columns = "avg_log2FC",
               background = styleColorBar( data = range( topMarkersDT[["avg_log2FC"]]), 'lightblue', angle = -90),
               backgroundSize = '95% 50%',      # Set horizontal and vertical span in cell
               backgroundRepeat = 'no-repeat',
               backgroundPosition = 'center') %>%
  # Add color from cluster
  formatStyle( columns = "cluster",
               backgroundColor = styleEqual( names(clustersColor),
                                             scales::alpha(clustersColor, 0.3)));



## @knitr heterogeneity_markerGenes_heatmap_mean

# Get the matrix of expression and associated clusters from Seurat object
# ........................................................................
expMat = as.matrix( GetAssayData( sc10x));
Idents( sc10x) = "seurat_clusters"
clusterID = Idents( sc10x);

# Select marker genes and reorder cells to group clusters together
topMarkersGenes = topMarkersDF[["gene"]];
clusterOrdering = order( clusterID);

expMat = expMat[topMarkersGenes, clusterOrdering];
clusterID = clusterID[clusterOrdering];

# Compute the mean of top markers in clusters and produce matrix with the result
# ........................................................................
all_mean_expression = vector()
cluster_set = levels( clusterID)
clusters_marker_genes = topMarkersDF[ which( topMarkersDF$cluster %in% cluster_set), "gene"]
for( cluster_id in cluster_set){
  clusters_cells = names( clusterID)[ which( clusterID == cluster_id)]
  all_mean_expression = c( all_mean_expression, BiocGenerics::rowMeans( expMat[ clusters_marker_genes, clusters_cells]))
}
meanExpMat = t( matrix( all_mean_expression, byrow = TRUE, ncol = FINDMARKERS_SHOWTOP*length( cluster_set)))
colnames( meanExpMat) = cluster_set
rownames( meanExpMat) = paste0( topMarkersDF[ which( topMarkersDF$cluster %in% cluster_set), "cluster"], '.', clusters_marker_genes)

# Plot a heatmap of the matrix of mean expression of top marker genes in clusters
# ............................................................................


cat("\n \n")
pheatmap( expMat,
          color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(80),
          breaks = seq( -4.2, 4.2, 0.1),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          scale = "row",
          annotation_row = data.frame(Markers = topMarkersDF[ rownames( meanExpMat), "cluster"], stringsAsFactors = FALSE, row.names = rownames( meanExpMat) ),
          annotation_col = data.frame(Cluster = clusterID, stringsAsFactors = FALSE, row.names = colnames( expMat)),
          annotation_colors = list( Markers = clustersColor,
                                    Cluster = clustersColor),
          show_colnames = FALSE,
          fontsize_row = 5,
          main = "Scaled normalized expression by cell\nof top markers");

cat("\n \n")
DotPlot( sc10x, 
         features = topMarkersGenes,
         cols = "RdBu") + 
  theme( axis.text.x = element_text( angle = 45, hjust=1)) +
  ggtitle( "Scaled mean normalized expression\nby cluster of top markers")

cat("\n \n")
pheatmap( meanExpMat,
          color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          scale = "row",
          annotation_row = data.frame(Markers = topMarkersDF[ rownames( meanExpMat), "cluster"], stringsAsFactors = FALSE, row.names = rownames( meanExpMat) ),
          annotation_col = data.frame(Cluster = cluster_set, stringsAsFactors = FALSE, row.names = cluster_set),
          annotation_colors = list( Markers = clustersColor,
                                    Cluster = clustersColor),
          show_colnames = TRUE,
          fontsize_row = 5,
          main = "Scaled mean normalized expression\nby cluster of top markers");


## @knitr heterogeneity_markerGenes_expression_projection
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Plot expression values of marker genes on dimreduc figures for each cluster (TODO: message if list empty)
invisible( lapply( names( topMarkers), function(clusterName)
{
  cat("##### Cl. <span style='border-radius: 3px; border: 3px solid ", clustersColor[clusterName], "; padding:0px 2px'>", clusterName, "</span>\n");

  # Highlight cells of current cluster on a dimreduc plot
  highlightClusterPlot(clusterName, seuratObject = sc10x, reduction = ifelse( exists("useReduction"), useReduction, "umap"));

  # Plots expression on projected cells
  invisible( lapply( topMarkers[[clusterName]][["gene"]], function(featureName)
    {
      print( FeaturePlot( sc10x, features = featureName, reduction = ifelse( exists("useReduction"), useReduction, "umap"), order = TRUE) +
               theme( axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      legend.position = "none"));
    }));

  cat(" \n \n"); # Required for '.tabset'
}));




## @knitr heterogeneity_markerGenes_expression_violin

# Plot expression values of marker genes as violinplot for each cluster (TODO: message if list empty)
invisible( lapply( names( topMarkers), function(clusterName)
{
  cat("##### Cl. <span style='border-radius: 3px; border: 3px solid ", clustersColor[clusterName], "; padding:0px 2px'>", clusterName, "</span>\n");

  # Remind cluster name in an empty figure to keep consistent alignment of panels between tabs
  plot( c( 0, 1), c( 0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n');
  text( x = 0.5, y = 0.5, paste( "Cluster", clusterName), cex = 2, col = clustersColor[clusterName]);

  # Violinplot for expression value of marker genes by cluster (+ number of 'zero' and 'not zero' cells)
  invisible( lapply( topMarkers[[clusterName]][["gene"]], violinFeatureByCluster, seuratObject = sc10x, clustersColor = clustersColor));

  cat(" \n \n"); # Required for '.tabset'
}));




# CLUSTER FUNCTIONAL ENRICHMENT
################################

## @knitr heterogeneity_markerGenes_functional_enrichment

for( clusterName in levels( markers$cluster)){

  # Add a sub-tab with the Cluster name
  # ...............................................................................
  cat( "\n \n")
  cat( "\n##### Cluster", clusterName, "{.tabset}")
  cat( "\n \n")
  
  
  # Add a sub-tab with the enrichment of the positive markers in the GO BP terms
  # ...............................................................................
  cat( "\n \n")
  cat( "\n###### Markers genes GO Biological Process")
  cat( "\n \n")
  
  positive_ego <- enrichGO(gene          = markers[ which( markers$cluster == clusterName), "gene"],
                           universe      = rownames( sc10x@assays$RNA),
                           OrgDb         = org.Mm.eg.db,
                           keyType       = "SYMBOL",
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = ENRICHMENT_GO_PVALUECUTOFF,
                           qvalueCutoff  = ENRICHMENT_GO_QVALUECUTOFF)
  
  selected_term_number = length( which( positive_ego@result$p.adjust <= positive_ego@pvalueCutoff & positive_ego@result$qvalue <= positive_ego@qvalueCutoff))
  
  if( selected_term_number> 0){  
    
    print( dotplot( positive_ego, showCategory=10) + ggtitle("Enrichment of marker genes in BP GO terms"))
    positive_ego_df = data.frame( positive_ego)
    
    positive_ego_simMatrix <- calculateSimMatrix( positive_ego_df$ID,
                                                  orgdb="org.Mm.eg.db",
                                                  ont="BP",
                                                  method="Rel")
    
    positive_ego_scores <- setNames(-log10( positive_ego_df$qvalue), positive_ego_df$ID)
    positive_ego_reducedTerms <- reduceSimMatrix(positive_ego_simMatrix,
                                                 positive_ego_scores,
                                                 threshold=0.7,
                                                 orgdb="org.Mm.eg.db")
    
    print( scatterPlot( positive_ego_simMatrix, positive_ego_reducedTerms))
    
    treemapPlot( positive_ego_reducedTerms)
  }
  else{
    cat("<BR>No result<BR>")
  }
  
  
  # Add a sub-tab with the enrichment of the positive markers in the GO MF terms
  # ...............................................................................
  cat( "\n \n")
  cat( "\n###### Markers genes GO Molecular Function")
  cat( "\n \n")
  
  positive_ego <- enrichGO(gene          = markers[ which( markers$cluster == clusterName), "gene"],
                           universe      = rownames( sc10x@assays$RNA),
                           OrgDb         = org.Mm.eg.db,
                           keyType       = "SYMBOL",
                           ont           = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = ENRICHMENT_GO_PVALUECUTOFF,
                           qvalueCutoff  = ENRICHMENT_GO_QVALUECUTOFF)
  
  selected_term_number = length( which( positive_ego@result$p.adjust <= positive_ego@pvalueCutoff & positive_ego@result$qvalue <= positive_ego@qvalueCutoff))
  
  if( selected_term_number> 0){  
    
    print( dotplot( positive_ego, showCategory=10) + ggtitle("Enrichment of positive markers in MF GO terms"))
    positive_ego_df = data.frame( positive_ego)
    
    positive_ego_simMatrix <- calculateSimMatrix( positive_ego_df$ID,
                                                  orgdb="org.Mm.eg.db",
                                                  ont="MF",
                                                  method="Rel")
    
    positive_ego_scores <- setNames(-log10( positive_ego_df$qvalue), positive_ego_df$ID)
    positive_ego_reducedTerms <- reduceSimMatrix(positive_ego_simMatrix,
                                                 positive_ego_scores,
                                                 threshold=0.7,
                                                 orgdb="org.Mm.eg.db")
    
    print( scatterPlot( positive_ego_simMatrix, positive_ego_reducedTerms))
    
    treemapPlot( positive_ego_reducedTerms)
  }
  else{
    cat("<BR>No result<BR>")
  }
}


# GROUP GENES
#################

## @knitr heterogeneity_chosenGroupGenes_heatmap_mean

# Compute the mean of chosen group of genes in clusters and produce matrix with the result
# ........................................................................

for( groupid in names( GROUP_GENES)){ 
  
  cat("\n \n")
  cat("#### ", groupid, "\n");
  
  # Get the genes of the current group
  group_genes = GROUP_GENES[[ groupid]]
  
  # Get the matrix of counts
  expMat = as.matrix( GetAssayData( sc10x));
  
  # Get the list of genes to analyse
  group_genes = intersect( group_genes, rownames( expMat))
  
  # Get the condition of the cells
  Idents( sc10x) = "HTO_classification"
  hto_classification = Idents( sc10x);
  
  # Get the cluster ID of the cells
  Idents( sc10x) = "seurat_clusters"
  clusterID = Idents( sc10x);
  
  # Reorder cells to group clusters together
  clusterOrdering = order( clusterID);
  clusterID = clusterID[clusterOrdering];
  hto_classification = hto_classification[ clusterOrdering]
  
  # Get the count matrix with the right genes and the right cell order
  expMat = expMat[ group_genes, names( clusterID)];
  
  # Compute the mean of cells by gene and cluster for all cells and build a matrix with the results
  # The genes with only zeros are removed to avoid issue during scaling in heatmap
  all_mean_expression = vector()
  cluster_set = levels( clusterID)
  for( cluster_id in cluster_set){
    clusters_cells = names( clusterID)[ which( clusterID == cluster_id)]
    all_mean_expression = c( all_mean_expression, BiocGenerics::rowMeans( expMat[ group_genes, clusters_cells]))
  }
  meanExpMat = t( matrix( all_mean_expression, byrow = TRUE, ncol = length( group_genes)))
  colnames( meanExpMat) = cluster_set
  rownames( meanExpMat) = group_genes
  meanExpMat = meanExpMat[ which( apply( meanExpMat, 1, sum) >0), ]
  
  # Compute the mean of cells by gene and cluster for WT cells only and build a matrix with the results
  # The genes with only zeros are removed to avoid issue during scaling in heatmap
  all_mean_expression = vector()
  cluster_set = levels( clusterID)
  for( cluster_id in cluster_set){
    clusters_cells = names( clusterID)[ which( clusterID == cluster_id & hto_classification == SAMPLE_WT)]
    all_mean_expression = c( all_mean_expression, BiocGenerics::rowMeans( expMat[ group_genes, clusters_cells]))
  }
  meanExpMat_WT = t( matrix( all_mean_expression, byrow = TRUE, ncol = length( group_genes)))
  colnames( meanExpMat_WT) = cluster_set
  rownames( meanExpMat_WT) = group_genes
  meanExpMat_WT = meanExpMat_WT[ which( apply( meanExpMat_WT, 1, sum) >0), ]
  
  # Compute the mean of cells by gene and cluster for WT cells only and build a matrix with the results
  # The genes with only zeros are removed to avoid issue during scaling in heatmap
  all_mean_expression = vector()
  cluster_set = levels( clusterID)
  for( cluster_id in cluster_set){
    clusters_cells = names( clusterID)[ which( clusterID == cluster_id & hto_classification == SAMPLE_KO)]
    all_mean_expression = c( all_mean_expression, BiocGenerics::rowMeans( expMat[ group_genes, clusters_cells]))
  }
  meanExpMat_KO = t( matrix( all_mean_expression, byrow = TRUE, ncol = length( group_genes)))
  colnames( meanExpMat_KO) = cluster_set
  rownames( meanExpMat_KO) = group_genes
  meanExpMat_KO = meanExpMat_KO[ which( apply( meanExpMat_KO, 1, sum) >0), ]

  # Plot the heatmap of mean counts
  # ................................
  
  cat("<p style='color:red;'><b>WARNING: the genes with only zero values along clusters are removed from the heatmap</b></p>")
  
  pheatmap( meanExpMat,
            color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
            cluster_rows = TRUE,
            cluster_cols = FALSE,
            scale = "row",
            annotation_col = data.frame(Cluster = cluster_set, stringsAsFactors = FALSE, row.names = cluster_set),
            annotation_colors = list( Cluster = clustersColor),
            show_colnames = TRUE,
            fontsize_row = 10,
            border_color = NA, 
            main = paste( "Scaled mean expression of genes in clusters \nfor group of genes:", groupid));
  
  pheatmap( meanExpMat_WT,
            color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
            cluster_rows = TRUE,
            cluster_cols = FALSE,
            scale = "row",
            annotation_col = data.frame(Cluster = cluster_set, stringsAsFactors = FALSE, row.names = cluster_set),
            annotation_colors = list( Cluster = clustersColor),
            show_colnames = TRUE,
            fontsize_row = 10,
            border_color = NA, 
            main = paste( "Scaled mean expression of genes in clusters \non WT cells only for group of genes:", groupid));
  
  pheatmap( meanExpMat_KO,
            color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
            cluster_rows = TRUE,
            cluster_cols = FALSE,
            scale = "row",
            annotation_col = data.frame(Cluster = cluster_set, stringsAsFactors = FALSE, row.names = cluster_set),
            annotation_colors = list( Cluster = clustersColor),
            show_colnames = TRUE,
            fontsize_row = 10,
            border_color = NA, 
            main = paste( "Scaled mean expression of genes in clusters \non KO cells only for group of genes:", groupid));
  
  
  # Compute the Module Score of group
  # ..................................
  
  sc10x = AddModuleScore( sc10x, features = list( group_genes), name = groupid)
  module_name = paste0( groupid, "1")
  
  module_df = sc10x@meta.data[ , c( "seurat_clusters", "HTO_classification", module_name)]
  module_df$HTO_classification = factor( module_df$HTO_classification, levels = c( SAMPLE_WT, SAMPLE_KO))
  
  # Plot the modules scores by cluster and condition
  print(
    ggplot( module_df, 
            aes_string( x="seurat_clusters", y=module_name, fill = "HTO_classification", color = "HTO_classification")) +
      geom_violin( position = "dodge") + 
      stat_summary( fun = mean, geom = "point", color = "black", position = position_dodge( width = 0.9)) +
      stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", size = 0.5, width = 0.3, position = position_dodge( width = 0.9)) +
      scale_fill_manual( values = SAMPLE_COLOR[ c( SAMPLE_WT, SAMPLE_KO)]) +
      scale_color_manual( values = SAMPLE_COLOR[ c( SAMPLE_WT, SAMPLE_KO)]) +
      theme_minimal() + ggtitle( paste( "Module scores of", groupid, "by cluster and condition", "\nwith mean and standard deviation"))
  )
  
  # Compute the comparison statistics of module score by cluster between conditions with wilcoxon test
  module_stats_df = data.frame()
  for( clusterid in levels( sc10x@meta.data$seurat_clusters)){
    data_cluster = sc10x@meta.data[ which( sc10x@meta.data$seurat_clusters == clusterid), c( "seurat_clusters", "HTO_classification", module_name)]
    stats_df = wilcox_test( data = data_cluster, formula = as.formula( paste( module_name, "~ HTO_classification")))
    stats_df$seurat_clusters = clusterid
    module_stats_df = rbind( module_stats_df,
                             stats_df
                            )
  }
  
  # Adjust the p-value for multitesting with BH method to get FDR
  module_stats_df$FDR = p.adjust( module_stats_df$p, method = "BH")
  module_stats_df$signif = sapply( module_stats_df$FDR, function( fdr){ if( fdr < 0.05) "*" else "" })
  module_stats_df = module_stats_df[ , c( "seurat_clusters", "group1", "group2", "n1", "n2", "p", "FDR", "signif")]
  
  # Display the statistics result table
  print(
    module_stats_df %>% 
      kable( caption = paste( groupid, " : Comparison between conditions by clusters with Wilcoxon test")) %>% 
      kable_styling( full_width = FALSE)
  )
  
}

# MONITORED GENES
#################

## @knitr heterogeneity_monitoredGenes

# Just remind the warning for genes names not in object
if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( monitoredGenesNotFound, collapse=" - ")));



## @knitr heterogeneity_monitoredGenes_heatmap

# DoHeatmap replaced by use of iHeatmapr after testing several options
#htmltools::tagList( ggplotly( DoHeatmap(object = sc10x, features = unlist(MONITORED_GENES))));


if( length( MONITORED_GENES) > 1){
  # Get the matrix of expression and associated clusters from Seurat object
  expMat = as.matrix( GetAssayData( sc10x));
  clusterID = Idents( sc10x);
  
  # Select monitored genes and reorder cells to group clusters together
  monitoredGenes = unlist( MONITORED_GENES);
  clusterOrdering = order( clusterID);
  
  expMat = expMat[monitoredGenes, clusterOrdering];
  clusterID = clusterID[ clusterOrdering];
  
  # Prepare rows and columns annotation bars (monitored group and cluster respectively)
  rowsAnnot = data.frame( Monitored = fct_rev(fct_inorder(rep( names( MONITORED_GENES), sapply( MONITORED_GENES, length)))));
  colsAnnot = data.frame( Cluster = clusterID);
  
  # Prepare unique rows and cols names for pheatmap (annotation rows) and match with rowAnnots and colAnnots row names
  originalRowNames = rownames( expMat);
  originalColNames = colnames( expMat);
  rownames( expMat) = make.unique( originalRowNames);
  colnames( expMat) = make.unique( originalColNames);
  rownames( rowsAnnot) = rownames( expMat);
  rownames( colsAnnot) = colnames( expMat);
  # Prepare colors of monitored genes groups for pheatmap (requires named vector matching data factor levels)
  monitoredColors = rainbow( nlevels( rowsAnnot[["Monitored"]]), s = 0.8);
  names( monitoredColors) = levels( rowsAnnot[["Monitored"]]);
  
  pheatmap( expMat,
            color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
            cluster_rows = FALSE,
            cluster_cols = FALSE,
            scale = "column",
            annotation_row = rowsAnnot,
            annotation_col = colsAnnot,
            labels_row = originalRowNames,
            annotation_colors = list( Monitored = monitoredColors,
                                      Cluster = clustersColor),
            show_colnames = FALSE);

}else{
  cat("<BR>Only one gene in list : not displaying heatmap")
}


## @knitr heterogeneity_monitoredGenes_expression_projection
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Plot expression values of monitored genes (TODO: message if list empty)
invisible( lapply( names( MONITORED_GENES), function(monitoredGroup)
{
  cat("#####", monitoredGroup, "\n");

  # Plots expression on projected cells (or error message if feature not found)
  invisible( lapply( MONITORED_GENES[[monitoredGroup]], function(featureName)
  {
    print(
      tryCatch( suppressWarnings( # Raise a warning # A warning is raised when all values are 0 (e.g. contamination genes after filtering)
                FeaturePlot( sc10x, features = featureName, reduction = ifelse(exists("useReduction"), useReduction, "umap"), order = TRUE)  +
                  theme( axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         legend.position = "none")),
                error = function(e){return(conditionMessage(e))})); # If gene not found in object (should have been removed earlier anyway)
  }));

  cat(" \n \n"); # Required for '.tabset'
}));




## @knitr heterogeneity_monitoredGenes_expression_violin

# Plot expression values of monitored genes as violinplot for each cluster (TODO: message if list empty)
invisible( lapply( names( MONITORED_GENES), function(monitoredGroup)
{
  cat("#####", monitoredGroup, "\n");

  # Violinplot for expression value of monitored genes by cluster (+ number of 'zero' and 'not zero' cells)
  invisible( lapply( MONITORED_GENES[[monitoredGroup]], violinFeatureByCluster, seuratObject = sc10x, clustersColor = clustersColor));

  cat(" \n \n"); # Required for '.tabset'
}));




# MODULES ANALYSIS
##################

## @knitr heterogeneity_modules

# Just remind the warning for genes names not in object, or modules that were transfered to individual monitoring of genes
if(any( is.na( matchModulesGenes)))
{
    warning( paste( "Following gene(s) from modules list could not be found in experimental data and will be ignored:",
                    paste( modulesGenesNotFound, collapse=" - ")));
}

if( exists( "modulesToTransfer") && any(modulesToTransfer))
{
    warning( paste0( "Following Module(s) contained very few genes (<",
                     MONITORED_GENES_SMALL_MODULE,
                     "). These genes were transfered to 'Monitored genes' to be analyzed individually: ",
                     paste( names(modulesToTransfer)[modulesToTransfer], collapse = " - ")));
}




## @knitr heterogeneity_modules_scoring

# Compute the score of the cells according to group of monitored genes
for( moduleName in names( MODULES_GENES))
{
  message(moduleName);
  if( length( MODULES_GENES[[moduleName]]) == 0)
  {
    warning( paste0( "List of genes in module '", moduleName, "' is empty, ignoring..."));
  } else
  {
    sc10x <- AddModuleScore( object = sc10x,
                             features = MODULES_GENES[ moduleName], # Must be a list
                             ctrl = MODULES_CONTROL_SIZE,           #length(MODULES_GENES[[ moduleName]]),
                             name = moduleName,
                             seed = SEED);
  }
}




## @knitr heterogeneity_modules_heatmap

# DoHeatmap replaced by use of iHeatmapr after testing several options
#htmltools::tagList( ggplotly( DoHeatmap(object = sc10x, features = topMarkersDF[["gene"]])));

# Get the matrix of module scores (not expression) for each cell and associated clusters from Seurat object
modulesScoreMat = t( as.matrix( sc10x[[ paste0(names( MODULES_GENES),1)]]));
clusterID = Idents( sc10x);

# Remove the extra numeric character added to modules names by Seurat
rownames( modulesScoreMat) = substr( rownames( modulesScoreMat), 1, nchar( rownames( modulesScoreMat))-1);

# Select genes in modules and reorder cells to group clusters together
clusterOrdering = order( clusterID);

modulesScoreMat = modulesScoreMat[, clusterOrdering];
clusterID = clusterID[ clusterOrdering];

# Prepare rows and columns annotation bars (module and cluster respectively)
#rowsAnnot = data.frame( Module = names( MODULES_GENES));
colsAnnot = data.frame( Cluster = clusterID);

# Prepare unique rows and cols names for pheatmap (annotation rows) and match with rowAnnots and colAnnots row names
originalRowNames = rownames( modulesScoreMat);
originalColNames = colnames( modulesScoreMat);
rownames( modulesScoreMat) = make.unique( originalRowNames);
colnames( modulesScoreMat) = make.unique( originalColNames);
#rownames( rowsAnnot) = rownames( modulesScoreMat);
rownames( colsAnnot) = colnames( modulesScoreMat);

# Plot the 'non-interactive' heatmap
pheatmap( modulesScoreMat,
          color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          scale = "column",
          annotation_col = colsAnnot,
          labels_row = originalRowNames,
          annotation_colors = list( Cluster = clustersColor),
          show_colnames = FALSE);


## @knitr heterogeneity_modules_expression_projection
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Plot scores of modules on dimensionality reduction figures
invisible( lapply( names(MODULES_GENES), function(featureName)
{
  print( FeaturePlot( sc10x, features = paste0(featureName, "1"), reduction = ifelse( exists("useReduction"), useReduction, "umap"), order = TRUE) +
           ggtitle( label = featureName) +
           theme( axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = "none"));
}));

cat(" \n \n"); # Required for '.tabset'




## @knitr heterogeneity_modules_expression_violin

# Violinplot for module values by cluster
invisible( lapply( paste0(names(MODULES_GENES), 1), violinFeatureByCluster, seuratObject = sc10x, clustersColor = clustersColor, yLabel = "Score", addStats = FALSE, trimTitle = 1));

cat(" \n \n"); # Required for '.tabset'




## @knitr heterogeneity_modules_expression_projection_pngFile
# Plot expression values of individual module genes as png files (not in report)
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Compute 'pixel' dimensions to match figures in report (based on chunk params)
defaultDim = c( 800, 800);  # Fallback values if chunk params not available
defaultDpi = 72;            # (e.g. executing directly from R)

chunkDim    = knitr::opts_current$get( "fig.dim");
chunkDpi    = knitr::opts_current$get( "dpi");

figDim = if( is.null( chunkDim) || is.null( chunkDpi) ) defaultDim else chunkDim * chunkDpi;
figDpi = if( is.null( chunkDpi) ) defaultDpi else chunkDpi;


invisible( lapply( names( MODULES_GENES), function(moduleName)
{
  message(moduleName); # Just for tracking progress in console

  # Create subfolder for current module to store all png files
  pathCurrentModule = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "modulesExpressionIndividual_", ifelse(exists("useReduction"), useReduction, "umap")), moduleName);
  dir.create(pathCurrentModule, showWarnings = FALSE, recursive = TRUE);

  # Plots expression on projected cells (or error message if feature not found)
  invisible( lapply( MODULES_GENES[[moduleName]], function(featureName)
  {
    # Create png file
    png( file.path( pathCurrentModule, paste0( featureName, ".png")), 
         width = figDim[1],
         height = figDim[2],
         res = figDpi);

    print(FeaturePlot( sc10x, features = featureName, reduction = ifelse(exists("useReduction"), useReduction, "umap"), order = TRUE)  +
      theme( axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none"))

    dev.off(); # Close file descriptor
  }));

}));




## @knitr heterogeneity_modules_expression_violin_pngFile
# Plot expression values of monitored genes as violinplot in a png file for each cluster (TODO: message if list empty)

# Compute 'pixel' dimensions to match figures in report (based on chunk params)
defaultDim = c( 800, 800);  # Fallback values if chunk params not available
defaultDpi = 72;            # (e.g. executing directly from R)

chunkDim    = knitr::opts_current$get( "fig.dim");
chunkDpi    = knitr::opts_current$get( "dpi");

figDim = if( is.null( chunkDim) || is.null( chunkDpi) ) defaultDim else chunkDim * chunkDpi;
figDpi = if( is.null( chunkDpi) ) defaultDpi else chunkDpi;


invisible( lapply( names( MODULES_GENES), function(moduleName)
{
  message(moduleName); # Just for tracking progress in console

  # Create subfolder for current module to store all png files
  pathCurrentModule = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "modulesExpressionIndividual_violin"), moduleName);
  dir.create(pathCurrentModule, showWarnings = FALSE, recursive = TRUE);

  # Violinplot for expression value of monitored genes by cluster (+ number of 'zero' and 'not zero' cells)
  invisible( lapply( MODULES_GENES[[moduleName]], function(featureName)
  {
    # Create png file
    png( file.path( pathCurrentModule, paste0( featureName, ".png")), 
         width = figDim[1],
         height = figDim[2],
         res = figDpi);

    violinFeatureByCluster(featureName, seuratObject = sc10x, clustersColor = clustersColor);

    dev.off(); # Close file descriptor
  }))

}));


