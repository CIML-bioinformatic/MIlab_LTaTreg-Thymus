# ##################################################
# Retrieve the archetype usage and study it in the
# Seurat embeddings
# ##################################################


## @knitr study_archetype_usage

# Get back from python the archetype definition
usage_df = py$usage

# For each archetype, plot the distribution of score on UMAP
for( archetype_index in 1:AUTOENCODER_LAYER_INTERNAL_NODES){
  archetype_name = paste0( "Archetype_", archetype_index)
  sc10x.seurat = AddMetaData( sc10x.seurat, metadata = usage_df[ Cells( sc10x.seurat), archetype_index], col.name = archetype_name)
  print( FeaturePlot( sc10x.seurat, features = archetype_name))
}

umap_df = data.frame( UMAP_1 = sc10x.seurat@reductions$umap@cell.embeddings[ , "UMAP_1"],
                      UMAP_2 = sc10x.seurat@reductions$umap@cell.embeddings[ , "UMAP_2"],
                      best_archetype = as.factor( apply (usage_df[ row.names( sc10x.seurat@reductions$umap@cell.embeddings), ], 1, which.max)),
                      alpha = apply (usage_df[ row.names( sc10x.seurat@reductions$umap@cell.embeddings), ], 1, function( usage){
                        return( (max(usage) - (1/AUTOENCODER_LAYER_INTERNAL_NODES))/ (1- (1/AUTOENCODER_LAYER_INTERNAL_NODES)))
                      })
                    )

# Add the best archetype as metadata in the seurat object
sc10x.seurat = AddMetaData( sc10x.seurat, metadata = umap_df[ Cells( sc10x.seurat), "best_archetype"], col.name = "BestArchetype")

# Plot the archetype on UMAP
Idents( sc10x.seurat) = "BestArchetype"
DimPlot( sc10x.seurat, label = TRUE, label.size = 8) + ggtitle( "Best Archetype")

DimPlot( sc10x.seurat, label = TRUE, label.size = 8, split.by = "BestArchetype", ncol = 3) + ggtitle( "Best Archetype")

# Plot the seurat clusters on UMAP
Idents( sc10x.seurat) = "seurat_clusters"
DimPlot( sc10x.seurat, label = TRUE, label.size = 8) + ggtitle( "Seurat Clusters")

# Compare the UMAP clusters with the archetypes with a confusion table
confusion_df = as.data.frame.matrix( table( c(Idents( sc10x.seurat)), umap_df[ names( Idents( sc10x.seurat)), "best_archetype" ]))
names( confusion_df) = paste( "Archetype", names( confusion_df))
row.names( confusion_df) = paste( "Cluster", row.names( confusion_df))

confusion_df %>% kbl() %>% kable_styling( full_width = FALSE)

