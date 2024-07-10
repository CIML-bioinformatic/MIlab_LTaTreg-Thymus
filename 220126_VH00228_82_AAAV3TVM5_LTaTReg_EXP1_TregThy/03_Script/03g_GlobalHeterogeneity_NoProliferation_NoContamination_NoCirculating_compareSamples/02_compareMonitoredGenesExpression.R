# ##################################################
# This script compare the expression of selected
# genes between the two conditions
# ##################################################

# MONITORED GENES
#################

## @knitr heterogeneity_monitoredGenes

## @knitr heterogeneity_monitoredGenes_expression_projection
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Plot expression values of monitored genes (TODO: message if list empty)
invisible( lapply( names( MONITORED_GENES), function(monitoredGroup)
{
  cat("#####", monitoredGroup, " {.tabset .tabset-fade}\n");
  
  for( current_sample in c( SAMPLE_WT, SAMPLE_KO)){
  
    cat("######", current_sample, "\n");
    
    # Plots expression on projected cells (or error message if feature not found)
    invisible( lapply( MONITORED_GENES[[monitoredGroup]], function(featureName)
    {
      selected_cells = Cells( sc10x.seurat)[ which( sc10x.seurat@meta.data$HTO_classification == current_sample)]
      print(
        tryCatch( suppressWarnings( # Raise a warning # A warning is raised when all values are 0 (e.g. contamination genes after filtering)
          FeaturePlot( sc10x.seurat, cells = selected_cells, features = featureName, reduction = ifelse(exists("useReduction"), useReduction, "umap"), order = TRUE)  +
            theme( axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   legend.position = "none")),
          error = function(e){return(conditionMessage(e))})); # If gene not found in object (should have been removed earlier anyway)
    }));
  
    cat(" \n \n"); # Required for '.tabset'
  }
  cat(" \n \n"); # Required for '.tabset'
}));


## @knitr heterogeneity_monitoredGenes_density_projection
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)
invisible( lapply( names( MONITORED_GENES), function(monitoredGroup)
{
  cat("#####", monitoredGroup, " {.tabset .tabset-fade}\n");
  
  for( current_sample in c( SAMPLE_WT, SAMPLE_KO)){
    
    cat("######", current_sample, "\n");
    
    # Plots expression on projected cells (or error message if feature not found)
    invisible( lapply( MONITORED_GENES[[monitoredGroup]], function(featureName)
    {
      selected_cells = Cells( sc10x.seurat)[ which( sc10x.seurat@meta.data$HTO_classification == current_sample)]
      print(
        tryCatch( suppressWarnings( # Raise a warning # A warning is raised when all values are 0 (e.g. contamination genes after filtering)
          plot_density( subset( sc10x.seurat, cells = selected_cells), features = featureName, reduction = ifelse(exists("useReduction"), useReduction, "umap"))  +
            theme( axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   legend.position = "none")),
          error = function(e){return(conditionMessage(e))})); # If gene not found in object (should have been removed earlier anyway)
    }));
    
    cat(" \n \n"); # Required for '.tabset'
  }
  cat(" \n \n"); # Required for '.tabset'
}));


## @knitr heterogeneity_monitoredGenes_expression_violin

# Plot expression values of monitored genes as violinplot for each cluster (TODO: message if list empty)
invisible( lapply( names( MONITORED_GENES), function(monitoredGroup)
{
  cat("#####", monitoredGroup, " {.tabset .tabset-fade}\n");
  
  for( current_sample in c( SAMPLE_WT, SAMPLE_KO)){
    
    cat("######", current_sample, "\n");
    
    selected_cells = Cells( sc10x.seurat)[ which( sc10x.seurat@meta.data$HTO_classification == current_sample)]
    
    # Violinplot for expression value of monitored genes by cluster (+ number of 'zero' and 'not zero' cells)
    invisible( lapply( MONITORED_GENES[[monitoredGroup]], violinFeatureByCluster, seuratObject = subset( sc10x.seurat, cells = selected_cells), clustersColor = clustersColor));
    
    cat(" \n \n"); # Required for '.tabset'
  }
  cat(" \n \n"); # Required for '.tabset'
}))


