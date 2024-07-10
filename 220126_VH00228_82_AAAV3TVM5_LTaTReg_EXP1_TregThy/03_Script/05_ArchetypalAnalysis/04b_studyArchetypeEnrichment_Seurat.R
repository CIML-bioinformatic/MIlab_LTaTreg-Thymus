# ######################################################
# For each archetype, study the most differentially expressed genes
# using Seurat FindMarkers
# ######################################################


## @knitr study_archetype_enrichment_seurat

# Define the list that will contains the lists of DEG for each archetype
deg_result_list = list()

# For each archetype, look at the DEG
# ..........................................
Idents( sc10x.seurat) = "BestArchetype"
all_markers = FindAllMarkers( object  = sc10x.seurat,
                test.use        = FINDMARKERS_METHOD,
                only.pos        = FALSE,
                min.pct         = FINDMARKERS_MINPCT,
                logfc.threshold = FINDMARKERS_LOGFC_THR,
                return.thresh   = FINDMARKERS_PVAL_THR,
                random.seed     = SEED,
                verbose         = .VERBOSE);

names( all_markers) = gsub( "cluster", "archetype", names( all_markers))

for( i in 1:ncol( usage_df)){
  
  # Split the results into several dataframe according to archetype
  deg_result_list[[ i]] = all_markers[ all_markers$archetype == i,]
  
  # Define the file to export the results.
  deg_file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "archetype_usage_K", AUTOENCODER_LAYER_INTERNAL_NODES, "_N", AUTOENCODER_LAYER_IN_OUT_NODES, "_SeuratDEG", i,".csv"))
  
  # Write the result to file
  write.csv( deg_result_list[[ i]], file = deg_file)
}


# For each archetype, filter the DEG according to defined threshold and display the result
# .........................................................................................
positive_filtered_deg_result_list = list()
negative_filtered_deg_result_list = list()

for( i in 1:ncol( usage_df)){

  # Extract the significant positive markers from global results  
  positive_filtered_deg_result_list[[ i]] = deg_result_list[[ i]][ which( deg_result_list[[i]]$p_val_adj < DEG_ALPHA_THRESHOLD &
                                                                            deg_result_list[[i]]$avg_log2FC > 0), ]
  
  # Reorder the positive markers by increasing p-values and decreasing estimate
  positive_filtered_deg_result_list[[ i]] = positive_filtered_deg_result_list[[ i]][ 
                                                          order( 
                                                              positive_filtered_deg_result_list[[ i]]$p_val_adj, 
                                                              positive_filtered_deg_result_list[[i]]$avg_log2FC, 
                                                              decreasing = c( FALSE, TRUE) ),
                                                          ] 
  
  # Extract the significant negative markers from global results  
  negative_filtered_deg_result_list[[ i]] = deg_result_list[[ i]][ which( deg_result_list[[i]]$p_val_adj < DEG_ALPHA_THRESHOLD &
                                                                            deg_result_list[[i]]$avg_log2FC < 0), ]
  
  # Reorder the negative markers by increasing p-values and increasing estimate
  negative_filtered_deg_result_list[[ i]] = negative_filtered_deg_result_list[[ i]][ 
                                                          order( 
                                                            negative_filtered_deg_result_list[[ i]]$p_val_adj, 
                                                            negative_filtered_deg_result_list[[i]]$avg_log2FC, 
                                                            decreasing = c( FALSE, FALSE) ),
                                                          ]
  
  # Display the result for the current archetype in a specific HTML tab
  # ...............................................................................
  cat( "\n \n")
  cat( "\n#### Archetype", i, " {.tabset .tabset-fade}")
  cat( "\n \n")
  
  # Add a sub-tab with the table of positive markers
  # ...............................................................................
  cat( "\n \n")
  cat( "\n##### Positive markers (top 50)")
  cat( "\n \n")
  
  print( 
    positive_filtered_deg_result_list[[ i]][1:50, ] %>% 
                                  dplyr::mutate_if(is.numeric, funs(as.character(signif(., 3))))  %>% 
                                  kbl() %>%
                                  kable_styling()
  )
  
  # Add a sub-tab with the enrichment of the positive markers in the modules genes
  # ...............................................................................
  cat( "\n \n")
  cat( "\n##### Positive markers Module Enrichment")
  cat( "\n \n")
  
  # Parse the modules genes to test their enrichment in the positive markers
  positive_enrichment_df = data.frame()
  for(module_name in names( MODULES_GENES)){
    
    # Compute a hypergeometric statistics from the modules genes and the positive markers
    universe = deg_result_list[[ i]]$gene
    white_balls = intersect( MODULES_GENES[[ module_name]], universe)
    black_balls = setdiff( universe, white_balls)
    sample = positive_filtered_deg_result_list[[ i]]$gene
    inter = intersect( sample, white_balls)
    
    pval = phyper( q = length( inter) -1,
            m = length( white_balls),
            n = length( black_balls),
            k = length( sample),
            lower.tail = FALSE
           )
    
    # Store the enrichment result in a global table
    positive_enrichment_df = rbind( positive_enrichment_df, data.frame( Module = module_name,
                                                      All_genes = length( universe),
                                                      Module_genes = length( MODULES_GENES[[ module_name]]),
                                                      Kept_module_genes = length( white_balls),
                                                      DEG_genes = length( sample),
                                                      DEG_in_module = length( inter),
                                                      p_value = pval
                                                      ))
  }
  
  # Adjust the pvalues for multi-test
  positive_enrichment_df$Adjusted_pval = p.adjust( positive_enrichment_df$p_value, method = "BH")
  
  # Reorder the result by significance
  positive_enrichment_df = positive_enrichment_df[ order( positive_enrichment_df$Adjusted_pval, decreasing = FALSE), ]
  
  # Display the result in a static table
  print( 
    positive_enrichment_df %>%
      dplyr::mutate_if(is.numeric, funs(as.character(signif(., 3))))  %>% 
      kbl() %>%
      kable_styling()
  )
  
  # Add a sub-tab with the table of positive markers
  # ...............................................................................
  cat( "\n \n")
  cat( "\n##### Negative markers (top 50)")
  cat( "\n \n")
  
  print( 
    negative_filtered_deg_result_list[[ i]][1:50, ] %>% 
                                  dplyr::mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>%
                                  kbl() %>%
                                  kable_styling()
  )
  
  # Add a sub-tab with the enrichment of the negative markers in the modules genes
  # ...............................................................................
  cat( "\n \n")
  cat( "\n##### Negative markers Module Enrichment")
  cat( "\n \n")
  
  # Parse the modules genes to test their enrichment in the negative markers
  negative_enrichment_df = data.frame()
  for(module_name in names( MODULES_GENES)){
    
    # Compute a hypergeometric statistics from the modules genes and the positive markers
    universe = deg_result_list[[ i]]$gene
    white_balls = intersect( MODULES_GENES[[ module_name]], universe)
    black_balls = setdiff( universe, white_balls)
    sample = negative_filtered_deg_result_list[[ i]]$gene
    inter = intersect( sample, white_balls)
    
    pval = phyper( q = length( inter) -1,
                   m = length( white_balls),
                   n = length( black_balls),
                   k = length( sample),
                   lower.tail = FALSE
    )

    # Store the enrichment result in a global table    
    negative_enrichment_df = rbind( negative_enrichment_df, data.frame( Module = module_name,
                                                      All_genes = length( universe),
                                                      Module_genes = length( MODULES_GENES[[ module_name]]),
                                                      Kept_module_genes = length( white_balls),
                                                      DEG_genes = length( sample),
                                                      DEG_in_module = length( inter),
                                                      p_value = pval
    ))
  }
  
  # Adjust the pvalues for multi-test
  negative_enrichment_df$Adjusted_pval = p.adjust( negative_enrichment_df$p_value, method = "BH")
  
  # Reorder the result by significance
  negative_enrichment_df = negative_enrichment_df[ order( negative_enrichment_df$Adjusted_pval, decreasing = FALSE), ]
  
  # Display the result in a static table
  print( 
    negative_enrichment_df %>% 
      dplyr::mutate_if(is.numeric, funs(as.character(signif(., 3))))  %>% 
      kbl() %>%
      kable_styling()
  )
}







