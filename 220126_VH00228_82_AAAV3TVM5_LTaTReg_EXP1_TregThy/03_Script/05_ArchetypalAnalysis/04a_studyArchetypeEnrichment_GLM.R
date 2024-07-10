# ######################################################
# Retrieve the archetype DEG and study their enrichment
# in gene modules
# ######################################################


## @knitr study_archetype_enrichment_glm

# ====================================================================
# For each archetype, study the most differentially expressed genes
# ====================================================================

# Define the list that will contains the lists of DEG for each archetype
deg_result_list = list()

# Define the cells with low expression
lib_size = apply( sc10x_counts_df, 1, sum)
min_counts_per_cell=1/20 * ncol( sc10x_counts_df)
cell_to_keep_index = which( lib_size >= min_counts_per_cell)

# Limit the counts to the selected cells
filtered_count = sc10x_counts_df[ cell_to_keep_index, ]
filtered_usage_df = usage_df[ cell_to_keep_index, ]
filtered_log_lib_size =  log( lib_size[ cell_to_keep_index])

min_cell_per_genes = max( 100, nrow( filtered_count)/500)

# Define the formula for the GLM regression
formula = "gene ~ usage + log_lib_size"



# For each archetype, look at the DEG
# ..........................................
for( i in 1:ncol( usage_df)){
  
  # Define the file to export the results. If the file exists, skip the analysis
  deg_file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "archetype_usage_K", AUTOENCODER_LAYER_INTERNAL_NODES, "_N", AUTOENCODER_LAYER_IN_OUT_NODES, "_DEG", i,".csv"))
  if( file.exists( deg_file)){
    deg_result_list[[ i]] = read.csv( deg_file, row.names = 1)
    next()
  }
  
  # Define the dataframe that will contain the DEG for the current archetype
  result = data.frame()
  excluded_genes = 0
  
  # For each genes, look at its expression relted to the archetype usage
  result = do.call( rbind, sapply( seq(1, ncol( filtered_count), 1), function( j){
  # result = do.call( rbind, lapply( seq(1, 100, 1), function( j){
    
    # cat("\ni/j=", i,"/",j)
    
    # Check if the current gene has 0 variation across cells
    if( sum( filtered_count[ ,j]>0) < min_cell_per_genes){
      excluded_genes = excluded_genes + 1
      return( data.frame( gene = names( filtered_count)[ j],
                          estimate = NA,
                          stderror = NA,
                          zvalue = NA,
                          pval = NA)
      )    
    }else{
      # Build the dataframe for regression
      data = data.frame( usage = filtered_usage_df[ , i],
                         gene = filtered_count[ ,j],
                         log_lib_size = filtered_log_lib_size
      )
      
      # Fit a Binomial GLM on data
      model = tryCatch( {
        MASS::glm.nb( formula = formula, data = data)
      },
      error=function(cond) {
        model = glm( formula, data = data, family=poisson)
        # cat("\n--POISSON-2")
        return( model)
      })
      
      # If the negative binomial GLM did not converge, use poisson GLM
      if ( !is.null( model$th.warn) && model$th.warn == "iteration limit reached") {
        model = glm( formula, data = data, family=poisson)
        # cat("\n--POISSON-1")
      }
      
      model_coef = coef( summary( model))
      
      # Add the result to result dataframe
      return( data.frame( gene = names( filtered_count)[ j],
                          estimate = model_coef[,1][ "usage"],
                          stderror = model_coef[,2][ "usage"],
                          zvalue = model_coef[,3][ "usage"],
                          pval = model_coef[,4][ "usage"])
      )    
    }
  }))
  
  # Remove all NA result from dataframe
  result = result[ which( !is.na( result$estimate)), ]
  
  # Adjust the p-values for multi-testing
  result$adjusted_pval = p.adjust( result$pval, method = "BH")
  
  # Write the result to file
  write.csv( result, file = deg_file)
  
  # Save the result to a global list
  deg_result_list[[ i]] = result
}


# For each archetype, filter the DEG according to defined threshold and display the result
# .........................................................................................
positive_filtered_deg_result_list = list()
negative_filtered_deg_result_list = list()

for( i in 1:ncol( usage_df)){

  # Extract the significant positive markers from global results  
  positive_filtered_deg_result_list[[ i]] = deg_result_list[[ i]][ which( deg_result_list[[i]]$adjusted_pval < DEG_ALPHA_THRESHOLD &
                                                                            deg_result_list[[i]]$estimate > 0), ]
  
  # Reorder the positive markers by increasing p-values and decreasing estimate
  positive_filtered_deg_result_list[[ i]] = positive_filtered_deg_result_list[[ i]][ 
                                                          order( 
                                                              positive_filtered_deg_result_list[[ i]]$adjusted_pval, 
                                                              positive_filtered_deg_result_list[[i]]$estimate, 
                                                              decreasing = c( FALSE, TRUE) ),
                                                          ] 
  
  # Change rownames to simple numbers
  row.names( positive_filtered_deg_result_list[[ i]]) = seq(1, nrow( positive_filtered_deg_result_list[[ i]]), 1)

  # Extract the significant negative markers from global results  
  negative_filtered_deg_result_list[[ i]] = deg_result_list[[ i]][ which( deg_result_list[[i]]$adjusted_pval < DEG_ALPHA_THRESHOLD &
                                                                            deg_result_list[[i]]$estimate < 0), ]
  
  # Reorder the negative markers by increasing p-values and increasing estimate
  negative_filtered_deg_result_list[[ i]] = negative_filtered_deg_result_list[[ i]][ 
                                                          order( 
                                                            negative_filtered_deg_result_list[[ i]]$adjusted_pval, 
                                                            negative_filtered_deg_result_list[[i]]$estimate, 
                                                            decreasing = c( FALSE, FALSE) ),
                                                          ]
  
  # Change rownames to simple numbers
  row.names( negative_filtered_deg_result_list[[ i]]) = seq(1, nrow( negative_filtered_deg_result_list[[ i]]), 1)
  
  
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
                                                      Deg_in_module = length( inter),
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
                                                      Deg_in_module = length( inter),
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







