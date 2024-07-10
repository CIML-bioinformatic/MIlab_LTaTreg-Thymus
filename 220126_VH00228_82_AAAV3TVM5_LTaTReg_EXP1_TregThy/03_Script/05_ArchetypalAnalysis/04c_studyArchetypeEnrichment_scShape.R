# ######################################################
# For each archetype, study the most differentially expressed genes
# using Seurat FindMarkers
# ######################################################


## @knitr study_archetype_enrichment_scshape


# We first filter the genes to keep only genes expressed in at least 10% of cells:
scData_filt <- filter_counts( as.matrix( sc10x.seurat@assays$RNA@counts), perc.zero = 0.1)

#sum( colnames( scData_filt) == Cells( sc10x.seurat))
lib_size = apply( scData_filt, 2, sum)

# Perform Kolmogorov-Smirnov test to select genes belonging to the family of ZINB distributions.
scData_KS <- ks_test( counts = as.matrix( sc10x.seurat@assays$RNA@counts), cexpr = sc10x.seurat@meta.data$BestArchetype, lib.size = lib_size, BPPARAM=SnowParam(workers=2,type="SOCK"))

# Select genes significant from the KS test.
# By default the 'ks_sig' function performs Benjamini-Hochberg correction for multiple   hypothese testing
# and selects genes significant at p-value of 0.01
scData_KS_sig <- ks_sig( scData_KS)

# Subset UMI counts corresponding to the genes significant from the KS test
scData.sig.genes <- scData_filt[rownames(scData_filt) %in% names(scData_KS_sig$genes),]

# Fit the 4 distributions P,NB,ZIP,ZINB for genes that belong to the ZINB family of distributions by fitting GLM with log of the library sizes as an offset 
# and cell types as a covariate in the GLM.
scData_models <- fit_models( counts = scData.sig.genes, cexpr = sc10x.seurat@meta.data$BestArchetype, lib.size = lib_size, BPPARAM = SnowParam(workers=2,type="SOCK"))

# Once the 4 distributions are fitted, we next calculate the BIC value for each model and select the model with the least BIC value.
scData_bicvals <- model_bic( scData_models)

# Select model with least bic value
scData_least.bic <- lbic_model( scData_bicvals, scData.sig.genes)

# To ensure the fit of the models selected based on the least BIC value, additionally we perform LRT to test for model adequacy and presence of zero-inflation.
scData_gof <- gof_model(scData_least.bic, cexpr=sc10x.seurat@meta.data$BestArchetype, lib.size=lib_size, BPPARAM=SnowParam(workers=2,type="SOCK"))

# Finally based on the results of the model adequacy tests, we can identify the distribution of best fit for each gene.
scData_fit <- select_model(scData_gof)

# Once the distribution of best fit is identified for genes of interest, it is also possible to extract parameters of interest for the models.
scData_params <- model_param(scData_models, scData_fit, model=NULL)




