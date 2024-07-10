# ########################################################
# This script wraps up analysis (clean, save, render, ...)
# ########################################################




## @knitr save_final_seurat_object


# Save a binary file of final Seurat object only (as RDS)
seuratObjectPath = file.path( PATH_ANALYSIS_OUTPUT, paste0( "archetype_K", AUTOENCODER_LAYER_INTERNAL_NODES, "_N", AUTOENCODER_LAYER_IN_OUT_NODES, "_seuratObject_final.RDS"))
saveRDS( object = sc10x.seurat, file = seuratObjectPath)



