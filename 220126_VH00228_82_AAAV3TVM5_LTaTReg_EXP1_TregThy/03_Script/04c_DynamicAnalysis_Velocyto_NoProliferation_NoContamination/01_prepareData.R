# ################################################
# This script aims to read and filter data
# ################################################

# This script has been inspired by the R Velocyto tutorial page:
# http://pklab.med.harvard.edu/velocyto/notebooks/R/DG1.nb.html

## @knitr prepare_data

# READ THE DATA
# ...............

## @knitr prepare_data

# Load the Seurat object from heterogenity analysis
sc10x.seurat = readRDS( SEURAT_OBJECT_HETEROGENEITY_ANALYSIS)

# Get the UMAP and Best Archetype from the Seurat object
umap_embedding = as.data.frame( sc10x.seurat@reductions$umap@cell.embeddings)

# Convert the names of the cells (they are different in the Seurat Object and in the loom object)
# RUN1_SPL1_ACGGGTCAGATGATTG-1 becomes mm10:ACGGGTCAGATGATTGx
row.names( umap_embedding) = gsub( "RUN1_THY1_", "mm10:", row.names( umap_embedding))
row.names( umap_embedding) = gsub( "-1", "x", row.names( umap_embedding))

# Read the loom file from velocyto pre-processing
input_file = file.path( VELOCYTO_LOOM_FILE)
ldat <- read.loom.matrices( input_file)

cat("<BR>Analyzed data are located on:", input_file)

# Get expression matrix
emat <- ldat$spliced

# Filter some very low values
emat = emat[ , row.names( umap_embedding)]
emat <- emat[,colSums(emat)>=1e3]


# Pagoda2 pre-processing
# .......................

## @knitr build_data
cat("<BR>Building the required data for velocity analysis")

# Create the Pagoda2 object
rownames( emat) <- make.unique( rownames( emat))
cat("<BR>Number of genes in expression matrix:", length( rownames( emat)))

r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)

# Adjust the variances
r$adjustVariance(plot=T,do.par=T,gam.k=10)

# Run basic analysis steps to generate cell embedding and clustering, visualize:

# - Run the PCA
set.seed( 123)
cat("<HR>")
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)
cat("<HR>")

# - Cluster the cells
r$makeKnnGraph(k=40,type='PCA',center=T,distance='cosine');
r$getKnnClusters(method=multilevel.community, type='PCA', name='multilevel')
r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=FALSE)

# - Look at the distribution of cells along cluster
cluster_df = data.frame( table( r$clusters[[ "PCA"]]$multilevel))
names( cluster_df) = c( "Cluster ID", "Cell Number")
cluster_df %>% kbl( caption="Distribution of cells along clusters", rownames=FALSE) %>%kable_styling( full_width = FALSE)

# Compute differentially expressed genes
r$getDifferentialGenes(type='PCA', verbose=T, clusterType='multilevel')
de = data.frame()
for( cluster_name in names( r$diffgenes[[ "PCA"]]$multilevel)){
  current_df = head( r$diffgenes[[ "PCA"]]$multilevel[[ cluster_name]], 10)
  current_df$cluster = rep( paste0( "cluster", cluster_name), nrow( current_df))
  de = rbind( de, current_df)
}
cat("<HR>")
#r$plotGeneHeatmap(genes = rownames( de), groups = r$clusters$PCA[[1]])
cat("<HR>")
de %>% kbl( caption = "Differentially expressed genes by cluster") %>% kable_styling()
cat("<HR>")

# Plot embedding, labeling clusters
# par( mfrow = c(1,1))
# cat("<HR>")
# r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,main='cell clusters')
# cat("<HR>")