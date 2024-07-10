# ####################################################
# This script aim to use the Velocyto R pipeline
# to analyse the RNA velocity data
# ####################################################

# This script has been inspired by the R Velocyto tutorial page:
# http://pklab.med.harvard.edu/velocyto/notebooks/R/DG1.nb.html


## @knitr estimate_velocity

# Get the raw data for spliced and unspliced RNA
emat <- ldat$spliced
nmat <- ldat$unspliced

# Restrict to cells that passed pagoda2 filter
emat <- emat[,rownames(r$counts)]
nmat <- nmat[,rownames(r$counts)]

# Take cluster labels from Pagoda2 pre-processing
cluster.label <- r$clusters$PCA[[ "multilevel"]]
cell.colors <- sccore:::fac2col(cluster.label)

# Take t-SNE embedding from Pagoda2 pre-processing
# emb <- r$embeddings$PCA$tSNE
seurat_emb = as.matrix( umap_embedding[ , c("UMAP_1", "UMAP_2")])

# In addition to clustering and the t-SNE embedding, from the pagoda2 processing 
# we will also take a cell-cell distance, which will be better than the default 
# whole-transcriptome correlation distance that velocyto.R would normally use.
cell.dist <- as.dist(1-armaCor(t(r$reductions$PCA)))

# Filter genes based on the minimum average expresion magnitude 
# (in at least one of the clusters), output total number of resulting valid genes
emat <- filter.genes.by.cluster.expression(emat,cluster.label,min.max.cluster.average = 0.5)
nmat <- filter.genes.by.cluster.expression(nmat,cluster.label,min.max.cluster.average = 0.05)
cat("<BR>Number of validated genes for analysis =",length(intersect(rownames(emat),rownames(nmat))), "<BR>")
write.csv( data.frame( gene_name = sort( intersect(rownames(emat),rownames(nmat)))),
           file = file.path( PATH_ANALYSIS_OUTPUT, "velocyto_validated_genes.csv"))

# Estimate RNA velocity (using gene-relative model with k=20 cell kNN pooling and 
# using top/bottom 2% quantiles for gamma fit)
fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat,
                                            nmat,
                                            deltaT=1,
                                            kCells=20,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile)

# Visualize velocity on the t-SNE embedding, using velocity vector fields
cat("<HR>")
par( mfrow = c(1,1))
show.velocity.on.embedding.cor( seurat_emb,
                               rvel.cd,
                               n=100,
                               scale='sqrt',
                               cell.colors=ac( cell.colors, alpha=0.5),
                               cex=0.8,
                               arrow.scale=1.2,
                               show.grid.flow=FALSE,
                               # min.grid.cell.mass=0.5,
                               # grid.n=40,
                               arrow.lwd=1,
                               do.par=F,
                               cell.border.alpha = 0.1)

par( mfrow = c(1,1))
show.velocity.on.embedding.cor( seurat_emb,
                                rvel.cd,
                                n=100,
                                scale='sqrt',
                                cell.colors=ac( cell.colors, alpha=0.5),
                                cex=0.8,
                                arrow.scale=1.2,
                                show.grid.flow=TRUE,
                                min.grid.cell.mass=0.5,
                                grid.n=40,
                                arrow.lwd=1,
                                do.par=F,
                                cell.border.alpha = 0.1)

# Visualize a fit for maker genes (we reuse rvel.cd to save on calcualtions here):
cat("<HR>")
par( mfrow = c(2,2))
excluded_genes = vector()
for( gene_name in MONITORED_GENES){
  cat("<H4>RNA for gene", gene_name, "</H4>")
  if( gene_name %in% rownames( emat) && gene_name %in% rownames( nmat)){
    
    gene.relative.velocity.estimates(emat, nmat, deltaT=1, kCells = 20, kGenes=1,
                                     fit.quantile=fit.quantile, cell.emb=seurat_emb, cell.colors=cell.colors,
                                     cell.dist=cell.dist, old.fit=rvel.cd, do.par=FALSE,
                                     show.gene= gene_name)
  }else{
    cat("<BR>The gene", gene_name, "is not present in the matrices : can't display graphs.")
  }
}
par( mfrow = c(1,1))