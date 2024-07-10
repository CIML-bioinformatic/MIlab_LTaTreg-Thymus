# ########################################################################################
# This script aims to perform KEGG enrichment analysis
# ########################################################################################

## @knitr kegg_functional_enrichment

# Get the biomart ensembl reference
httr::set_config(httr::config(ssl_verifypeer = FALSE))
mart <- biomaRt::useDataset("mmusculus_gene_ensembl", biomaRt::useMart("ensembl"))

for( clusterid in unique( sample_cluster_markergenes_df$cluster)){
  
  current_sample_cluster_markergenes_df = sample_cluster_markergenes_df[ which( sample_cluster_markergenes_df$cluster == clusterid),]
  all_genes = current_sample_cluster_markergenes_df$gene
  
  # Convert genes names to EntrezID and order the genes by decreasing LogFC
  
  # Get the entrezID of the genes
  all_genes_entrezid <- biomaRt::getBM( filters="external_gene_name",
                                        attributes=c("entrezgene_id", "external_gene_name"),
                                        values= all_genes,
                                        mart=mart)
  
  # Look at NA in the entrezid and remove them
  indNA = which(is.na(all_genes_entrezid$entrezgene_id))
  if( length( indNA) > 0 ){
    all_genes_entrezid_noNA <- all_genes_entrezid[-indNA,]
  }else{
    all_genes_entrezid_noNA = all_genes_entrezid
  }
  
  # Look at the duplicates in entrezID and remove them
  indnodup = which(duplicated( all_genes_entrezid_noNA$ entrezgene_id) == F)
  all_genes_entrezid_noNA_nodup <- all_genes_entrezid_noNA[ indnodup,]
  
  # Keep only the entrezid that are no NA and not replicated in the DEseq result
  all_genes_entrezid_noNA_nodup_merge_deg = merge( all_genes_entrezid_noNA_nodup, current_sample_cluster_markergenes_df, by.x = "external_gene_name", by.y = "gene", all.y = FALSE)
  
  
  # Sort fold changes in decreasing order
  lFC_entrezid = all_genes_entrezid_noNA_nodup_merge_deg$avg_log2FC
  names( lFC_entrezid) = all_genes_entrezid_noNA_nodup_merge_deg$entrezgene_id
  lFC_entrezid <- sort( lFC_entrezid, decreasing = TRUE)
  
  # Enrichment analysis using GSEA on KEGG Pathways
  # ______________________________________________________________
  
  # Execute GSEA KEGG analysis
  gseaKEGG <- gseKEGG(geneList = lFC_entrezid,
                      organism = "mmu",
                      # nPerm = 1000, # default number permutations
                      minGSSize = 5, # minimum gene set size
                      pvalueCutoff = 0.2, # padj cutoff value
                      verbose = FALSE)
  # Extract the GSEA results
  gseaKEGG_results <- gseaKEGG@result
  
  if( nrow( gseaKEGG_results) > 0){
    gseaKEGG_results = merge( gseaKEGG_results, all_genes_entrezid_noNA_nodup_merge_deg, by.x="ID", by.y = "entrezgene_id", all.y = FALSE)
    print( DT::datatable( gseaKEGG_results, caption = "Enriched KEGG pathways"))
  }else{
    cat("<BR><b>No KEGG enrichment for the list of DEG of cluster", clusterid, "</b>")
  }
}
