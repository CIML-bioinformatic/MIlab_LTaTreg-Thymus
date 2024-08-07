#' Harmony single cell integration
#'
#' Run Harmony algorithm with Seurat and SingleCellAnalysis pipelines.
#'
#' @param object Pipeline object. Must have PCA computed.
#' @param group.by.vars Which variable(s) to remove (character vector).
#' @param dims.use Which PCA dimensions to use for Harmony. By default, use all
#' @param theta Diversity clustering penalty parameter. Specify for each
#' variable in group.by.vars. Default theta=2. theta=0 does not encourage any
#'  diversity. Larger values of theta result in more diverse clusters.
#' @param lambda Ridge regression penalty parameter. Specify for each variable
#' in group.by.vars. Default lambda=1. Lambda must be strictly positive.
#' Smaller values result in more aggressive correction.
#' @param sigma Width of soft kmeans clusters. Default sigma=0.1. Sigma scales
#' the distance from a cell to cluster centroids. Larger values of sigma result
#'  in cells assigned to more clusters. Smaller values of sigma make soft
#'  kmeans cluster approach hard clustering.
#' @param nclust Number of clusters in model. nclust=1 equivalent to simple
#'  linear regression.
#' @param tau Protection against overclustering small datasets with large ones.
#'  tau is the expected number of cells per cluster.
#' @param block.size What proportion of cells to update during clustering.
#'  Between 0 to 1, default 0.05. Larger values may be faster but less accurate
#' @param max.iter.cluster Maximum number of rounds to run clustering at each
#' round of Harmony.
#' @param epsilon.cluster Convergence tolerance for clustering round of Harmony
#'  Set to -Inf to never stop early.
#' @param max.iter.harmony Maximum number of rounds to run Harmony. One round
#'  of Harmony involves one clustering and one correction step.
#' @param epsilon.harmony Convergence tolerance for Harmony. Set to -Inf to
#'  never stop early.
#' @param plot_convergence Whether to print the convergence plot of the
#' clustering objective function. TRUE to plot, FALSE to suppress. This can be
#'  useful for debugging.
#' @param verbose Whether to print progress messages. TRUE to print, FALSE to
#'  suppress.
#' @param reference_values (Advanced Usage) Defines reference dataset(s). Cells
#' that have batch variables values matching reference_values will not be moved
#' @param reduction.save Keyword to save Harmony reduction. Useful if you want
#' to try Harmony with multiple parameters and save them as e.g.
#' 'harmony_theta0', 'harmony_theta1', 'harmony_theta2'
#' @param assay.use (Seurat V3 only) Which assay to run PCA on if no PCA present?
#' @param ... other parameters
#'
#'
#' @rdname Run_Harmony
#' @export
Run_Harmony <- function(object, group.by.vars, ...) {
  UseMethod("Run_Harmony")
}



#' @rdname Run_Harmony
#' @param reduction Name of dimension reduction to use. Default is PCA.
#' @param project.dim Project dimension reduction loadings. Default TRUE.
#' @return Seurat (version 3) object. Harmony dimensions placed into
#' dimensional reduction object harmony. For downstream Seurat analyses,
#' use reduction='harmony'.
#' @export
Run_Harmony.Seurat <- function(
  object,
  group.by.vars,
  reduction = 'pca',
  dims.use = NULL,
  theta = NULL,
  lambda = NULL,
  sigma = 0.1,
  nclust = NULL,
  tau = 0,
  block.size = 0.05,
  max.iter.harmony = 10,
  max.iter.cluster = 20,
  epsilon.cluster = 1e-5,
  epsilon.harmony = 1e-4,
  plot_convergence = FALSE,
  verbose = TRUE,
  reference_values = NULL,
  reduction.save = "harmony",
  assay.use = NULL,
  project.dim = TRUE,
  ...
) {
  if (!requireNamespace('Seurat', quietly = TRUE)) {
    stop("Running Harmony on a Seurat object requires Seurat")
  }
  assay.use <- assay.use %||% Seurat::DefaultAssay(object)
  if (reduction == "pca" && !reduction %in% Seurat::Reductions(object = object)) {
    if (isTRUE(x = verbose)) {
      message("Harmony needs PCA. Trying to run PCA now.")
    }
    object <- tryCatch(
      expr = Seurat::RunPCA(
        object = object,
        assay = assay.use,
        verbose = verbose,
        reduction.name = reduction
      ),
      error = function(...) {
        stop("Harmony needs PCA. Tried to run PCA and failed.")
      }
    )
  }
  if (!reduction %in% Seurat::Reductions(object = object)) {
    stop("Requested dimension reduction is not present in the Seurat object")
  }
  embedding <- Seurat::Embeddings(object, reduction = reduction)
  if (is.null(dims.use)) {
    dims.use <- seq_len(ncol(embedding))
  }
  dims_avail <- seq_len(ncol(embedding))
  if (!all(dims.use %in% dims_avail)) {
    stop("trying to use more dimensions than computed. Rereun dimension reduction
         with more dimensions or run Harmony with fewer dimensions")
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }
  metavars_df <- Seurat::FetchData(
    object,
    group.by.vars,
    cells = Seurat::Cells(x = object[[reduction]])
  )
  
  harmonyEmbed <- HarmonyMatrix(
    embedding,
    metavars_df,
    group.by.vars,
    FALSE,
    0,
    theta,
    lambda,
    sigma,
    nclust,
    tau,
    block.size,
    max.iter.harmony,
    max.iter.cluster,
    epsilon.cluster,
    epsilon.harmony,
    plot_convergence,
    FALSE,
    verbose,
    reference_values
  )
  
  # This is the line to change to fix the issue 
  # Error in UseMethod(generic = "Key", object = object) :
  # no applicable method for 'Key' applied to an object of class "character"
  # ------------------
  # reduction.key <- Seurat::Key(reduction.save, quiet = TRUE)
  reduction.key = reduction.save
  # ------------------
  rownames(harmonyEmbed) <- rownames(embedding)
  colnames(harmonyEmbed) <- paste0(reduction.key, seq_len(ncol(harmonyEmbed)))
  
  object[[reduction.save]] <- Seurat::CreateDimReducObject(
    embeddings = harmonyEmbed,
    stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
    assay = Seurat::DefaultAssay(object = object[[reduction]]),
    key = reduction.key
  )
  if (project.dim) {
    object <- Seurat::ProjectDim(
      object,
      reduction = reduction.save,
      overwrite = TRUE,
      verbose = FALSE
    )
  }
  return(object)
}



#' @rdname Run_Harmony
#' @return SingleCellExperiment object. After running Run_Harmony, the corrected
#' cell embeddings can be accessed with reducedDim(object, "Harmony").
#' @export
Run_Harmony.SingleCellExperiment <- function(
  object,
  group.by.vars,
  dims.use = NULL,
  theta = NULL,
  lambda = NULL,
  sigma = 0.1,
  nclust = NULL,
  tau = 0,
  block.size = 0.05,
  max.iter.harmony = 10,
  max.iter.cluster = 20,
  epsilon.cluster = 1e-5,
  epsilon.harmony = 1e-4,
  plot_convergence = FALSE,
  verbose = TRUE,
  reference_values = NULL,
  reduction.save = "HARMONY",
  ...
) {
  
  ## Get PCA embeddings
  if (!"PCA" %in% SingleCellExperiment::reducedDimNames(object)) {
    stop("PCA must be computed before running Harmony.")
  }
  pca_embedding <- SingleCellExperiment::reducedDim(object, "PCA")
  if (is.null(dims.use)) {
    dims.use <- seq_len(ncol(pca_embedding))
  }
  
  if (is.null(dims.use)) {
    dims.use <- seq_len(ncol(pca_embedding))
  }
  dims_avail <- seq_len(ncol(pca_embedding))
  if (!all(dims.use %in% dims_avail)) {
    stop("trying to use more dimensions than computed with PCA. Rereun
            PCA with more dimensions or use fewer PCs")
  }
  
  metavars_df <- SingleCellExperiment::colData(object)
  if (!all(group.by.vars %in% colnames(metavars_df))) {
    stop('Trying to integrate over variables missing in colData')
  }
  
  harmonyEmbed <- HarmonyMatrix(
    pca_embedding,
    metavars_df,
    group.by.vars,
    FALSE,
    0,
    theta,
    lambda,
    sigma,
    nclust,
    tau,
    block.size,
    max.iter.harmony,
    max.iter.cluster,
    epsilon.cluster,
    epsilon.harmony,
    plot_convergence,
    FALSE,
    verbose,
    reference_values
  )
  
  rownames(harmonyEmbed) <- row.names(metavars_df)
  colnames(harmonyEmbed) <- paste0(reduction.save, "_", seq_len(ncol(harmonyEmbed)))
  SingleCellExperiment::reducedDim(object, reduction.save) <- harmonyEmbed
  
  return(object)
}