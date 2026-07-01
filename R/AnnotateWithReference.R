#' Reference-based cell-type annotation (CellTypist / scANVI / scmap)
#'
#' Projects cell-type labels onto a query Seurat object using one of three
#' backends. All three write a predicted label column and (where the method
#' provides one) a confidence-score column into \code{obj@meta.data},
#' returning the modified Seurat object.
#'
#' \describe{
#'   \item{\code{method = "celltypist"} (default)}{Runs the Python
#'     \href{https://www.celltypist.org/}{CellTypist} logistic-regression
#'     classifier against one of its pre-trained models. Best default for
#'     broad tissue coverage — no user-supplied reference required, models
#'     for PBMC, tumor, lung, gut, kidney, developmental, and many others
#'     ship in the model zoo.
#'     Requires \code{reticulate} and a Python environment with
#'     \code{celltypist} installed (\code{pip install celltypist}).
#'     Uses the Seurat -> AnnData bridge from \code{ToAnnData()}.}
#'   \item{\code{method = "scanvi"}}{Runs
#'     \href{https://scvi-tools.org/}{scvi-tools} scANVI, a
#'     semi-supervised VAE that trains on a labeled reference and
#'     transfers labels to the query. Best when reference and query come
#'     from different technologies / protocols and you already have a
#'     high-quality reference.
#'     Requires \code{reticulate} + \code{scvi-tools} in Python, and a
#'     \code{reference} Seurat / AnnData object with cell-type labels.}
#'   \item{\code{method = "scmap"}}{Runs the Bioconductor package
#'     \code{scmap} (nearest-cluster or nearest-cell projection).
#'     R-native, no Python required. Requires a labeled \code{reference}
#'     (Seurat or SingleCellExperiment).}
#' }
#'
#' @param obj A Seurat object (query).
#' @param method Backend to use; see Details. Default \code{"celltypist"}.
#' @param assay Assay on \code{obj} to project from. Default DefaultAssay.
#' @param new_col Metadata column to write the predicted label to. Default
#'   \code{"predicted_cell_type"}.
#' @param min_score Optional confidence threshold. Cells whose backend
#'   score is below this get relabelled to \code{unassigned_label}.
#'   \code{NULL} (default) keeps every call.
#' @param unassigned_label Label applied to low-confidence cells when
#'   \code{min_score} is set. Default \code{"Unknown"}.
#'
#' @section CellTypist parameters:
#' \describe{
#'   \item{\code{model}}{Model name or file path. Default
#'     \code{"Immune_All_Low.pkl"} — good for PBMC-style data. Any model
#'     from \code{celltypist.models.get_all_models()} is valid. If not
#'     already downloaded, CellTypist fetches it on first use.}
#'   \item{\code{majority_voting}}{If TRUE (default), CellTypist first
#'     over-clusters the query and takes a majority vote within each
#'     cluster. Smoothes noisy per-cell calls.}
#'   \item{\code{python}}{Optional path to the Python executable used by
#'     \code{reticulate}. If \code{NULL}, uses the current
#'     \code{reticulate::use_python()} setting.}
#' }
#'
#' @section scANVI parameters:
#' \describe{
#'   \item{\code{reference}}{A Seurat object (or path to an \code{.h5ad}
#'     file) containing labeled reference cells.}
#'   \item{\code{ref_label_col}}{Metadata column on \code{reference}
#'     holding cell-type labels. Cells with \code{NA} labels are treated
#'     as unlabeled during training (the semi-supervised setting).}
#'   \item{\code{batch_col}}{Metadata column identifying batches on the
#'     \emph{combined} reference + query. Optional — pass \code{NULL} to
#'     assume a single batch.}
#'   \item{\code{n_latent}}{scVI latent dimensions. Default 30.}
#'   \item{\code{max_epochs_scvi}}{Max epochs for the underlying scVI
#'     model. Default 400.}
#'   \item{\code{max_epochs_scanvi}}{Max epochs for scANVI. Default 20.}
#'   \item{\code{python}}{As for CellTypist.}
#' }
#'
#' @section scmap parameters:
#' \describe{
#'   \item{\code{reference}}{A Seurat object or \code{SingleCellExperiment}
#'     with cell-type labels.}
#'   \item{\code{ref_label_col}}{Metadata column with the labels.}
#'   \item{\code{scmap_method}}{\code{"cluster"} (default; scmap-cluster)
#'     or \code{"cell"} (scmap-cell, product-quantization nearest neighbor).}
#'   \item{\code{threshold}}{Similarity threshold; scmap returns
#'     \code{"unassigned"} below this. Default 0.5 (scmap-cluster).}
#'   \item{\code{n_features}}{Number of features to select for the
#'     projection. Default 500.}
#' }
#'
#' @param ... Backend-specific parameters — see the sections above.
#' @return The Seurat object with a new \code{new_col} column (and, where
#'   available, a companion \code{<new_col>_score} numeric column).
#' @examples
#' \dontrun{
#' # 1. CellTypist (default) — no user-supplied reference
#' obj <- AnnotateWithReference(obj,
#'                              method          = "celltypist",
#'                              model           = "Immune_All_Low.pkl",
#'                              majority_voting = TRUE,
#'                              min_score       = 0.5)
#'
#' # 2. scANVI — from a labeled reference Seurat object
#' obj <- AnnotateWithReference(obj,
#'                              method        = "scanvi",
#'                              reference     = pbmc_ref_seurat,
#'                              ref_label_col = "cell_type",
#'                              batch_col     = "orig.ident")
#'
#' # 3. scmap-cluster — R-native
#' obj <- AnnotateWithReference(obj,
#'                              method        = "scmap",
#'                              reference     = pbmc_ref_seurat,
#'                              ref_label_col = "cell_type",
#'                              scmap_method  = "cluster")
#' }
#' @importFrom Seurat DefaultAssay
#' @export
AnnotateWithReference <- function(obj,
                                  method            = c("celltypist",
                                                        "scanvi",
                                                        "scmap"),
                                  assay             = NULL,
                                  new_col           = "predicted_cell_type",
                                  min_score         = NULL,
                                  unassigned_label  = "Unknown",
                                  ...) {

  method <- match.arg(method)
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay
  args <- list(...)

  message(sprintf("--- AnnotateWithReference (method = '%s') ---", method))
  labeled <- switch(
    method,
    celltypist = .annotate_celltypist(obj, a, new_col, args),
    scanvi     = .annotate_scanvi   (obj, a, new_col, args),
    scmap      = .annotate_scmap    (obj, a, new_col, args)
  )

  # Optional confidence threshold
  score_col <- paste0(new_col, "_score")
  if (!is.null(min_score) && score_col %in% colnames(labeled@meta.data)) {
    low <- labeled@meta.data[[score_col]] < min_score
    low[is.na(low)] <- TRUE
    n_low <- sum(low)
    if (n_low > 0) {
      labeled@meta.data[[new_col]][low] <- unassigned_label
      message(sprintf("  Applied min_score = %.2f: %d cell(s) -> '%s'.",
                      min_score, n_low, unassigned_label))
    }
  } else if (!is.null(min_score)) {
    warning("`min_score` set but backend '", method,
            "' did not return a score column; threshold ignored.",
            call. = FALSE)
  }

  message(sprintf("  Wrote '%s' with %d unique label(s).", new_col,
                  length(unique(labeled@meta.data[[new_col]]))))
  labeled
}


# ============================================================================
# CellTypist backend
# ============================================================================
#' @keywords internal
#' @noRd
.annotate_celltypist <- function(obj, assay, new_col, args) {

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("'reticulate' is required for method = 'celltypist'.")
  }
  if (!is.null(args$python)) reticulate::use_python(args$python,
                                                    required = TRUE)
  if (!reticulate::py_module_available("celltypist")) {
    stop("Python module 'celltypist' not available. Install with ",
         "`reticulate::py_install(\"celltypist\")` or ",
         "`pip install celltypist` in the active environment.")
  }

  model            <- args$model            %||% "Immune_All_Low.pkl"
  majority_voting  <- args$majority_voting  %||% TRUE

  # Convert query Seurat -> AnnData via ToAnnData if available, else via
  # zellkonverter / a direct reticulate write.
  h5ad_path <- tempfile(pattern = "sctools_query_", fileext = ".h5ad")
  on.exit(unlink(h5ad_path), add = TRUE)
  if (exists("ToAnnData", mode = "function")) {
    ToAnnData(obj, file = h5ad_path, assay = assay, layer = "counts")
  } else {
    stop("ToAnnData() not available; cannot convert Seurat to AnnData ",
         "for CellTypist. Ensure SeuratAnnDataInterop.R is loaded.")
  }

  celltypist <- reticulate::import("celltypist", delay_load = FALSE)
  anndata    <- reticulate::import("anndata",    delay_load = FALSE)
  message(sprintf("  Reading query into AnnData and running CellTypist ",
                  "(model = '%s', majority_voting = %s)...",
                  model, majority_voting))
  adata <- anndata$read_h5ad(h5ad_path)

  # CellTypist wants log-normalized counts on the X matrix. If ToAnnData
  # wrote raw counts, normalize here (same recipe scanpy uses).
  sc <- tryCatch(reticulate::import("scanpy", delay_load = FALSE),
                 error = function(e) NULL)
  if (!is.null(sc)) {
    sc$pp$normalize_total(adata, target_sum = 1e4)
    sc$pp$log1p(adata)
  }

  # CellTypist auto-downloads the requested model if missing.
  celltypist$models$download_models(model = model,
                                     force_update = FALSE)
  pred <- celltypist$annotate(adata,
                              model            = model,
                              majority_voting  = majority_voting)

  # `pred$predicted_labels` is a pandas DataFrame with columns
  # 'predicted_labels', 'over_clustering', 'majority_voting' (when MV=T),
  # and 'conf_score'.
  labels_df <- reticulate::py_to_r(pred$predicted_labels)
  labels_df$cell <- rownames(labels_df)

  # Choose the majority-vote column when it exists, otherwise per-cell.
  final_label <- if (isTRUE(majority_voting) &&
                     "majority_voting" %in% colnames(labels_df)) {
    labels_df$majority_voting
  } else {
    labels_df$predicted_labels
  }
  scores <- labels_df$conf_score
  if (is.null(scores)) scores <- rep(NA_real_, length(final_label))

  ix <- match(colnames(obj), labels_df$cell)
  obj@meta.data[[new_col]]                  <- as.character(final_label[ix])
  obj@meta.data[[paste0(new_col, "_score")]] <- as.numeric(scores[ix])
  obj
}


# ============================================================================
# scANVI backend
# ============================================================================
#' @keywords internal
#' @noRd
.annotate_scanvi <- function(obj, assay, new_col, args) {

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("'reticulate' is required for method = 'scanvi'.")
  }
  if (!is.null(args$python)) reticulate::use_python(args$python,
                                                    required = TRUE)
  if (!reticulate::py_module_available("scvi")) {
    stop("Python module 'scvi' (scvi-tools) not available. Install with ",
         "`reticulate::py_install(\"scvi-tools\")` or `pip install scvi-tools`.")
  }

  reference     <- args$reference     %||% stop("`reference` is required.")
  ref_label_col <- args$ref_label_col %||% stop("`ref_label_col` is required.")
  batch_col     <- args$batch_col     %||% NULL
  n_latent      <- args$n_latent      %||% 30L
  max_epochs_scvi   <- args$max_epochs_scvi   %||% 400L
  max_epochs_scanvi <- args$max_epochs_scanvi %||% 20L

  # Write query + reference to h5ad, load in Python, concatenate, train,
  # transfer labels.
  h5ad_q <- tempfile(pattern = "sctools_query_", fileext = ".h5ad")
  h5ad_r <- tempfile(pattern = "sctools_ref_",   fileext = ".h5ad")
  on.exit(unlink(c(h5ad_q, h5ad_r)), add = TRUE)

  if (!exists("ToAnnData", mode = "function")) {
    stop("ToAnnData() not available; cannot convert Seurat to AnnData ",
         "for scANVI.")
  }
  ToAnnData(obj, file = h5ad_q, assay = assay, layer = "counts")
  if (inherits(reference, "Seurat")) {
    ToAnnData(reference, file = h5ad_r,
              assay = Seurat::DefaultAssay(reference), layer = "counts")
  } else if (is.character(reference) && file.exists(reference)) {
    file.copy(reference, h5ad_r, overwrite = TRUE)
  } else {
    stop("`reference` must be a Seurat object or a path to an .h5ad file.")
  }

  scvi    <- reticulate::import("scvi",    delay_load = FALSE)
  ad_mod  <- reticulate::import("anndata", delay_load = FALSE)
  np      <- reticulate::import("numpy",   delay_load = FALSE)

  message("  Loading query and reference AnnData...")
  q <- ad_mod$read_h5ad(h5ad_q)
  r <- ad_mod$read_h5ad(h5ad_r)

  # Tag data source
  q$obs["_source"] <- "query"
  r$obs["_source"] <- "reference"

  # Ensure the label column exists on both; query is unlabeled.
  # scANVI convention: use `unlabeled_category = "Unknown"`.
  if (!(ref_label_col %in% colnames(r$obs))) {
    stop("`ref_label_col` '", ref_label_col, "' not found on reference.")
  }
  q$obs[ref_label_col] <- "Unknown"

  # Concatenate on shared genes
  combined <- ad_mod$concat(
    reticulate::dict(reference = r, query = q),
    join = "inner", label = "_dataset"
  )

  # Batch key
  if (!is.null(batch_col) && (batch_col %in% colnames(combined$obs))) {
    scvi$model$SCVI$setup_anndata(combined, batch_key = batch_col,
                                  labels_key = ref_label_col)
  } else {
    scvi$model$SCVI$setup_anndata(combined, labels_key = ref_label_col)
  }

  message("  Training scVI...")
  vae <- scvi$model$SCVI(combined, n_latent = as.integer(n_latent))
  vae$train(max_epochs = as.integer(max_epochs_scvi))

  message("  Training scANVI on top of scVI...")
  lvae <- scvi$model$SCANVI$from_scvi_model(vae,
                                            unlabeled_category = "Unknown",
                                            labels_key = ref_label_col)
  lvae$train(max_epochs = as.integer(max_epochs_scanvi))

  # Predict on the query subset
  preds       <- lvae$predict(soft = FALSE)
  soft        <- lvae$predict(soft = TRUE)   # a pandas DataFrame
  is_query    <- reticulate::py_to_r(combined$obs["_source"] == "query")

  preds_r <- reticulate::py_to_r(preds)
  soft_r  <- as.data.frame(reticulate::py_to_r(soft))
  # Confidence = predicted-class probability
  max_prob <- apply(soft_r, 1, max, na.rm = TRUE)

  # Match cell order back to Seurat columns
  cell_names <- as.character(reticulate::py_to_r(combined$obs_names$values))
  df <- data.frame(cell   = cell_names,
                   label  = as.character(preds_r),
                   score  = as.numeric(max_prob),
                   source = as.character(reticulate::py_to_r(
                     combined$obs["_source"]$values)),
                   stringsAsFactors = FALSE)
  df <- df[df$source == "query", ]
  ix <- match(colnames(obj), df$cell)

  obj@meta.data[[new_col]]                  <- df$label[ix]
  obj@meta.data[[paste0(new_col, "_score")]] <- df$score[ix]
  obj
}


# ============================================================================
# scmap backend
# ============================================================================
#' @keywords internal
#' @noRd
.annotate_scmap <- function(obj, assay, new_col, args) {

  if (!requireNamespace("scmap", quietly = TRUE)) {
    stop("'scmap' is required. Install with BiocManager::install('scmap').")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("'SingleCellExperiment' is required.")
  }

  reference     <- args$reference     %||% stop("`reference` is required.")
  ref_label_col <- args$ref_label_col %||% stop("`ref_label_col` is required.")
  scmap_method  <- args$scmap_method  %||% "cluster"
  threshold     <- args$threshold     %||% 0.5
  n_features    <- args$n_features    %||% 500L

  # Build SCE for reference
  if (inherits(reference, "Seurat")) {
    ref_a <- Seurat::DefaultAssay(reference)
    ref_counts <- SeuratObject::LayerData(reference, assay = ref_a,
                                          layer = "counts")
    ref_data   <- tryCatch(
      SeuratObject::LayerData(reference, assay = ref_a, layer = "data"),
      error = function(e) NULL
    )
    if (is.null(ref_data)) ref_data <- log1p(ref_counts)
    ref_sce <- SingleCellExperiment::SingleCellExperiment(
      assays  = list(counts = ref_counts, logcounts = ref_data),
      colData = reference@meta.data
    )
    SummarizedExperiment::rowData(ref_sce)$feature_symbol <- rownames(ref_sce)
    ref_sce$cell_type1 <- as.character(reference@meta.data[[ref_label_col]])
  } else if (inherits(reference, "SingleCellExperiment")) {
    ref_sce <- reference
    if (is.null(SummarizedExperiment::rowData(ref_sce)$feature_symbol)) {
      SummarizedExperiment::rowData(ref_sce)$feature_symbol <- rownames(ref_sce)
    }
    ref_sce$cell_type1 <- as.character(
      SummarizedExperiment::colData(ref_sce)[[ref_label_col]])
  } else {
    stop("`reference` must be a Seurat or SingleCellExperiment.")
  }

  # Build SCE for query
  q_counts <- SeuratObject::LayerData(obj, assay = assay, layer = "counts")
  q_data   <- tryCatch(
    SeuratObject::LayerData(obj, assay = assay, layer = "data"),
    error = function(e) NULL
  )
  if (is.null(q_data)) q_data <- log1p(q_counts)
  q_sce <- SingleCellExperiment::SingleCellExperiment(
    assays  = list(counts = q_counts, logcounts = q_data),
    colData = obj@meta.data
  )
  SummarizedExperiment::rowData(q_sce)$feature_symbol <- rownames(q_sce)

  # Restrict to shared genes (scmap does this internally too, but be explicit)
  shared <- intersect(rownames(ref_sce), rownames(q_sce))
  if (length(shared) < 100) {
    warning("Only ", length(shared), " gene(s) shared between ",
            "reference and query; scmap projection may be unreliable.",
            call. = FALSE)
  }
  ref_sce <- ref_sce[shared, ]
  q_sce   <- q_sce[shared, ]

  message(sprintf("  Selecting %d informative features for scmap...",
                  n_features))
  ref_sce <- scmap::selectFeatures(ref_sce, n_features = as.integer(n_features),
                                   suppress_plot = TRUE)

  if (scmap_method == "cluster") {
    ref_sce <- scmap::indexCluster(ref_sce, cluster_col = "cell_type1")
    proj <- scmap::scmapCluster(
      projection      = q_sce,
      index_list      = list(ref = SingleCellExperiment::metadata(ref_sce)$scmap_cluster_index),
      threshold       = threshold
    )
    labels <- proj$scmap_cluster_labs[, 1]
    scores <- proj$scmap_cluster_siml[, 1]
  } else {
    ref_sce <- scmap::indexCell(ref_sce)
    proj <- scmap::scmapCell(
      projection = q_sce,
      index_list = list(ref = SingleCellExperiment::metadata(ref_sce)$scmap_cell_index)
    )
    # Vote across the returned nearest neighbors
    # (proj$ref$cells is a nearest-neighbor index matrix)
    nn <- proj$ref$cells
    votes <- apply(nn, 2, function(inds) {
      labs <- as.character(ref_sce$cell_type1)[inds]
      tab <- sort(table(labs), decreasing = TRUE)
      names(tab)[1]
    })
    sims <- apply(proj$ref$similarities, 2, max)
    below <- sims < threshold
    votes[below] <- "unassigned"
    labels <- votes
    scores <- sims
  }

  obj@meta.data[[new_col]]                  <- as.character(labels)
  obj@meta.data[[paste0(new_col, "_score")]] <- as.numeric(scores)
  obj
}


# Small "or default" utility so the switch branches read cleanly.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b
