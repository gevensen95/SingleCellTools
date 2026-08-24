#' Call doublets on a Seurat object with DoubletFinder
#'
#' Runs the standard DoubletFinder workflow on a single Seurat object:
#' normalize → PCA → choose # PCs → clusters → pK sweep → doubletFinder.
#'
#' On return, the normalized \code{data} and \code{scale.data} layers, the
#' variable features, and the \code{pca} reduction created during the
#' workflow are stripped from the object so the result carries only counts
#' plus the new \code{doublet_finder} (and \code{seurat_clusters}) metadata
#' columns.
#'
#' @param obj A Seurat object.
#' @param samplenameIndex Passed through unchanged from the original API
#'   (currently unused inside the function; retained for call-site
#'   compatibility).
#' @param normalization One of \code{"LogNormalize"} (default) or \code{"SCT"}.
#'   Selects the \code{NormalizeData / FindVariableFeatures / ScaleData}
#'   pipeline vs \code{SCTransform} AND flips the \code{sct} flag passed to
#'   DoubletFinder's \code{paramSweep} and \code{doubletFinder}.
#' @param vars.to.regress Optional character vector of variables to regress
#'   out during normalization (passed to either \code{SCTransform} or
#'   \code{ScaleData}). Default \code{NULL}.
#' @param cluster_resolution Resolution passed to \code{FindClusters}.
#'   Default \code{0.1}.
#' @param doublet_rate Assumed doublet \emph{formation} rate for this library
#'   -- the fraction of barcodes expected to be multiplets given how many
#'   cells were loaded. Default \code{0.075} (the value DoubletFinder's own
#'   vignette uses). This sets \code{nExp}, the number of doublets the
#'   classifier is told to find, via
#'   \code{round(doublet_rate * ncol(obj) * (1 - homotypic_proportion))}.
#'
#'   For 10X Chromium the rate scales roughly linearly with recovered cells
#'   (about 0.8\% per 1,000 cells recovered), so a run recovering ~5,000
#'   cells is nearer \code{0.04} and one recovering ~10,000 nearer
#'   \code{0.08}; consult the loading table for your chemistry rather than
#'   accepting the default if your recoveries vary a lot between samples.
#'
#'   Note this is deliberately \emph{not} \code{pK}. \code{pK} is the
#'   neighborhood size the pANN statistic is computed over -- it is chosen
#'   from the data by \code{find.pK()} and says nothing about how many
#'   doublets the run actually produced. Earlier versions of this function
#'   used \code{optimal.pk} here, which made each sample's assumed doublet
#'   rate an artifact of where its BCmvn curve happened to peak: across one
#'   8-sample experiment that produced called-doublet rates ranging from
#'   0.4\% to 18.4\%, rank-ordered by pK rather than by anything biological.
#'   Objects created before this change carry doublet calls made that way
#'   and should be re-run.
#' @param pk_sweep_max_cells Maximum number of cells used to estimate pK via
#'   \code{DoubletFinder::paramSweep}. \code{paramSweep} reprocesses a fresh
#'   real+artificial-doublet matrix from scratch for each of its 6 fixed pN
#'   values, which is normally the dominant cost of this function -- that
#'   cost scales worse than linearly with cell count. When the object has
#'   more than \code{pk_sweep_max_cells} cells, pK is estimated from a
#'   random subsample of that size (subset from the already-normalized/
#'   scaled/PCA'd object, so no preprocessing is redone); homotypic
#'   proportion, \code{nExp}, and the final doublet classification all still
#'   use every cell in \code{obj}. This trades a small amount of pK-estimate
#'   precision for speed -- set to \code{Inf} (or higher than your largest
#'   sample) to disable and always sweep on the full object. Default
#'   \code{4000}.
#' @param sweep_cores Number of cores \code{DoubletFinder::paramSweep} uses
#'   internally to parallelize across its 6 pN values (its own
#'   \code{num.cores} argument, fork-based -- Mac/Linux only). Default
#'   \code{1} (sequential). If \code{calldoublet} is itself being called
#'   from several \code{workers} in parallel (e.g. via
#'   \code{CreateRNAObjects}'s \code{workers} argument), raising this without
#'   lowering the outer worker count oversubscribes the CPU (\code{workers *
#'   sweep_cores} concurrent processes) and can make things slower, not
#'   faster.
#' @return The input Seurat object with a \code{doublet_finder} metadata
#'   column ("Doublet" / "Singlet").
#' @export
calldoublet <- function(obj,
                        samplenameIndex,
                        normalization      = c("LogNormalize", "SCT"),
                        vars.to.regress    = NULL,
                        cluster_resolution = 0.1,
                        doublet_rate       = 0.075,
                        pk_sweep_max_cells = 4000,
                        sweep_cores        = 1) {
  normalization <- match.arg(normalization)
  use_sct       <- identical(normalization, "SCT")

  if (!is.numeric(doublet_rate) || length(doublet_rate) != 1L ||
      is.na(doublet_rate) || doublet_rate < 0 || doublet_rate >= 1) {
    stop("`doublet_rate` must be a single number in [0, 1) -- got: ",
         paste(utils::capture.output(dput(doublet_rate)), collapse = ""))
  }
  if (doublet_rate > 0.3) {
    warning("`doublet_rate` = ", doublet_rate, " is very high; 10X runs are ",
            "typically 0.01-0.10. Check you haven't passed a pK value here.")
  }

  # ---- Normalize ----------------------------------------------------------
  if (use_sct) {
    obj <- Seurat::SCTransform(obj, vars.to.regress = vars.to.regress,
                               verbose = FALSE)
  } else {
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE)
    # Explicitly restrict to variable features rather than relying on
    # ScaleData()'s default (which, depending on Seurat version, may scale
    # every gene) -- scaling ~2000 genes instead of the whole transcriptome
    # is a meaningful speedup with no effect on the PCA used downstream.
    obj <- Seurat::ScaleData(obj, vars.to.regress = vars.to.regress,
                             features = SeuratObject::VariableFeatures(obj),
                             verbose = FALSE)
  }

  obj <- Seurat::RunPCA(obj, verbose = FALSE)

  # ---- Find significant PCs -----------------------------------------------
  stdv         <- obj[["pca"]]@stdev
  percent.stdv <- (stdv / sum(stdv)) * 100
  cumulative   <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:(length(percent.stdv) - 1)] -
                       percent.stdv[2:length(percent.stdv)]) > 0.1),
              decreasing = TRUE)[1] + 1
  # Fall back to whichever heuristic DID produce an answer. `min(co1, co2)`
  # returns NA the moment either is NA, so a sample where (say) the
  # cumulative-variance rule found a cutoff but no consecutive stdev drop
  # cleared 0.1 used to abort outright -- despite having a perfectly usable
  # co1 in hand. Only genuinely erroring when BOTH rules come up empty
  # matches what the error message below already claims to distinguish.
  min.pc <- suppressWarnings(min(c(co1, co2), na.rm = TRUE))
  if (!is.finite(min.pc)) {
    stop(
      "calldoublet(): could not determine a number of significant PCs -- ",
      if (is.na(co1)) "no PC reached >90% cumulative variance with <5% individual contribution" else "",
      if (is.na(co1) && is.na(co2)) " and " else "",
      if (is.na(co2)) "no consecutive PC-to-PC stdev drop exceeded 0.1" else "",
      ". This usually means too few informative PCs were computed for this ",
      "sample (check cell count vs. RunPCA()'s npcs) -- inspect obj[['pca']]@stdev directly."
    )
  }

  # ---- Finish pre-processing ----------------------------------------------
  # NB: irlba and RSpectra are pulled in via the package Imports so that the
  # PCA code that needs them can find them — no library() calls here.
  # (No RunUMAP() here: DoubletFinder only needs the PCA space and cluster
  # labels below, and nothing else in this function ever reads a UMAP
  # embedding -- computing one was pure wasted work.)
  obj <- Seurat::FindNeighbors(obj, dims = 1:min.pc, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = cluster_resolution, verbose = FALSE)

  # ---- pK identification (no ground-truth) --------------------------------
  # paramSweep() rebuilds a real+artificial-doublet matrix and reruns
  # Normalize/FindVariableFeatures/Scale/PCA from scratch for each of its 6
  # fixed pN values -- that cost scales worse than linearly with cell count,
  # so estimating pK from a capped-size random subsample (rather than every
  # cell) is a meaningful speedup. subset() on an object that's already been
  # normalized/scaled/PCA'd just slices the existing layers/reduction, it
  # doesn't redo any of that work -- so building obj_sweep here is cheap.
  # Everything downstream of the pK estimate (homotypic proportion, nExp,
  # and the final doubletFinder() classification) still uses every cell in
  # `obj`, unaffected by this subsampling.
  n_cells <- ncol(obj)
  if (is.finite(pk_sweep_max_cells) && n_cells > pk_sweep_max_cells) {
    sweep_cells <- sample(colnames(obj), pk_sweep_max_cells)
    obj_sweep   <- subset(obj, cells = sweep_cells)
  } else {
    obj_sweep <- obj
  }

  sweep.list  <- DoubletFinder::paramSweep(obj_sweep, PCs = 1:min.pc, sct = use_sct,
                                           num.cores = sweep_cores)

  # summarizeSweep() runs a gaussian KDE (KernSmooth::bkde) over each grid
  # point's pANN vector to score its bimodality. If the pANN values at some
  # (pN, pK) are all identical the KDE bandwidth is undefined, and the
  # failure surfaces from deep inside KernSmooth as the thoroughly
  # unhelpful "missing value where TRUE/FALSE needed" -- naming neither
  # DoubletFinder, nor pANN, nor the sample being processed. Translate it.
  sweep.stats <- tryCatch(
    DoubletFinder::summarizeSweep(sweep.list),
    error = function(e) {
      stop("calldoublet(): DoubletFinder's pK sweep failed while scoring ",
           "bimodality -- underlying error: ", conditionMessage(e), "\n",
           "This usually means the pANN distribution was degenerate at one ",
           "or more (pN, pK) grid points: every cell had the same fraction ",
           "of artificial-doublet neighbours, so the kernel density estimate ",
           "has no bandwidth to work with. It shows up when a sample has very ",
           "few distinct populations, unusually sharp separation between ",
           "them, or too few cells in the sweep. Try raising ",
           "pk_sweep_max_cells (currently ",
           format(pk_sweep_max_cells), ") so the sweep sees more of the ",
           "sample's heterogeneity.",
           call. = FALSE)
    }
  )

  bcmvn       <- DoubletFinder::find.pK(sweep.stats)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  bcmvn.max   <- bcmvn[which.max(bcmvn$BCmetric), ]
  optimal.pk  <- bcmvn.max$pK
  optimal.pk  <- as.numeric(levels(optimal.pk))[optimal.pk]

  # ---- Homotypic doublet proportion estimate ------------------------------
  # nExp is driven by `doublet_rate` (how many multiplets the LOADING is
  # expected to have produced), NOT by optimal.pk. pK is the pANN
  # neighborhood size that find.pK() selected from the data; it carries no
  # information about doublet abundance, so using it here made each sample's
  # assumed doublet rate a side effect of its own BCmvn curve shape. See the
  # `doublet_rate` parameter docs.
  annotations    <- obj@meta.data$seurat_clusters
  homotypic.prop <- DoubletFinder::modelHomotypic(annotations)
  nExp.poi       <- round(doublet_rate * nrow(obj@meta.data))
  nExp.poi.adj   <- round(nExp.poi * (1 - homotypic.prop))

  # ---- Run DoubletFinder --------------------------------------------------
  obj <- DoubletFinder::doubletFinder(seu  = obj,
                                      PCs  = 1:min.pc,
                                      pK   = optimal.pk,
                                      nExp = nExp.poi.adj,
                                      sct  = use_sct)

  # doubletFinder() adds TWO columns: the pANN score and the classification,
  # both with the run's pN/pK/nExp values baked into their names. Rename both
  # by PATTERN rather than by position.
  #
  # Renaming the classification by position ("the last column") happened to
  # work, but left the pANN column named e.g. "pANN_0.25_0.24_881" -- and
  # since pK and nExp differ per sample, every object in a multi-sample list
  # ended up with a uniquely-named metadata column. Any downstream function
  # that stacks metadata across samples (QCComparePlots, MergeSeurat,
  # CellSuiteSummary) then fails in rbind.data.frame() with the distinctly
  # unhelpful "names do not match previous names".
  metadata <- obj@meta.data

  df_col <- grep("^DF\\.classifications", colnames(metadata))
  if (length(df_col) != 1L) {
    stop("calldoublet(): expected exactly one 'DF.classifications*' column ",
         "after doubletFinder(), found ", length(df_col),
         ". DoubletFinder's output naming may have changed -- got: ",
         paste(colnames(metadata), collapse = ", "))
  }
  colnames(metadata)[df_col] <- "doublet_finder"

  pann_col <- grep("^pANN_", colnames(metadata))
  if (length(pann_col) == 1L) {
    colnames(metadata)[pann_col] <- "doublet_pANN"
  }

  obj@meta.data <- metadata

  cat("Doublet estimation with DoubletFinder() ...\n")
  cat("Normalization: ", normalization, "\n")
  # Named lookup rather than positional [1]/[2] -- table()'s alphabetical
  # default ordering happens to put "Doublet" before "Singlet" today, but
  # relying on that silently mislabels these counts if it ever doesn't.
  doublet_tab <- table(metadata[["doublet_finder"]])
  cat("Doublet: ", if ("Doublet" %in% names(doublet_tab)) doublet_tab[["Doublet"]] else 0, "\n")
  cat("Singlet: ", if ("Singlet" %in% names(doublet_tab)) doublet_tab[["Singlet"]] else 0, "\n")

  # ---- Cleanup: strip everything we created during the workflow -----------
  # For SCT, drop the entire SCT assay — all of its derived layers go with it.
  if (use_sct && "SCT" %in% SeuratObject::Assays(obj)) {
    SeuratObject::DefaultAssay(obj) <- 'RNA'
    obj[['SCT']]                    <- NULL
  }

  # Drop the normalized 'data' and 'scale.data' layers from the RNA assay,
  # including v5 split variants (data.1, scale.data.2, ...).
  #
  # The deletion MUST go through LayerData()<-. The obvious-looking
  # `obj[["RNA"]][["data"]] <- NULL` is a SILENT NO-OP on a v5 Assay5 -- it
  # neither errors nor removes anything, so this whole cleanup step did
  # nothing and every object this function returned still carried its
  # normalized layers, directly contradicting the promise in the docs above.
  # Those layers are dense-ish and typically dwarf the counts they sit
  # beside, so this was the difference between a lean counts-only object and
  # one several times larger, silently, on every sample.
  drop_lyrs <- grep("^(data|scale\\.data)(\\.|$)",
                    SeuratObject::Layers(obj[["RNA"]]), value = TRUE)
  for (lyr in drop_lyrs) {
    SeuratObject::LayerData(obj, assay = "RNA", layer = lyr) <- NULL
  }

  # Clear the variable-features set chosen by FindVariableFeatures /
  # SCTransform. tryCatch keeps this safe across Seurat versions where the
  # setter signature has varied.
  tryCatch({
    SeuratObject::VariableFeatures(obj, assay = "RNA") <- character(0)
  }, error = function(e) invisible(NULL))

  # Drop the PCA reduction we computed for the pK sweep.
  for (red in intersect("pca", SeuratObject::Reductions(obj))) {
    obj[[red]] <- NULL
  }

  return(obj)
}
