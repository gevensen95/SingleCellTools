#!/usr/bin/env Rscript
# ============================================================================
# validate_combine_rounds.R
#
# Decide, with evidence, whether combining two sequencing rounds is worth it
# -- and whether the summed-counts approximation used by
# CombineCellRangerRounds() / CombineParseRounds() is accurate enough to trust.
#
# It answers three separate questions. They are independent; read them as
# three verdicts, not one.
#
#   Q1  Is round 2 buying anything?          -> saturation + genes/cell gain
#   Q2  How wrong is the summed combine?     -> exact molecule-level overlap
#   Q3  Does it change the biology?          -> pseudobulk + cluster stability
#
# NOTE ON THE GOLD STANDARD
# If you still have the FASTQs, none of this is necessary for Q2. Run
#     cellranger count --fastqs=/round1,/round2 --sample=X ...
# (or `cat` both rounds' FASTQs and run split-pipe once, for Parse) and you
# get an exact combine -- one alignment, one dedup, one cell-calling pass.
# This script exists for when that isn't available, and to answer Q1/Q3,
# which the FASTQ route does not.
#
# Usage:  Rscript validate_combine_rounds.R
#         (edit the CONFIG block first)
# ============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

# ---- CONFIG ----------------------------------------------------------------

PLATFORM <- "cellranger"        # "cellranger" or "parse"

SAMPLES <- c("Alt1", "Alt2", "Alt3", "PBS1", "PBS2", "PN1", "PN2", "PN3")

ROUND1_DIR   <- "/path/to/round1"    # <SAMPLE> subdirectory per sample
ROUND2_DIR   <- "/path/to/round2"
COMBINED_DIR <- "/path/to/combined"  # output_dir you passed to Combine*Rounds()

OUT_PREFIX      <- "combine_validation"
RUN_CLUSTERING  <- TRUE   # Q3 cluster stability; the slow part. FALSE to skip.
MIN_CELLS_CLUST <- 200    # skip clustering for samples smaller than this

r1_path  <- function(s) file.path(ROUND1_DIR, s)
r2_path  <- function(s) file.path(ROUND2_DIR, s)
cmb_path <- function(s) file.path(COMBINED_DIR, s)

# ---- Readers ---------------------------------------------------------------

read_counts <- function(path) {
  if (PLATFORM == "parse") {
    dge <- file.path(path, "DGE_filtered")
    if (!dir.exists(dge)) return(NULL)
    m <- Matrix::t(Matrix::readMM(file.path(dge, "count_matrix.mtx")))
    g <- utils::read.csv(file.path(dge, "all_genes.csv"), stringsAsFactors = FALSE)
    c <- utils::read.csv(file.path(dge, "cell_metadata.csv"), stringsAsFactors = FALSE)
    dimnames(m) <- list(if ("gene_id" %in% names(g)) g$gene_id else g$gene_name,
                        c$bc_wells)
    return(methods::as(m, "CsparseMatrix"))
  }
  # Explicit ReadMtx() rather than Read10X(): Seurat 5's Read10X() insists on
  # gzipped files, and CombineCellRangerRounds() writes its combined triplet
  # UNCOMPRESSED. Read10X() on the combined output fails with
  # "Barcode file missing. Expecting barcodes.tsv.gz". Matching either form
  # here mirrors what the package's own .read_10x_triplet() does.
  pick <- function(d, pat) {
    f <- list.files(d, pattern = pat, full.names = TRUE)
    if (length(f)) f[1] else NA_character_
  }
  for (d in c(path,
              file.path(path, "filtered_feature_bc_matrix"),
              file.path(path, "outs", "filtered_feature_bc_matrix"))) {
    if (!dir.exists(d)) next
    mtx  <- pick(d, "matrix\\.mtx(\\.gz)?$")
    bcs  <- pick(d, "barcodes\\.tsv(\\.gz)?$")
    fts  <- pick(d, "(features|genes)\\.tsv(\\.gz)?$")
    if (anyNA(c(mtx, bcs, fts))) next
    return(tryCatch(
      Seurat::ReadMtx(mtx = mtx, cells = bcs, features = fts,
                      feature.column = 2),
      error = function(e) {
        message("  reader failed at ", d, ": ", conditionMessage(e)); NULL
      }))
  }
  NULL
}

as_obj <- function(m, project) {
  if (is.null(m)) return(NULL)
  o <- suppressWarnings(
    Seurat::CreateSeuratObject(counts = methods::as(m, "CsparseMatrix"),
                               project = project))
  mt <- grep("^(MT|mt)-", rownames(o), value = TRUE)
  o$percent.mt <- if (length(mt)) Seurat::PercentageFeatureSet(o, features = mt) else 0
  o
}

# ============================================================================
# Q1  Is round 2 buying anything?
# ============================================================================
# Two signals. Sequencing saturation says how much of the library each round
# had already seen -- a round-1 saturation of 80%+ means round 2 is mostly
# re-reading molecules you already had. Median genes/cell says what the extra
# depth actually recovered. A big saturation number with a small gene gain is
# the clearest possible "don't bother".

read_saturation <- function(path) {
  f <- c(file.path(path, "outs", "metrics_summary.csv"),
         file.path(path, "metrics_summary.csv"))
  f <- f[file.exists(f)][1]
  if (is.na(f)) return(NA_real_)
  m <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  col <- grep("Sequencing Saturation", colnames(m), value = TRUE)[1]
  if (is.na(col)) return(NA_real_)
  as.numeric(sub("%", "", as.character(m[[col]][1])))
}

message("\n=== Q1: is round 2 buying anything? ===")

q1 <- do.call(rbind, lapply(SAMPLES, function(s) {
  m1 <- read_counts(r1_path(s)); m2 <- read_counts(r2_path(s))
  mc <- read_counts(cmb_path(s))
  if (is.null(m1) || is.null(mc)) {
    message(sprintf("  [%s] missing round1 or combined matrix -- skipping", s))
    return(NULL)
  }
  g1 <- Matrix::colSums(m1 > 0); gc_ <- Matrix::colSums(mc > 0)
  u1 <- Matrix::colSums(m1);     uc  <- Matrix::colSums(mc)
  data.frame(
    sample        = s,
    sat_r1_pct    = read_saturation(r1_path(s)),
    sat_r2_pct    = read_saturation(r2_path(s)),
    cells_r1      = ncol(m1),
    cells_r2      = if (is.null(m2)) NA_integer_ else ncol(m2),
    cells_comb    = ncol(mc),
    med_genes_r1  = median(g1),
    med_genes_cb  = median(gc_),
    gene_gain_pct = round(100 * (median(gc_) / median(g1) - 1), 1),
    med_umi_r1    = median(u1),
    med_umi_cb    = median(uc),
    umi_gain_pct  = round(100 * (median(uc) / median(u1) - 1), 1),
    stringsAsFactors = FALSE
  )
}))
if (is.null(q1)) message("  no sample yielded both a round-1 and a combined matrix.") else
  print(q1, row.names = FALSE)

# ============================================================================
# Q2  How wrong is the summed combine?
# ============================================================================
# Summing two post-dedup matrices counts twice any molecule that was
# sequenced in BOTH rounds. molecule_info.h5 records every
# (barcode, UMI, feature) tuple, so the overlap is directly measurable --
# no re-alignment needed. Restricted to barcodes present in both rounds'
# filtered sets, which is both the population that matters and a large
# memory saving.

message("\n=== Q2: how much does summing double-count? ===")

molecules_for <- function(path, keep_bcs) {
  f <- c(file.path(path, "outs", "molecule_info.h5"),
         file.path(path, "molecule_info.h5"))
  f <- f[file.exists(f)][1]
  if (is.na(f)) return(NULL)

  bcs <- as.character(rhdf5::h5read(f, "barcodes"))
  fts <- as.character(rhdf5::h5read(f, "features/id"))
  bi  <- as.integer(rhdf5::h5read(f, "barcode_idx")) + 1L   # h5 is 0-based
  fi  <- as.integer(rhdf5::h5read(f, "feature_idx")) + 1L
  umi <- rhdf5::h5read(f, "umi")

  bc_chr <- sub("-\\d+$", "", bcs[bi])
  sel    <- bc_chr %in% keep_bcs
  paste(bc_chr[sel], umi[sel], fts[fi[sel]], sep = "|")
}

if (PLATFORM != "cellranger") {
  message("  Q2 needs CellRanger molecule_info.h5; skipping for Parse.")
  message("  Parse equivalent: pull CB/UB/gene tags from each round's BAM.")
  q2 <- NULL
} else if (!requireNamespace("rhdf5", quietly = TRUE)) {
  message("  {rhdf5} not installed -- skipping Q2.")
  message("  BiocManager::install('rhdf5')")
  q2 <- NULL
} else {
  q2 <- do.call(rbind, lapply(SAMPLES, function(s) {
    m1 <- read_counts(r1_path(s)); m2 <- read_counts(r2_path(s))
    if (is.null(m1) || is.null(m2)) return(NULL)
    shared <- intersect(sub("-\\d+$", "", colnames(m1)),
                        sub("-\\d+$", "", colnames(m2)))
    if (!length(shared)) return(NULL)

    mo1 <- molecules_for(r1_path(s), shared)
    mo2 <- molecules_for(r2_path(s), shared)
    if (is.null(mo1) || is.null(mo2)) {
      message(sprintf("  [%s] no molecule_info.h5 -- skipping", s)); return(NULL)
    }
    dup    <- length(intersect(mo1, mo2))
    summed <- length(mo1) + length(mo2)
    rm(mo1, mo2); invisible(gc())
    data.frame(sample = s, shared_cells = length(shared),
               summed_molecules = summed, double_counted = dup,
               inflation_pct = round(100 * dup / (summed - dup), 2),
               stringsAsFactors = FALSE)
  }))
  if (!is.null(q2)) print(q2, row.names = FALSE)
}

# ============================================================================
# Q3  Does combining change the biology?
# ============================================================================
# Two checks. Pseudobulk concordance asks whether the expression profile
# moved at all. Cluster stability asks the question that actually matters:
# would you draw different conclusions? If round-1-only and combined give the
# same clusters over the same cells, round 2 is not changing your answer and
# round 1 alone is the simpler, assumption-free choice.

message("\n=== Q3: does it change the biology? ===")

quick_cluster <- function(o, seed = 1) {
  o <- Seurat::NormalizeData(o, verbose = FALSE)
  o <- Seurat::FindVariableFeatures(o, verbose = FALSE)
  o <- Seurat::ScaleData(o, features = Seurat::VariableFeatures(o), verbose = FALSE)
  o <- Seurat::RunPCA(o, npcs = 30, verbose = FALSE, seed.use = seed)
  o <- Seurat::FindNeighbors(o, dims = 1:30, verbose = FALSE)
  o <- Seurat::FindClusters(o, resolution = 0.5, verbose = FALSE, random.seed = seed)
  o
}

q3 <- do.call(rbind, lapply(SAMPLES, function(s) {
  m1 <- read_counts(r1_path(s)); mc <- read_counts(cmb_path(s))
  if (is.null(m1) || is.null(mc)) return(NULL)

  # --- pseudobulk concordance on shared genes -------------------------------
  g  <- intersect(rownames(m1), rownames(mc))
  p1 <- log1p(Matrix::rowSums(m1[g, , drop = FALSE]) /
                sum(m1[g, , drop = FALSE]) * 1e6)
  pc <- log1p(Matrix::rowSums(mc[g, , drop = FALSE]) /
                sum(mc[g, , drop = FALSE]) * 1e6)
  r  <- stats::cor(p1, pc, method = "pearson")

  # --- round-2-only barcodes: real cells, or debris? ------------------------
  b1  <- sub("-\\d+$", "", colnames(m1))
  bc  <- sub("-\\d+$", "", colnames(mc))
  new <- setdiff(bc, b1)
  new_med_umi    <- if (length(new)) median(Matrix::colSums(mc)[match(new, bc)]) else NA_real_
  shared_med_umi <- median(Matrix::colSums(mc)[match(intersect(bc, b1), bc)])

  # --- cluster stability over the cells both analyses contain ---------------
  ari <- NA_real_
  if (isTRUE(RUN_CLUSTERING) && ncol(m1) >= MIN_CELLS_CLUST &&
      requireNamespace("mclust", quietly = TRUE)) {
    o1 <- quick_cluster(as_obj(m1, "r1"))
    oc <- quick_cluster(as_obj(mc, "cb"))
    common <- intersect(sub("-\\d+$", "", colnames(o1)),
                        sub("-\\d+$", "", colnames(oc)))
    if (length(common) > 50) {
      l1 <- as.character(Seurat::Idents(o1))[match(common, sub("-\\d+$", "", colnames(o1)))]
      lc <- as.character(Seurat::Idents(oc))[match(common, sub("-\\d+$", "", colnames(oc)))]
      ari <- round(mclust::adjustedRandIndex(l1, lc), 3)
    }
    rm(o1, oc); invisible(gc())
  }

  data.frame(sample = s,
             pseudobulk_r      = round(r, 4),
             n_new_barcodes    = length(new),
             new_med_umi       = new_med_umi,
             shared_med_umi    = shared_med_umi,
             cluster_ARI       = ari,
             stringsAsFactors  = FALSE)
}))
if (is.null(q3)) message("  no sample yielded both a round-1 and a combined matrix.") else
  print(q3, row.names = FALSE)

# ============================================================================
# Verdict
# ============================================================================
# Deliberately prints the evidence rather than a yes/no -- the thresholds
# below are rules of thumb, not decision procedure. Look at the numbers.

message("\n=== READ THIS ===")

if (!is.null(q1)) {
  if (all(!is.na(q1$sat_r1_pct)) && median(q1$sat_r1_pct, na.rm = TRUE) > 80) {
    message("* Round 1 saturation is already high (median ",
            round(median(q1$sat_r1_pct, na.rm = TRUE), 1),
            "%). Round 2 is mostly re-reading molecules you had --")
    message("  small gain, and MAXIMAL double-counting. Both point at: don't combine.")
  }
  if (median(q1$gene_gain_pct, na.rm = TRUE) < 5) {
    message("* Median genes/cell gain is <5%. Round 2 is not buying detection.")
  }
}
if (!is.null(q2) && any(q2$inflation_pct > 5, na.rm = TRUE)) {
  message("* Some samples inflate >5% from double-counted molecules. Summing is")
  message("  too lossy there -- re-run from FASTQ, or union at the molecule level.")
}
if (!is.null(q3)) {
  if (any(q3$pseudobulk_r < 0.99, na.rm = TRUE)) {
    message("* Pseudobulk r < 0.99 somewhere -- the profile moved more than depth alone")
    message("  should explain. Check those samples are really the same library.")
  }
  if (any(q3$cluster_ARI < 0.8, na.rm = TRUE)) {
    message("* Cluster ARI < 0.8 somewhere: combining changes your conclusions.")
    message("  That justifies the approximation IF the combined result is the better one --")
    message("  confirm the extra structure is real biology, not depth artifact.")
  } else if (all(!is.na(q3$cluster_ARI))) {
    message("* Clusters are stable (all ARI >= 0.8): round 2 is not changing the")
    message("  biology. Round 1 alone is the simpler, assumption-free choice.")
  }
  bad <- q3$n_new_barcodes > 0 & q3$new_med_umi < 0.5 * q3$shared_med_umi
  if (any(bad, na.rm = TRUE)) {
    message("* Round-2-only barcodes have <half the UMIs of shared ones in: ",
            paste(q3$sample[which(bad)], collapse = ", "))
    message("  The union of independently-filtered cell lists is importing debris.")
  }
}

written <- character(0)
for (nm in c("q1", "q2", "q3")) {
  d <- get(nm)
  if (!is.null(d) && nrow(d)) {
    f <- sprintf("%s_%s.csv", OUT_PREFIX, nm)
    utils::write.csv(d, f, row.names = FALSE)
    written <- c(written, f)
  }
}
if (length(written)) {
  message("\nWrote: ", paste(written, collapse = ", "))
} else {
  message("\nNothing written -- no section produced results. Check the CONFIG paths.")
}
