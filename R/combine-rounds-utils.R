# Shared internal helpers for CombineParseRounds() and
# CombineCellRangerRounds().
#
# Everything in this file is format-agnostic: it operates on a generic
# per-round representation --
#   list(counts = <genes x cells sparse matrix>,
#        genes  = <data.frame with a gene-key column>,
#        cells  = <data.frame with a cell-id column>,
#        path   = <source path, for messages only>)
# -- with no assumption about which sequencing platform produced it. Each
# platform-specific Combine*Rounds() function is responsible for reading its
# own file layout into this shape (see .read_round() in CombineParseRounds.R
# and .read_round_cellranger() in CombineCellRangerRounds.R) and for writing
# the result back out in its own native format; this file only implements
# the actual "union genes/cells, sum overlapping barcodes, reconcile
# non-count metadata" logic once, so it isn't duplicated per platform.

# ---- Remap a sparse matrix's rows or columns to a new key set -------------
# Uses the sparse triplet representation so genes/cells absent from a round
# simply have no entries (i.e. are structural zeros) rather than
# materializing a dense intermediate.
#' @keywords internal
#' @noRd
.reindex <- function(mat, old_keys, new_keys, margin = c("rows", "cols")) {
  margin <- match.arg(margin)
  map <- match(old_keys, new_keys)
  # `new_keys` is always built as a union that contains every `old_key`, so an
  # NA here means a caller passed mismatched key sets. Catch it explicitly:
  # left alone, Matrix::sparseMatrix() would fail further down with an opaque
  # message about NAs in `i`/`j` that gives no hint which key set is wrong.
  if (anyNA(map)) {
    stop("Internal error in .reindex(): ", sum(is.na(map)),
         " value(s) in `old_keys` are absent from `new_keys`.")
  }
  trip <- methods::as(mat, "TsparseMatrix")
  if (margin == "rows") {
    new_i <- map[trip@i + 1L]
    Matrix::sparseMatrix(i = new_i, j = trip@j + 1L, x = trip@x,
                         dims = c(length(new_keys), ncol(mat)),
                         dimnames = list(new_keys, colnames(mat)))
  } else {
    new_j <- map[trip@j + 1L]
    Matrix::sparseMatrix(i = trip@i + 1L, j = new_j, x = trip@x,
                         dims = c(nrow(mat), length(new_keys)),
                         dimnames = list(rownames(mat), new_keys))
  }
}

# ---- Combine one sample's two rounds (format-agnostic core) ---------------
# `gene_key_col` / `cell_id_col` name the columns in r1$genes / r1$cells (and
# r2$genes / r2$cells) to join on -- e.g. "gene_id"/"bc_wells" for Parse,
# "gene_id"/"cell_id" for CellRanger. Callers are responsible for resolving
# which column that should be (Parse's gene_id-vs-gene_name fallback logic,
# for instance, is platform-specific and stays in CombineParseRounds.R) and
# for confirming both rounds actually have it before calling this.
#
# NB the same caveat documented on both public functions about summing two
# already-UMI-deduplicated rounds applies here identically regardless of
# platform -- this core makes no attempt to detect or correct for molecules
# sequenced (and counted) in both rounds' read sets.
#' @keywords internal
#' @noRd
.combine_one_core <- function(r1, r2, sample_name, gene_key_col, cell_id_col,
                              metadata_conflict = c("warn", "error")) {
  metadata_conflict <- match.arg(metadata_conflict)

  genes1_key <- r1$genes[[gene_key_col]]
  genes2_key <- r2$genes[[gene_key_col]]
  if (anyDuplicated(genes1_key) || anyDuplicated(genes2_key)) {
    warning(sample_name, ": '", gene_key_col, "' has duplicate values within ",
            "a round's gene table; gene-level alignment between rounds may ",
            "be unreliable for the affected genes.")
  }

  # Union genes (round 1 order first, then round-2-only genes appended).
  union_genes <- union(genes1_key, genes2_key)
  # Union cells (round 1 order first, then round-2-only barcodes appended).
  union_cells <- union(colnames(r1$counts), colnames(r2$counts))

  m1 <- .reindex(r1$counts, genes1_key, union_genes, margin = "rows")
  m1 <- .reindex(m1, colnames(m1), union_cells, margin = "cols")
  m2 <- .reindex(r2$counts, genes2_key, union_genes, margin = "rows")
  m2 <- .reindex(m2, colnames(m2), union_cells, margin = "cols")

  # This single addition IS the "sum overlapping, keep non-overlapping" rule:
  # a barcode/gene present in only one round has structural zeros in the
  # other round's reindexed matrix, so addition reproduces that round's
  # value exactly; a barcode/gene present in both has its counts summed.
  combined_counts <- m1 + m2

  shared      <- intersect(colnames(r1$counts), colnames(r2$counts))
  round1_only <- setdiff(colnames(r1$counts), colnames(r2$counts))
  round2_only <- setdiff(colnames(r2$counts), colnames(r1$counts))

  # ---- Build combined gene table --------------------------------------------
  # Rows for round-2-only genes are filled in COLUMN BY NAME, never by whole-
  # row positional assignment: the two rounds' gene tables are not guaranteed
  # to carry the same columns in the same order, and `df[rows, ] <- other_df`
  # matches by position, so a differing column order silently transposes
  # values between columns and a differing column count silently recycles.
  rows1 <- match(union_genes, genes1_key)
  rows2 <- match(union_genes, genes2_key)
  genes_out <- r1$genes[rows1, , drop = FALSE]
  need_r2 <- is.na(rows1) & !is.na(rows2)
  if (any(need_r2)) {
    for (gcol in colnames(genes_out)) {
      if (gcol %in% colnames(r2$genes)) {
        src <- r2$genes[[gcol]][rows2[need_r2]]
        # Factors from the two rounds can carry different level sets; compare
        # and store as character rather than emitting "invalid factor level".
        if (is.factor(genes_out[[gcol]]) || is.factor(src)) {
          genes_out[[gcol]] <- as.character(genes_out[[gcol]])
          src               <- as.character(src)
        }
        genes_out[[gcol]][need_r2] <- src
      } else {
        # Column exists only in round 1: round-2-only genes have no value.
        genes_out[[gcol]][need_r2] <- NA
      }
    }
    dropped <- setdiff(colnames(r2$genes), colnames(genes_out))
    if (length(dropped)) {
      warning(sample_name, ": round 2's gene table has column(s) absent from ",
              "round 1's (", paste(dropped, collapse = ", "), "); they are ",
              "not carried into the combined gene table.")
    }
  }
  genes_out[[gene_key_col]] <- union_genes
  rownames(genes_out) <- NULL

  # ---- Build combined cell table --------------------------------------------
  # Count-derived columns can't be trusted from either round individually
  # once counts have been summed, so recompute the two we can define
  # unambiguously straight from the combined matrix. All other columns are
  # passed through: round 1's value where round 1 has that barcode, else
  # round 2's; on a shared barcode where both rounds have the column and
  # disagree, round 1 wins (warn or error per `metadata_conflict`).
  # `nCount`/`nFeature` are matched as PREFIXES, not exact names: Seurat writes
  # these per-assay ("nCount_RNA", "nFeature_RNA"), so an exact-name test never
  # fires and a stale round-1 nCount_RNA would survive into the combined table
  # alongside a recomputed tscp_count that disagrees with it.
  count_derived_exact <- c("tscp_count", "gene_count", "umi_count", "reads",
                           "n_genes", "n_counts")
  all_cols   <- union(colnames(r1$cells), colnames(r2$cells))
  is_derived <- all_cols %in% count_derived_exact |
    grepl("^(nCount|nFeature)([._]|$)", all_cols)
  other_cols <- setdiff(all_cols[!is_derived], cell_id_col)

  idx1 <- match(union_cells, r1$cells[[cell_id_col]])
  idx2 <- match(union_cells, r2$cells[[cell_id_col]])

  cells_out <- data.frame(union_cells, stringsAsFactors = FALSE)
  names(cells_out) <- cell_id_col
  n_out <- length(union_cells)
  for (col in other_cols) {
    in1 <- col %in% colnames(r1$cells)
    in2 <- col %in% colnames(r2$cells)
    v1  <- if (in1) r1$cells[[col]][idx1] else rep(NA, n_out)
    v2  <- if (in2) r2$cells[[col]][idx2] else rep(NA, n_out)

    # A round can supply a value only where it has BOTH the column and the
    # barcode. Testing the barcode alone (`!is.na(idx1)`) drops round 2's
    # value on every shared barcode for any column round 1 doesn't have --
    # silent data loss for columns that exist in only one round.
    have1 <- in1 & !is.na(idx1)
    have2 <- in2 & !is.na(idx2)

    # Factor columns with different level sets make `v1 != v2` throw
    # ("level sets of factors are different"), so compare as character. The
    # same coercion is applied to the output when either side is a factor,
    # since combining two different level sets isn't well defined anyway.
    if (is.factor(v1) || is.factor(v2)) {
      v1 <- as.character(v1)
      v2 <- as.character(v2)
    }

    mismatch <- have1 & have2 & !is.na(v1) & !is.na(v2) &
      (as.character(v1) != as.character(v2))
    if (any(mismatch)) {
      msg <- sprintf(
        "%s: column '%s' disagrees between rounds for %d shared barcode(s); using round 1's value.",
        sample_name, col, sum(mismatch)
      )
      if (metadata_conflict == "error") stop(msg) else warning(msg)
    }

    # Round 2 as the base, round 1 overwriting wherever it actually has a
    # value -- this keeps "round 1 wins on conflict" while still falling back
    # to round 2 whenever round 1 has nothing to offer. Assigning into a
    # vector (rather than ifelse()) also preserves Date/POSIXct classes,
    # which ifelse() strips.
    out_v <- v2
    out_v[have1] <- v1[have1]
    cells_out[[col]] <- out_v
  }
  cells_out$gene_count <- Matrix::colSums(combined_counts > 0)
  cells_out$tscp_count <- Matrix::colSums(combined_counts)

  list(
    counts = combined_counts,
    genes  = genes_out,
    cells  = cells_out,
    summary = data.frame(
      sample_name          = sample_name,
      n_cells_round1       = ncol(r1$counts),
      n_cells_round2       = ncol(r2$counts),
      n_cells_shared       = length(shared),
      n_cells_round1_only  = length(round1_only),
      n_cells_round2_only  = length(round2_only),
      n_cells_combined     = ncol(combined_counts),
      n_genes_round1       = nrow(r1$counts),
      n_genes_round2       = nrow(r2$counts),
      n_genes_combined     = nrow(combined_counts),
      stringsAsFactors     = FALSE
    )
  )
}
