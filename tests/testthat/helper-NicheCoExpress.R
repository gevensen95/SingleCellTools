# Shared fixtures for NicheCoExpress() / plotNicheCoExpress() tests.
#
# Layout: 2 niches x 4 samples x `n_per_group` cells, 2 conditions (2 samples
# each: S1/S2 = "healthy", S3/S4 = "tumor"), optionally 2 cell types.
# Genes G1/G2 share a latent factor (co-expressed); G3-G6 are independent
# filler genes that serve as the abundance-matched background pool.

.skip_if_no_seurat <- function() {
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")
}

.make_coexpr_object <- function(seed = 1, with_celltypes = TRUE,
                                n_per_group = 40, cond_shift = 0) {
  set.seed(seed)
  samples     <- c("S1", "S2", "S3", "S4")
  sample_cond <- c(S1 = "healthy", S2 = "healthy", S3 = "tumor", S4 = "tumor")
  niches      <- c("N1", "N2")
  celltypes   <- c("Tcell", "Myeloid")

  # Build the full cell-id vector and expression matrix up front (indexed by
  # position), rather than rbind/cbind-ing a named list of per-group pieces
  # -- that pattern is fragile: do.call(rbind, ...) on a *named* list of
  # data.frames can rewrite row names when reconciling the list's own tags
  # against each piece's row.names, which then silently desyncs from the
  # colnames of the (separately built) expression matrix. Filling
  # pre-allocated structures by explicit column index avoids that entirely.
  groups   <- expand.grid(niche = niches, sample = samples,
                          stringsAsFactors = FALSE)
  n_groups <- nrow(groups)
  n_total  <- n_groups * n_per_group
  all_ids  <- paste0("cell", seq_len(n_total))

  niche_v     <- character(n_total)
  sample_v    <- character(n_total)
  condition_v <- character(n_total)
  celltype_v  <- character(n_total)
  expr <- matrix(NA_real_, nrow = 6, ncol = n_total,
                dimnames = list(paste0("G", 1:6), all_ids))

  for (i in seq_len(n_groups)) {
    nz  <- groups$niche[i]
    sp  <- groups$sample[i]
    cols <- ((i - 1) * n_per_group + 1):(i * n_per_group)

    niche_v[cols]    <- nz
    sample_v[cols]   <- sp
    condition_v[cols] <- sample_cond[[sp]]
    celltype_v[cols] <- if (with_celltypes) {
      sample(celltypes, n_per_group, replace = TRUE, prob = c(0.5, 0.5))
    } else {
      NA_character_
    }

    shift <- if (sample_cond[[sp]] == "tumor") cond_shift else 0
    base <- rgamma(n_per_group, shape = 3, rate = 1) + shift
    g1 <- pmax(base + rnorm(n_per_group, sd = 0.15), 0.01)
    g2 <- pmax(base + rnorm(n_per_group, sd = 0.15), 0.01)      # co-expressed w/ g1
    g3 <- pmax(rgamma(n_per_group, shape = 3, rate = 1), 0.01)   # independent
    g4 <- pmax(rgamma(n_per_group, shape = 1, rate = 1), 0.01)
    g5 <- pmax(rgamma(n_per_group, shape = 1, rate = 1), 0.01)
    g6 <- pmax(rgamma(n_per_group, shape = 1, rate = 1), 0.01)

    expr[, cols] <- rbind(g1, g2, g3, g4, g5, g6)
  }

  meta <- data.frame(
    niche            = niche_v,
    sample           = sample_v,
    condition        = condition_v,
    celltype         = celltype_v,
    row.names        = all_ids,
    stringsAsFactors = FALSE
  )
  if (!with_celltypes) meta$celltype <- NULL

  counts <- round(expr * 10)
  storage.mode(counts) <- "double"

  obj <- SeuratObject::CreateSeuratObject(counts = methods::as(counts, "CsparseMatrix"), meta.data = meta)
  SeuratObject::LayerData(obj, assay = "RNA", layer = "data") <- expr
  obj
}
