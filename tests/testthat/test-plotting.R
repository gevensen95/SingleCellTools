# Tests for PlotVolcano() (pure data.frame, no Seurat needed), MarkerHeatmap(),
# MarkerPlot(), MarkerPctPlot(), StackedViolinPlot(), PlotFeatureDensity(),
# and AddGenePositivity()/PlotGenePositivity() (tested together since the
# latter consumes exactly the metadata columns the former produces).

# ============================================================================
# PlotVolcano() -- pure data.frame, no Seurat object needed
# ============================================================================

.make_volcano_df <- function() {
  data.frame(
    gene       = paste0("G", 1:30),
    avg_log2FC = c(seq(3, 0.2, length.out = 10), seq(-0.2, -3, length.out = 10),
                  runif(10, -0.5, 0.5)),
    p_val_adj  = c(rep(0.001, 10), rep(0.001, 10), runif(10, 0.2, 0.9)),
    stringsAsFactors = FALSE
  )
}

test_that("PlotVolcano auto-detects columns and returns a ggplot", {
  df <- .make_volcano_df()
  p <- PlotVolcano(df)
  expect_s3_class(p, "ggplot")
})

test_that("PlotVolcano works with PseudobulkDE-style column names too", {
  df <- .make_volcano_df()
  names(df)[names(df) == "avg_log2FC"] <- "log2FC"
  names(df)[names(df) == "p_val_adj"]  <- "padj"
  p <- PlotVolcano(df)
  expect_s3_class(p, "ggplot")
})

test_that("PlotVolcano falls back to rownames when no gene column exists", {
  df <- .make_volcano_df()
  rownames(df) <- df$gene
  df$gene <- NULL
  p <- PlotVolcano(df)
  expect_s3_class(p, "ggplot")
})

test_that("PlotVolcano validates its inputs", {
  expect_error(PlotVolcano(list(1)), "data frame")
  expect_error(PlotVolcano(.make_volcano_df(), colors = c("a", "b")), "length-3")
  bad <- data.frame(gene = "G1", weird = 1)
  expect_error(PlotVolcano(bad), "auto-detect")
})

test_that("PlotVolcano top_n = 0 disables labeling without erroring", {
  df <- .make_volcano_df()
  p <- PlotVolcano(df, top_n = 0)
  expect_s3_class(p, "ggplot")
})

test_that("PlotVolcano auto-generates a title from contrast columns", {
  df <- .make_volcano_df()
  df$contrast_num <- "treated"
  df$contrast_den <- "control"
  p <- PlotVolcano(df)
  expect_true(grepl("treated vs control", p$labels$title))
})


# ============================================================================
# MarkerHeatmap()
# ============================================================================

test_that("MarkerHeatmap builds a heatmap from a supplied markers table", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60, n_clusters = 3)
  markers <- data.frame(
    gene       = rownames(obj)[1:10],
    cluster    = rep(levels(obj$seurat_clusters), length.out = 10),
    avg_log2FC = runif(10, 1, 3),
    p_val_adj  = rep(0.001, 10),
    stringsAsFactors = FALSE
  )
  p <- MarkerHeatmap(obj, markers = markers, n = 3)
  expect_s3_class(p, "ggplot")
})

test_that("MarkerHeatmap validates inputs", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(MarkerHeatmap(list(1)), "Seurat object")
  expect_error(
    MarkerHeatmap(obj, markers = data.frame(gene = "a")),
    "missing required columns"
  )
})

test_that("MarkerHeatmap errors when nothing passes the significance filter", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 10, n_cells = 20)
  markers <- data.frame(gene = rownames(obj)[1:5], cluster = "0",
                        avg_log2FC = 1, p_val_adj = 0.9,  # all non-significant
                        stringsAsFactors = FALSE)
  expect_error(MarkerHeatmap(obj, markers = markers), "No genes passed")
})

test_that("MarkerHeatmap pseudobulk = TRUE returns a ggplot with finite z-scores", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60, n_clusters = 3)
  markers <- data.frame(
    gene       = rownames(obj)[1:10],
    cluster    = rep(levels(obj$seurat_clusters), length.out = 10),
    avg_log2FC = runif(10, 1, 3),
    p_val_adj  = rep(0.001, 10),
    stringsAsFactors = FALSE
  )
  p <- suppressMessages(MarkerHeatmap(obj, markers = markers, n = 3, pseudobulk = TRUE))
  expect_s3_class(p, "ggplot")
  expect_true(all(is.finite(p$data$value)))
  expect_equal(p$labels$fill, "Z-score")
})

test_that("MarkerHeatmap pseudobulk = TRUE, scale_rows = FALSE labels the fill as pseudobulk log-CPM", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60, n_clusters = 3)
  markers <- data.frame(
    gene       = rownames(obj)[1:10],
    cluster    = rep(levels(obj$seurat_clusters), length.out = 10),
    avg_log2FC = runif(10, 1, 3),
    p_val_adj  = rep(0.001, 10),
    stringsAsFactors = FALSE
  )
  p <- suppressMessages(MarkerHeatmap(obj, markers = markers, n = 3,
                                      pseudobulk = TRUE, scale_rows = FALSE))
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$fill, "Pseudobulk log-CPM")
  expect_true(all(p$data$value >= 0))  # log1p(CPM) is never negative
})

test_that("MarkerHeatmap pseudobulk = TRUE messages and drops marker genes absent from the assay", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60, n_clusters = 3)
  markers <- data.frame(
    gene       = c(rownames(obj)[1:9], "NotAGene"),
    cluster    = rep(levels(obj$seurat_clusters), length.out = 10),
    avg_log2FC = runif(10, 1, 3),
    p_val_adj  = rep(0.001, 10),
    stringsAsFactors = FALSE
  )
  expect_message(
    p <- MarkerHeatmap(obj, markers = markers, n = 10, pseudobulk = TRUE),
    "NotAGene"
  )
  expect_false("NotAGene" %in% levels(p$data$gene))
})

test_that("MarkerHeatmap genes= plots a caller-supplied gene list, bypassing markers", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60, n_clusters = 3)
  want <- rownames(obj)[c(3, 1, 5)]
  p <- MarkerHeatmap(obj, genes = want)
  expect_s3_class(p, "ggplot")
  expect_setequal(levels(p$data$gene), want)
})

test_that("MarkerHeatmap errors when both markers and genes are supplied", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60, n_clusters = 3)
  markers <- data.frame(gene = rownames(obj)[1:3], cluster = "0",
                        avg_log2FC = 1, p_val_adj = 0.001,
                        stringsAsFactors = FALSE)
  expect_error(
    MarkerHeatmap(obj, markers = markers, genes = rownames(obj)[1:3]),
    "either .markers. or .genes."
  )
})

test_that("MarkerHeatmap genes= validates and drops genes absent from the assay", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60, n_clusters = 3)
  expect_error(MarkerHeatmap(obj, genes = character(0)), "non-empty character vector")
  expect_message(
    p <- MarkerHeatmap(obj, genes = c(rownames(obj)[1:3], "NotAGene")),
    "NotAGene"
  )
  expect_false("NotAGene" %in% levels(p$data$gene))
  expect_error(MarkerHeatmap(obj, genes = "NotAGene"), "None of the requested genes")
})


# ============================================================================
# MarkerPlot() / MarkerPctPlot()
# ============================================================================

.make_annotated_genes_df <- function(obj) {
  g <- rownames(obj)[1:6]
  data.frame(Gene = g, Details = rep(c("TypeA", "TypeB"), each = 3),
            stringsAsFactors = FALSE)
}

test_that("MarkerPlot returns a ggplot with annotation labels", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60, n_clusters = 3)
  genes_df <- .make_annotated_genes_df(obj)
  p <- suppressMessages(MarkerPlot(obj, genes_df))
  expect_s3_class(p, "ggplot")
})

test_that("MarkerPlot validates inputs", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(MarkerPlot(list(1), data.frame(a = 1, b = 2)), "Seurat object")
  expect_error(MarkerPlot(obj, "not_a_df"), "data frame")
  expect_error(MarkerPlot(obj, .make_annotated_genes_df(obj), assay = "nope"), "nope")
})

test_that("MarkerPlot drops genes not present in the assay and errors if none remain", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60)
  genes_df <- .make_annotated_genes_df(obj)
  genes_df$Gene[1] <- "NotAGene"
  expect_message(
    p <- MarkerPlot(obj, genes_df),
    "not in assay"
  )
  expect_s3_class(p, "ggplot")

  expect_error(
    suppressMessages(MarkerPlot(obj, data.frame(Gene = "NotAGene", Details = "x"))),
    "No genes remain"
  )
})

test_that("MarkerPlot auto-sizes axis text / annotation label font from gene count", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60, n_clusters = 3)
  genes_df <- .make_annotated_genes_df(obj)  # 6 genes -> the <=30 bucket

  p <- suppressMessages(MarkerPlot(obj, genes_df))
  expect_equal(p$theme$axis.text.y$size, 9)

  annot_layer <- p$layers[[length(p$layers)]]
  expect_equal(annot_layer$aes_params$size, 3.5)
})

test_that("MarkerPlot respects explicit label.fontsize/axis_text_size overrides", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60, n_clusters = 3)
  genes_df <- .make_annotated_genes_df(obj)

  p <- suppressMessages(MarkerPlot(obj, genes_df,
                                   axis_text_size = 11, label.fontsize = 6))
  expect_equal(p$theme$axis.text.y$size, 11)

  annot_layer <- p$layers[[length(p$layers)]]
  expect_equal(annot_layer$aes_params$size, 6)
})

test_that("MarkerPlot attaches suggested_width/suggested_height attributes", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60, n_clusters = 3)
  genes_df <- .make_annotated_genes_df(obj)

  p <- suppressMessages(MarkerPlot(obj, genes_df))
  expect_true(is.numeric(attr(p, "suggested_width")))
  expect_true(is.numeric(attr(p, "suggested_height")))
  expect_gt(attr(p, "suggested_width"), 0)
  expect_gt(attr(p, "suggested_height"), 0)

  # Explicit width/height override the auto-computed suggestion used for
  # save_path, but the attributes always reflect the auto-computed values
  # (they're a suggestion for your OWN ggsave() call, not an echo of what
  # you passed in).
  p2 <- suppressMessages(MarkerPlot(obj, genes_df, width = 12, height = 12))
  expect_equal(attr(p2, "suggested_width"), attr(p, "suggested_width"))
  expect_equal(attr(p2, "suggested_height"), attr(p, "suggested_height"))
})

test_that("MarkerPlot's save_path writes a file at the requested/auto-computed size", {
  .skip_if_missing("Seurat", "SeuratObject")
  testthat::skip_if_not_installed("ggplot2")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60, n_clusters = 3)
  genes_df <- .make_annotated_genes_df(obj)

  out <- tempfile(fileext = ".pdf")
  on.exit(unlink(out), add = TRUE)
  expect_message(
    p <- MarkerPlot(obj, genes_df, save_path = out, width = 5, height = 5),
    "Saving to"
  )
  expect_true(file.exists(out))
  expect_s3_class(p, "ggplot")
})

test_that("MarkerPlot messages a suggested ggsave() call when save_path isn't set", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60, n_clusters = 3)
  genes_df <- .make_annotated_genes_df(obj)

  expect_message(
    MarkerPlot(obj, genes_df),
    "if labels look cramped"
  )
})

# ============================================================================
# TopMarkerPlot() -- FindAllMarkers -> top N per cluster -> MarkerPlot()
# ============================================================================

# A fabricated FindAllMarkers()-shaped table over real genes of `obj`. Most
# tests below use this rather than calling FindAllMarkers(), which is slow and
# finds nothing on the uniform-Poisson fixture; the ranking / de-duplication /
# label logic is what TopMarkerPlot() actually owns, and it is fully
# exercisable from a supplied table. The FindAllMarkers path gets its own
# signal-bearing fixture at the end.
.fake_markers <- function(obj, n_clusters = 3, per = 5) {
  g <- rownames(obj)
  stopifnot(length(g) >= n_clusters * per)
  do.call(rbind, lapply(seq_len(n_clusters) - 1L, function(k) {
    data.frame(gene       = g[k * per + seq_len(per)],
               cluster    = as.character(k),
               avg_log2FC = seq(3, 1, length.out = per),
               p_val_adj  = 1e-5,
               stringsAsFactors = FALSE)
  }))
}

test_that("TopMarkerPlot returns a ggplot and attaches the marker tables", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 30, n_cells = 90, n_clusters = 3)
  mk  <- .fake_markers(obj, n_clusters = 3, per = 5)

  p <- suppressMessages(TopMarkerPlot(obj, n = 3, markers = mk))
  expect_s3_class(p, "ggplot")

  top <- attr(p, "top_markers")
  expect_equal(nrow(top), 9)                      # 3 clusters x top 3
  expect_setequal(unique(top$cluster), c("0", "1", "2"))
  expect_identical(attr(p, "markers"), mk)        # full table passed through
})

test_that("TopMarkerPlot assigns a gene topping two clusters to its best one", {
  # Duplicated genes would render as repeated, identical dot-plot rows.
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 30, n_cells = 90, n_clusters = 3)
  mk  <- .fake_markers(obj, n_clusters = 3, per = 5)

  shared <- mk$gene[1]                            # cluster 0's best gene
  mk <- rbind(mk, data.frame(gene = shared, cluster = "2",
                             avg_log2FC = 99, p_val_adj = 1e-9,
                             stringsAsFactors = FALSE))

  p   <- suppressMessages(TopMarkerPlot(obj, n = 5, markers = mk))
  top <- attr(p, "top_markers")
  expect_equal(sum(top$gene == shared), 1L)
  expect_equal(top$cluster[top$gene == shared], "2")  # the higher log2FC wins
})

test_that("TopMarkerPlot zero-pads numeric cluster labels past 9", {
  # Regression: MarkerPlot() orders its annotation groups ALPHABETICALLY
  # (order(Details) / table(Details)), so unpadded numeric IDs put cluster 10
  # between 1 and 2. Padding makes alphabetical order match numeric order.
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 80, n_cells = 240, n_clusters = 12)
  mk  <- .fake_markers(obj, n_clusters = 12, per = 5)

  p   <- suppressMessages(TopMarkerPlot(obj, n = 2, markers = mk))
  labs <- sort(unique(attr(p, "top_markers")$cluster))
  expect_identical(labs, formatC(0:11, width = 2, flag = "0"))
  expect_identical(labs, labs[order(labs)])       # alphabetical == numeric

  # Fewer than 10 clusters needs no padding -- labels stay as-is.
  obj3 <- .make_small_seurat(seed = 1, n_genes = 30, n_cells = 90, n_clusters = 3)
  p3   <- suppressMessages(TopMarkerPlot(obj3, n = 2,
                                         markers = .fake_markers(obj3, 3, 5)))
  expect_setequal(unique(attr(p3, "top_markers")$cluster), c("0", "1", "2"))
})

test_that("TopMarkerPlot validates inputs and reports an empty result clearly", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 30, n_cells = 90, n_clusters = 3)
  mk  <- .fake_markers(obj, n_clusters = 3, per = 5)

  expect_error(TopMarkerPlot(list(1)), "Seurat object")
  expect_error(TopMarkerPlot(obj, n = 0), "positive number")
  expect_error(TopMarkerPlot(obj, group_by = "nope"), "not found in metadata")
  expect_error(TopMarkerPlot(obj, markers = data.frame(a = 1)),
               "missing column")
  expect_error(suppressMessages(TopMarkerPlot(obj, markers = mk, max_padj = 0)),
               "No markers passed")
})

test_that("TopMarkerPlot group_by relabels clusters and restores Idents", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 30, n_cells = 90,
                            n_clusters = 3, n_celltypes = 3)
  before <- Seurat::Idents(obj)

  # Label the fabricated markers with this fixture's actual celltype values
  # so the annotation groups correspond to the identities being plotted.
  ct <- sort(unique(as.character(obj$celltype)))
  mk <- .fake_markers(obj, n_clusters = length(ct), per = 5)
  mk$cluster <- ct[as.integer(mk$cluster) + 1L]

  p <- suppressMessages(TopMarkerPlot(obj, n = 2, markers = mk,
                                      group_by = "celltype"))
  expect_setequal(unique(attr(p, "top_markers")$cluster), ct)
  expect_identical(Seurat::Idents(obj), before)   # caller's object untouched
})

test_that("TopMarkerPlot runs FindAllMarkers end-to-end when markers = NULL", {
  # The only test that exercises the FindAllMarkers path, so it needs a
  # fixture with real per-cluster signal -- uniform Poisson counts yield no
  # markers and the call would find nothing to plot.
  .skip_if_missing("Seurat", "SeuratObject")
  testthat::skip_on_cran()
  set.seed(1)
  ng <- 60; nc <- 180; ngrp <- 3
  grp <- rep(seq_len(ngrp), length.out = nc)
  lam <- matrix(0.6, ng, nc)
  for (g in seq_len(ngrp)) lam[((g - 1) * 8 + 1):(g * 8), grp == g] <- 8
  lam <- sweep(lam, 2, stats::rlnorm(nc, 0, 0.3), "*")
  m <- matrix(stats::rpois(ng * nc, lam), ng,
              dimnames = list(paste0("Gene", seq_len(ng)),
                              paste0("c", seq_len(nc))))
  storage.mode(m) <- "double"
  obj <- suppressWarnings(
    SeuratObject::CreateSeuratObject(counts = methods::as(m, "CsparseMatrix")))
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  Seurat::Idents(obj) <- factor(grp - 1L)

  p <- suppressWarnings(suppressMessages(TopMarkerPlot(obj, n = 3)))
  expect_s3_class(p, "ggplot")
  top <- attr(p, "top_markers")
  expect_gt(nrow(top), 0)
  expect_false(anyDuplicated(top$gene) > 0)
  expect_true(all(top$p_val_adj < 0.05))
})


test_that("MarkerPctPlot returns a ggplot in both tile and dot styles", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60, n_clusters = 3)
  genes_df <- .make_annotated_genes_df(obj)
  p1 <- suppressMessages(MarkerPctPlot(obj, genes_df, style = "tile"))
  p2 <- suppressMessages(MarkerPctPlot(obj, genes_df, style = "dot"))
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("MarkerPctPlot validates its pct-range and color arguments", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  genes_df <- .make_annotated_genes_df(obj)
  expect_error(MarkerPctPlot(obj, genes_df, min_pct = 50, max_pct = 10), "min_pct")
  expect_error(MarkerPctPlot(obj, genes_df, colors = "onlyone"), "length-2")
})


# ============================================================================
# StackedViolinPlot()
# ============================================================================

test_that("StackedViolinPlot returns a ggplot faceted by gene", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 15, n_cells = 40, n_clusters = 3)
  p <- StackedViolinPlot(obj, features = rownames(obj)[1:5])
  expect_s3_class(p, "ggplot")
})

test_that("StackedViolinPlot errors when no requested feature is in the assay", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(StackedViolinPlot(obj, features = "NotAGene"), "None of the requested")
})

test_that("StackedViolinPlot group.by validates the metadata column", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 10, n_cells = 20)
  expect_error(
    StackedViolinPlot(obj, features = rownames(obj)[1:3], group.by = "nope"),
    "nope"
  )
})


# ============================================================================
# PlotFeatureDensity()
# ============================================================================

.add_fake_umap <- function(obj) {
  emb <- matrix(rnorm(ncol(obj) * 2), nrow = ncol(obj),
               dimnames = list(colnames(obj), c("UMAP_1", "UMAP_2")))
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "UMAP_",
                                                      assay = "RNA")
  obj
}

test_that("PlotFeatureDensity returns a ggplot for gene features", {
  .skip_if_missing("Seurat", "SeuratObject", "ks")
  obj <- .make_small_seurat(seed = 1, n_genes = 10, n_cells = 60)
  obj <- .add_fake_umap(obj)
  p <- PlotFeatureDensity(obj, features = rownames(obj)[1:2])
  expect_s3_class(p, "ggplot")
})

test_that("PlotFeatureDensity works with a metadata-column feature and joint = TRUE", {
  .skip_if_missing("Seurat", "SeuratObject", "ks")
  obj <- .make_small_seurat(seed = 1, n_genes = 10, n_cells = 60)
  obj <- .add_fake_umap(obj)
  p <- PlotFeatureDensity(obj, features = c("nCount_RNA", "nFeature_RNA"), joint = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("PlotFeatureDensity validates the reduction and feature name", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(PlotFeatureDensity(list(1), "G1"), "Seurat object")
  expect_error(PlotFeatureDensity(obj, "G1", reduction = "umap"), "not found")
  obj <- .add_fake_umap(obj)
  expect_error(PlotFeatureDensity(obj, "NotAFeature"), "not found")
})


# ============================================================================
# AddGenePositivity() + PlotGenePositivity()
# ============================================================================

test_that("AddGenePositivity adds logical <gene>_pos columns", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 10, n_cells = 40)
  genes <- rownames(obj)[1:3]
  out <- AddGenePositivity(obj, genes = genes, layer = "counts")
  expect_true(all(paste0(genes, "_pos") %in% colnames(out@meta.data)))
  expect_true(is.logical(out@meta.data[[paste0(genes[1], "_pos")]]))
})

test_that("AddGenePositivity drops genes missing from any object in a list, with a warning", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj1 <- .make_small_seurat(seed = 1, n_genes = 10, n_cells = 20, gene_prefix = "GeneA")
  obj2 <- .make_small_seurat(seed = 2, n_genes = 10, n_cells = 20, gene_prefix = "GeneB")
  expect_warning(
    out <- AddGenePositivity(list(a = obj1, b = obj2),
                             genes = c("GeneA1", "GeneB1")),
    "not present in every object"
  )
  # Neither gene is shared across both objects, so nothing gets added and a
  # second warning fires; both objects are returned unchanged.
  expect_type(out, "list")
})

test_that("PlotGenePositivity renders bar/heatmap/combo styles from real AddGenePositivity output", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 10, n_cells = 60, n_clusters = 3)
  genes <- rownames(obj)[1:2]
  obj <- AddGenePositivity(obj, genes = genes, layer = "counts")

  p_bar <- PlotGenePositivity(obj, genes = genes, group.by = "seurat_clusters", style = "bar")
  p_hm  <- PlotGenePositivity(obj, genes = genes, group.by = "seurat_clusters", style = "heatmap")
  p_combo <- PlotGenePositivity(obj, genes = genes, group.by = "seurat_clusters", style = "combo")
  expect_s3_class(p_bar, "ggplot")
  expect_s3_class(p_hm, "ggplot")
  expect_s3_class(p_combo, "ggplot")
})

test_that("PlotGenePositivity errors when no positivity data is available", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(
    suppressWarnings(PlotGenePositivity(obj, genes = "NotAGene", drop_absent = TRUE)),
    "No positivity data"
  )
})

test_that(".collect_positivity attributes prefix-colliding gene names correctly (e.g. Gene1 vs Gene10)", {
  # Regression test: .one()'s gene-name recovery used to reverse-engineer
  # the gene from the column name via startsWith() prefix-matching against
  # `genes`, which picks the FIRST matching gene in `genes` order -- so with
  # genes = c("Gene1", "Gene10"), the column "Gene10_pos" incorrectly
  # matched "Gene1" first (since "Gene10_pos" also startsWith("Gene1")) and
  # got mislabeled as gene "Gene1" instead of "Gene10". The fix threads the
  # actual `suffix` parameter through directly instead of re-deriving it.
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 15, n_cells = 40, n_clusters = 2)
  genes <- c("Gene1", "Gene10")
  obj <- AddGenePositivity(obj, genes = genes, layer = "counts")

  df <- SingleCellTools:::.collect_positivity(
    obj, pos_cols = paste0(genes, "_pos"), genes = genes,
    group.by = "seurat_clusters", sample_col = NULL, drop_absent = TRUE,
    suffix = "_pos"
  )
  expect_setequal(unique(df$gene), genes)
  # Every row's gene should trace back to a column that actually contains
  # its own name as the exact prefix before "_pos" -- not silently
  # collapsed onto "Gene1" for both.
  expect_equal(sum(df$gene == "Gene10"), sum(df$gene == "Gene1"))
})
