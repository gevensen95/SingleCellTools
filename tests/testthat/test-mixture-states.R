# Tests for call_mixture_states() and call_stress_states() -- Gaussian
# mixture model (via mclust) state calling. Both accept a Seurat object or a
# plain data.frame; most tests here use a data.frame directly since fitting
# only needs numeric metadata columns, not a real Seurat object. mclust is a
# hard Imports dependency (not Suggests), so these run unconditionally except
# where a Seurat fixture is used.

# Two well-separated 1D clusters (mean 0 / mean 10, sd 1) so BIC reliably
# selects G = 2 and state assignment is unambiguous -- not meant to represent
# realistic biology, just a deterministic mixture-fitting target.
.make_bimodal_df <- function(seed = 1, n_per_state = 30, celltype = "TypeA") {
  set.seed(seed)
  low  <- stats::rnorm(n_per_state, mean = 0,  sd = 1)
  high <- stats::rnorm(n_per_state, mean = 10, sd = 1)
  ids  <- paste0("cell", seq_len(2 * n_per_state))
  data.frame(
    score                  = c(low, high),
    celltype               = celltype,
    annotation_first_pass  = celltype,
    row.names              = ids,
    stringsAsFactors       = FALSE
  )
}

# ============================================================================
# call_mixture_states()
# ============================================================================

test_that("call_mixture_states recovers 2 well-separated states and ranks them correctly", {
  df  <- .make_bimodal_df(seed = 1)
  out <- call_mixture_states(df, score_col = "score", label = "test")

  expect_s3_class(out, "data.frame")
  expect_true(all(c("id", "test_state", "test_confidence", "test_prob_high") %in% colnames(out)))
  expect_equal(attr(out, "G"), 2)

  low_ids  <- rownames(df)[df$score < 5]
  high_ids <- rownames(df)[df$score >= 5]
  expect_true(all(out$test_state[out$id %in% low_ids]  == 1))
  expect_true(all(out$test_state[out$id %in% high_ids] == 2))
})

test_that("call_mixture_states with decreasing = TRUE flips the ranking", {
  df  <- .make_bimodal_df(seed = 2)
  out <- call_mixture_states(df, score_col = "score", decreasing = TRUE, label = "test")

  low_ids  <- rownames(df)[df$score < 5]
  high_ids <- rownames(df)[df$score >= 5]
  expect_true(all(out$test_state[out$id %in% high_ids] == 1))
  expect_true(all(out$test_state[out$id %in% low_ids]  == 2))
})

test_that("call_mixture_states stratifies by group_col/group_value", {
  df1 <- .make_bimodal_df(seed = 3, celltype = "TypeA")
  df2 <- .make_bimodal_df(seed = 4, celltype = "TypeB")
  rownames(df2) <- paste0("b_", rownames(df2))
  df <- rbind(df1, df2)

  out <- call_mixture_states(df, score_col = "score", group_col = "celltype",
                             group_value = "TypeA", label = "test")
  expect_true(all(out$id %in% rownames(df1)))
  expect_false(any(out$id %in% rownames(df2)))
})

test_that("call_mixture_states fits a joint multivariate mixture across several score columns", {
  df <- .make_bimodal_df(seed = 5)
  df$score2 <- df$score * 2 + stats::rnorm(nrow(df), sd = 0.1)

  out <- call_mixture_states(df, score_col = c("score", "score2"), label = "combo")
  expect_true(all(c("combo_state", "combo_confidence", "combo_prob_high") %in% colnames(out)))
  expect_equal(nrow(out), nrow(df))
})

test_that("call_mixture_states returns NULL (with a message) below min_n", {
  df <- .make_bimodal_df(seed = 6, n_per_state = 2)
  expect_message(
    out <- call_mixture_states(df, score_col = "score", min_n = 10),
    "Skipping"
  )
  expect_null(out)
})

test_that("call_mixture_states errors when score_col is missing from the data", {
  df <- .make_bimodal_df(seed = 7)
  expect_error(call_mixture_states(df, score_col = "nope"), "nope")
})

test_that("call_mixture_states errors when group_col is set without group_value", {
  df <- .make_bimodal_df(seed = 8)
  expect_error(
    call_mixture_states(df, score_col = "score", group_col = "celltype"),
    "group_value"
  )
})

test_that("call_mixture_states errors when group_col itself isn't found", {
  df <- .make_bimodal_df(seed = 8)
  expect_error(
    call_mixture_states(df, score_col = "score", group_col = "nope", group_value = "x"),
    "nope"
  )
})

test_that("call_mixture_states errors on input that's neither a Seurat object nor a data.frame", {
  expect_error(call_mixture_states(list(1, 2), score_col = "score"), "data.frame")
})

test_that("call_mixture_states works directly on a Seurat object's metadata", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 80)
  set.seed(1)
  obj$score <- c(stats::rnorm(40, 0, 1), stats::rnorm(40, 10, 1))

  out <- call_mixture_states(obj, score_col = "score", label = "test")
  expect_setequal(out$id, colnames(obj))
})


# ============================================================================
# call_stress_states() -- thin wrapper fixing group_col = "annotation_first_pass"
# and score_col = "stress_composite" by default
# ============================================================================

test_that("call_stress_states wraps call_mixture_states with the fixed column names", {
  df <- .make_bimodal_df(seed = 9, celltype = "Hepatocytes")
  colnames(df)[colnames(df) == "score"] <- "stress_composite"

  res <- call_stress_states(df, cell_type = "Hepatocytes")
  expect_named(res, c("calls", "G", "bic"))
  expect_true(all(c("cell", "stress_state", "stress_prob", "prob_stressed") %in% colnames(res$calls)))
  expect_equal(res$G, 2)
  expect_equal(nrow(res$calls), nrow(df))
})

test_that("call_stress_states returns NULL when the underlying fit is skipped (too few cells)", {
  df <- .make_bimodal_df(seed = 10, n_per_state = 2, celltype = "Hepatocytes")
  colnames(df)[colnames(df) == "score"] <- "stress_composite"

  expect_message(
    res <- call_stress_states(df, cell_type = "Hepatocytes", min_n = 10),
    "Skipping"
  )
  expect_null(res)
})

test_that("call_stress_states only fits within the requested cell_type", {
  df1 <- .make_bimodal_df(seed = 11, celltype = "Hepatocytes")
  df2 <- .make_bimodal_df(seed = 12, celltype = "Stellate")
  rownames(df2) <- paste0("s_", rownames(df2))
  df <- rbind(df1, df2)
  colnames(df)[colnames(df) == "score"] <- "stress_composite"

  res <- call_stress_states(df, cell_type = "Hepatocytes")
  expect_true(all(res$calls$cell %in% rownames(df1)))
})
