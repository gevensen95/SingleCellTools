#!/usr/bin/env Rscript
# ============================================================================
# build_docs.R — render every vignette into doc/ as .Rmd, .md and .html
#
# LAYOUT
#   vignettes/   authoritative SOURCE (.Rmd). R CMD build only looks here, so
#                this is where sources must live for vignette() and
#                browseVignettes() to work and for the tarball to ship them.
#   doc/         GENERATED output, three formats per vignette. Tracked in git
#                so the rendered docs are browsable on GitHub, but listed in
#                .Rbuildignore so the tarball doesn't ship them twice.
#
# Never hand-edit anything in doc/ — this script overwrites it.
#
# Usage:  Rscript dev/build_docs.R          # from the package root
#         source("dev/build_docs.R"); build_docs()
# ============================================================================

build_docs <- function(vig_dir  = "vignettes",
                       out_dir  = "doc",
                       formats  = c("rmd", "md", "html"),
                       vignettes = NULL,
                       quiet    = TRUE) {

  for (p in c("rmarkdown", "knitr")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Package '", p, "' is required. install.packages('", p, "')")
    }
  }
  # pandoc discovery. rmarkdown searches a fixed list of locations plus
  # RSTUDIO_PANDOC -- that list does NOT necessarily include wherever your
  # package manager put it (Homebrew on Apple Silicon uses /opt/homebrew/bin,
  # which rmarkdown does not look in). So if rmarkdown can't see pandoc but
  # it is on PATH, point rmarkdown at it rather than giving up.
  if (!rmarkdown::pandoc_available()) {
    found <- Sys.which("pandoc")
    if (nzchar(found)) {
      Sys.setenv(RSTUDIO_PANDOC = dirname(found))
      try(rmarkdown::find_pandoc(cache = FALSE), silent = TRUE)
      message("Found pandoc on PATH at ", found,
              " and pointed rmarkdown at it.")
    } else {
      # Deliberately not a hard stop: render() reports a missing pandoc
      # perfectly well on its own, and a pre-check that is stricter than the
      # thing it guards just invents a failure mode.
      message("Note: rmarkdown cannot find pandoc, and it is not on PATH.\n",
              "  Sys.which('pandoc'):  ", "'", Sys.which("pandoc"), "'\n",
              "  RSTUDIO_PANDOC:       ",
              Sys.getenv("RSTUDIO_PANDOC", unset = "<unset>"), "\n",
              "  rmarkdown::pandoc_version(): ",
              tryCatch(as.character(rmarkdown::pandoc_version()),
                       error = function(e) "<none>"), "\n",
              "If pandoc IS installed, tell R where:\n",
              "  Sys.setenv(RSTUDIO_PANDOC = \"/opt/homebrew/bin\")   # or your path\n",
              "or add RSTUDIO_PANDOC=/opt/homebrew/bin to ~/.Renviron.\n",
              "Continuing -- render() will report the real error if it is missing.")
    }
  }
  if (!dir.exists(vig_dir)) {
    stop("No '", vig_dir, "' directory here. Run from the package root.")
  }

  rmds <- list.files(vig_dir, pattern = "\\.Rmd$", full.names = TRUE)
  if (!is.null(vignettes)) {
    keep <- tools::file_path_sans_ext(basename(rmds)) %in%
      tools::file_path_sans_ext(vignettes)
    if (!any(keep)) {
      stop("None of the requested vignettes were found in '", vig_dir, "'.")
    }
    rmds <- rmds[keep]
  }
  if (!length(rmds)) stop("No .Rmd files found in '", vig_dir, "'.")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  results <- lapply(rmds, function(f) {
    base <- tools::file_path_sans_ext(basename(f))
    message("--- ", base, " ---")
    made <- character(0)

    # .Rmd: a verbatim copy of the source, so doc/ is self-contained.
    if ("rmd" %in% formats) {
      file.copy(f, file.path(out_dir, basename(f)), overwrite = TRUE)
      made <- c(made, "Rmd")
    }

    # .md: github_document, so GitHub renders it as a page rather than
    # showing raw source the way it does for .Rmd. html_preview = FALSE
    # stops it also emitting a throwaway preview .html next to the source.
    if ("md" %in% formats) {
      rmarkdown::render(
        f,
        output_format = rmarkdown::github_document(html_preview = FALSE),
        output_file   = paste0(base, ".md"),
        output_dir    = out_dir,
        quiet         = quiet
      )
      made <- c(made, "md")
    }

    # .html: the same html_vignette styling the installed package uses.
    if ("html" %in% formats) {
      rmarkdown::render(
        f,
        output_format = rmarkdown::html_vignette(),
        output_file   = paste0(base, ".html"),
        output_dir    = out_dir,
        quiet         = quiet
      )
      made <- c(made, "html")
    }

    data.frame(vignette = base, formats = paste(made, collapse = ", "),
               stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, results)
  message("\n--- Wrote ", nrow(out), " vignette(s) to ", out_dir, "/ ---")
  print(out, row.names = FALSE)

  # github_document leaves a *_files/ directory behind when a document
  # produces figures. Everything here is eval = FALSE so there should be
  # none; clean up any that do appear rather than committing them.
  strays <- list.files(out_dir, pattern = "_files$", full.names = TRUE)
  strays <- strays[dir.exists(strays)]
  if (length(strays)) {
    unlink(strays, recursive = TRUE)
    message("Removed ", length(strays), " empty *_files/ directory(ies).")
  }

  invisible(out)
}

# ---------------------------------------------------------------------------
# Also refresh Meta/vignette.rds and the tangled .R files, so
# vignette("...") resolves against a source checkout after load_all().
# Skipped if devtools isn't installed -- the three formats above don't
# depend on it.
# ---------------------------------------------------------------------------
build_docs_full <- function(...) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    message("--- devtools::build_vignettes() ---")
    devtools::build_vignettes(dependencies = FALSE, upgrade = "never")
  } else {
    message("devtools not installed; skipping build_vignettes(). ",
            "doc/ will have .Rmd/.md/.html but no .R tangles or Meta/.")
  }
  build_docs(...)
}

if (!interactive() && identical(environment(), globalenv())) {
  build_docs()
}
