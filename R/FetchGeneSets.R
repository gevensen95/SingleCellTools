#' Fetch a named gene-set list (GO, KEGG, Reactome, MSigDB collections, ...)
#'
#' Companion to \code{\link{RunPathwayEnrichment}} and
#' \code{\link{RunSingleCellGSEA}}: both take a \code{gene_sets} argument --
#' a named list of gene-set name -> gene symbols -- and this builds that
#' list from a named collection instead of you having to hand-assemble it
#' from \code{msigdbr} yourself each time.
#'
#' \strong{Collections} (case-insensitive, pass one or more):
#' \describe{
#'   \item{\code{"hallmark"}}{MSigDB Hallmark (H) -- 50 broad, non-redundant
#'     processes.}
#'   \item{\code{"go_bp"}, \code{"go_mf"}, \code{"go_cc"}}{Gene Ontology
#'     Biological Process / Molecular Function / Cellular Component.}
#'   \item{\code{"kegg"}}{KEGG pathways. Via \code{source = "msigdbr"} this
#'     is MSigDB's \code{CP:KEGG_MEDICUS} collection -- \strong{not} the
#'     original KEGG gene sets, which MSigDB removed for licensing reasons
#'     starting with release 2023.1 and replaced with this re-derived
#'     version. If you specifically need the original KEGG pathways, fetch
#'     them yourself via \code{KEGGREST} instead.}
#'   \item{\code{"reactome"}, \code{"biocarta"}, \code{"wikipathways"}}{The
#'     matching MSigDB C2:CP sub-collections.}
#'   \item{\code{"c7_immunologic"}, \code{"c8_celltype"}}{MSigDB C7 / C8.}
#' }
#'
#' \strong{Source.} \code{"msigdbr"} (default) covers every collection
#' above via the \code{msigdbr} package, and supports every species
#' \code{msigdbr} does (see \code{msigdbr::msigdbr_species()} for the
#' authoritative list for your installed version -- roughly two dozen,
#' including the common model organisms). \code{"godb"} instead builds GO
#' gene sets (\code{go_bp}/\code{go_mf}/\code{go_cc} only) directly from
#' \code{GO.db} + \code{AnnotationDbi} + an \code{org.*.eg.db} annotation
#' package -- \code{GO.db}/\code{AnnotationDbi} are already Imports of this
#' package, so this avoids adding \code{msigdbr} as a dependency if all you
#' want is GO. \code{"godb"} supports: Homo sapiens, Mus musculus, Rattus
#' norvegicus, Danio rerio, Drosophila melanogaster, Caenorhabditis
#' elegans, Bos taurus, Sus scrofa, Gallus gallus, Canis lupus familiaris,
#' Macaca mulatta, Pan troglodytes -- each backed by its own
#' \code{org.*.eg.db} Suggests. Gene-set names for \code{"godb"} are just
#' \code{"<GO ID>_<term>"} rather than MSigDB's curated naming.
#'
#' @param collection One or more of the collection names listed above.
#' @param species Full species name, e.g. \code{"Homo sapiens"},
#'   \code{"Mus musculus"} -- \strong{required, no default}. For
#'   \code{source = "msigdbr"}, passed straight to
#'   \code{msigdbr::msigdbr(species = ...)}; for \code{source = "godb"},
#'   must be one of the species listed in Details. There is deliberately no
#'   default and no auto-detection: gene symbol casing alone can't reliably
#'   tell some species apart (e.g. mouse and rat both use title case), so
#'   guessing risks silently fetching the wrong species' gene sets.
#' @param source \code{"msigdbr"} (default) or \code{"godb"}. See Details.
#' @param min_size,max_size Optional post-fetch size filter on each
#'   resulting gene set (e.g. to trim a huge collection like
#'   \code{"c2_curated"}-scale collections down before scoring).
#'   \code{NULL} (default, both) skips this -- \code{\link{RunPathwayEnrichment}}/
#'   \code{\link{RunSingleCellGSEA}} apply their own size filtering anyway,
#'   after restricting to genes actually present in your data.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return A named list of character vectors (gene-set name -> gene
#'   symbols), ready to pass as \code{gene_sets} to
#'   \code{\link{RunPathwayEnrichment}} or \code{\link{RunSingleCellGSEA}}
#'   (both of which also accept \code{collection}/\code{species} directly
#'   and call this internally, so you don't have to call it yourself
#'   unless you want the list for something else too).
#' @examples
#' \dontrun{
#' go_bp <- FetchGeneSets("go_bp", species = "Homo sapiens")
#' kegg  <- FetchGeneSets("kegg", species = "Homo sapiens")
#' both  <- FetchGeneSets(c("go_bp", "reactome"), species = "Mus musculus")
#'
#' # GO only, no msigdbr dependency
#' go_bp2 <- FetchGeneSets("go_bp", species = "Rattus norvegicus", source = "godb")
#'
#' de <- PseudobulkDE(obj, sample_col = "orig.ident", condition_col = "treatment",
#'                    ident_1 = "drug", ident_2 = "vehicle")
#' RunPathwayEnrichment(de, gene_sets = FetchGeneSets("hallmark", species = "Homo sapiens"))
#' }
#' @export
FetchGeneSets <- function(collection,
                          species  = NULL,
                          source   = c("msigdbr", "godb"),
                          min_size = NULL,
                          max_size = NULL,
                          verbose  = TRUE) {

  source <- match.arg(source)
  if (is.null(species)) {
    stop("`species` is required (e.g. 'Homo sapiens', 'Mus musculus') -- ",
         "there is no default. Gene symbol casing alone can't reliably ",
         "distinguish some species (mouse vs. rat, for instance), so this ",
         "is never guessed for you.")
  }
  collection <- tolower(collection)
  valid <- names(.msigdbr_collection_map())
  bad <- setdiff(collection, valid)
  if (length(bad) > 0) {
    stop("Unknown collection(s): ", paste(bad, collapse = ", "), ". ",
         "Valid collections: ", paste(valid, collapse = ", "))
  }
  if (source == "godb") {
    bad_godb <- setdiff(collection, c("go_bp", "go_mf", "go_cc"))
    if (length(bad_godb) > 0) {
      stop("source = 'godb' only supports GO collections (go_bp/go_mf/go_cc); ",
           "got: ", paste(bad_godb, collapse = ", "), ". Use source = 'msigdbr' ",
           "for the others.")
    }
  }

  out <- list()
  for (col in collection) {
    if (isTRUE(verbose)) {
      message(sprintf("--- Fetching '%s' (source = %s, species = %s) ---",
                      col, source, species))
    }
    sets <- if (source == "msigdbr") {
      .fetch_msigdbr_one(.msigdbr_collection_map()[[col]], species)
    } else {
      .fetch_go_via_godb(toupper(sub("^go_", "", col)), species)
    }
    if (isTRUE(verbose)) message(sprintf("  %d gene set(s).", length(sets)))
    out <- c(out, sets)
  }

  if (!is.null(min_size) || !is.null(max_size)) {
    sizes <- vapply(out, length, integer(1))
    keep <- rep(TRUE, length(out))
    if (!is.null(min_size)) keep <- keep & sizes >= min_size
    if (!is.null(max_size)) keep <- keep & sizes <= max_size
    if (isTRUE(verbose) && sum(!keep) > 0) {
      message(sprintf("  Dropping %d/%d gene set(s) outside [%s, %s].",
                      sum(!keep), length(out),
                      if (is.null(min_size)) "-Inf" else min_size,
                      if (is.null(max_size)) "Inf" else max_size))
    }
    out <- out[keep]
  }

  out
}


# ============================================================================
# Internal helpers for FetchGeneSets() -- also used directly by
# RunPathwayEnrichment() / RunSingleCellGSEA() when their `collection`
# argument is supplied instead of `gene_sets`.
# ============================================================================

#' @keywords internal
#' @noRd
.msigdbr_collection_map <- function() {
  list(
    hallmark       = list(category = "H",  subcategory = NULL),
    go_bp          = list(category = "C5", subcategory = "GO:BP"),
    go_mf          = list(category = "C5", subcategory = "GO:MF"),
    go_cc          = list(category = "C5", subcategory = "GO:CC"),
    kegg           = list(category = "C2", subcategory = "CP:KEGG_MEDICUS"),
    reactome       = list(category = "C2", subcategory = "CP:REACTOME"),
    biocarta       = list(category = "C2", subcategory = "CP:BIOCARTA"),
    wikipathways   = list(category = "C2", subcategory = "CP:WIKIPATHWAYS"),
    c2_curated     = list(category = "C2", subcategory = NULL),
    c7_immunologic = list(category = "C7", subcategory = NULL),
    c8_celltype    = list(category = "C8", subcategory = NULL)
  )
}

#' @keywords internal
#' @noRd
.fetch_msigdbr_one <- function(cat_info, species) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("'msigdbr' is required for source = 'msigdbr'. Install with ",
         "install.packages('msigdbr').")
  }
  args <- list(species = species, category = cat_info$category)
  if (!is.null(cat_info$subcategory)) args$subcategory <- cat_info$subcategory
  df <- do.call(msigdbr::msigdbr, args)

  # msigdbr's exact column names have shifted across versions -- match
  # defensively rather than hard-coding one version's schema (same approach
  # RunRCTD.R takes for spacexr's weights_doublet columns).
  name_col <- intersect(c("gs_name", "gene_set_name"), colnames(df))[1]
  gene_col <- intersect(c("gene_symbol", "human_gene_symbol"), colnames(df))[1]
  if (is.na(name_col) || is.na(gene_col)) {
    stop("Could not find the expected gene-set-name / gene-symbol columns ",
         "in msigdbr's output -- this may indicate an incompatible msigdbr ",
         "version. Columns found: ", paste(colnames(df), collapse = ", "))
  }
  sets <- split(df[[gene_col]], df[[name_col]])
  lapply(sets, unique)
}

#' @keywords internal
#' @noRd
.godb_species_map <- function() {
  c(
    "Homo sapiens"           = "org.Hs.eg.db",
    "Mus musculus"           = "org.Mm.eg.db",
    "Rattus norvegicus"      = "org.Rn.eg.db",
    "Danio rerio"            = "org.Dr.eg.db",
    "Drosophila melanogaster" = "org.Dm.eg.db",
    "Caenorhabditis elegans" = "org.Ce.eg.db",
    "Bos taurus"             = "org.Bt.eg.db",
    "Sus scrofa"             = "org.Ss.eg.db",
    "Gallus gallus"          = "org.Gg.eg.db",
    "Canis lupus familiaris" = "org.Cf.eg.db",
    "Macaca mulatta"         = "org.Mmu.eg.db",
    "Pan troglodytes"        = "org.Pt.eg.db"
  )
}

#' @keywords internal
#' @noRd
.org_db_pkg_for_species <- function(species) {
  m <- .godb_species_map()
  if (!species %in% names(m)) {
    stop("source = 'godb' does not support species '", species, "'. ",
         "Supported: ", paste(names(m), collapse = ", "), ". Use ",
         "source = 'msigdbr' for other species.")
  }
  unname(m[species])
}

#' @keywords internal
#' @noRd
.fetch_go_via_godb <- function(ontology, species) {
  pkg <- .org_db_pkg_for_species(species)
  if (!requireNamespace(pkg, quietly = TRUE) ||
      !requireNamespace("GO.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("'", pkg, "', 'GO.db', and 'AnnotationDbi' are required for ",
         "source = 'godb'. Install with BiocManager::install(c('", pkg,
         "', 'GO.db', 'AnnotationDbi')).")
  }
  org_db <- getExportedValue(pkg, pkg)

  mapping <- AnnotationDbi::select(
    org_db, keys = AnnotationDbi::keys(org_db, keytype = "SYMBOL"),
    keytype = "SYMBOL", columns = c("GO", "ONTOLOGY")
  )
  mapping <- mapping[!is.na(mapping$GO) & mapping$ONTOLOGY == ontology, ]
  sets <- lapply(split(mapping$SYMBOL, mapping$GO), unique)

  terms <- AnnotationDbi::select(
    GO.db::GO.db, keys = names(sets), keytype = "GOID", columns = "TERM"
  )
  term_by_id <- setNames(terms$TERM, terms$GOID)
  names(sets) <- paste0(names(sets), "_", make.names(term_by_id[names(sets)]))
  sets
}
