#' Build a Signac-compatible gene annotation from a raw GTF, and save it
#'
#' Runs the GTF -> \code{ensembldb::ensDbFromGRanges()} -> \code{EnsDb} ->
#' \code{Signac::GetGRangesFromEnsDb()} pipeline \code{\link{CreateATACObjects}}'s
#' \code{cellranger_ref} path uses internally, but standalone -- so you can
#' build the annotation once from your own genome-wide GTF, save it, and
#' reuse the saved file across every \code{CreateATACObjects()} call for
#' that genome (passing either the saved path or the returned \code{GRanges}
#' directly to its \code{annotation} argument) instead of relying on a
#' cellranger reference bundle's own \code{genes/genes.gtf.gz}.
#'
#' \strong{Why this exists, not just \code{cellranger_ref}.} A
#' cellranger-atac/cellranger-arc \code{mkref} reference bundle's
#' \code{genes/genes.gtf.gz} is not guaranteed to carry a \code{gene}-type
#' feature line for every chromosome that has \code{transcript}/\code{exon}
#' lines. One specific custom T2T macaque reference was found to have
#' \code{gene} rows on only 2 of 23 chromosomes despite having full
#' transcript/exon data on all of them -- \code{ensembldb::ensDbFromGRanges()}
#' reflected that faithfully, silently building a database with genes on
#' only those 2 chromosomes, no error and no obvious warning (the missing
#' genes simply weren't there to build from). That surfaced three steps
#' downstream as an apparently-unrelated dropped-seqlevel mismatch on the
#' final Seurat object, and took a full bisection (checking raw GTF import,
#' the EnsDb's own gene table, and Signac's wrapper on top of it
#' separately) to track back to the actual cause. This function runs the
#' same per-chromosome gene-count check up front -- see the
#' \code{transcript}/\code{exon}-vs-\code{gene} comparison in the source --
#' so that failure mode is caught immediately, with a warning naming the
#' affected chromosomes, instead of being discovered downstream.
#'
#' @param gtf_path Path to a genome-wide GTF file (e.g. from Ensembl, RefSeq,
#'   or your own annotation/liftover pipeline) -- NOT necessarily a
#'   cellranger reference's bundled \code{genes/genes.gtf.gz}, which may not
#'   carry a \code{gene}-type line for every chromosome (see Details).
#' @param fai_path Path to the \code{.fai} index of the genome FASTA this
#'   GTF's coordinates are relative to (e.g. \code{samtools faidx genome.fa}
#'   output, or the one bundled in a cellranger reference's \code{fasta/}
#'   directory). Used to fill in real seqlengths -- without them,
#'   \code{ensembldb::ensDbFromGRanges()} tries to fetch seqlengths from
#'   Ensembl over FTP, which fails on a custom assembly (there's no real
#'   Ensembl "release" for it) and/or on a compute node without outbound
#'   internet.
#' @param organism Organism name for the \code{EnsDb} metadata, e.g.
#'   \code{"Macaca_mulatta"}. Metadata only -- not used for matching.
#' @param genome_label A short genome-build label, e.g.
#'   \code{"T2T_MMU8v2_mtDNA"}. Stored as the resulting \code{GRanges}'
#'   \code{genome()} tag; must match whatever you pass
#'   \code{CreateATACObjects()}'s own \code{genome_label} later --
#'   \code{Annotation<-}'s underlying \code{SetAssayData()} errors on a
#'   genome-tag mismatch between an object and the annotation you attach to
#'   it.
#' @param save_path File path to save the resulting \code{GRanges} to via
#'   \code{saveRDS()} (e.g. \code{"macaque_t2t_annotation.rds"}) -- pick a
#'   name you'll recognize later. There's no default: this is meant to be a
#'   reusable, deliberately-named file you point \code{CreateATACObjects()}
#'   at going forward, not a throwaway temp file.
#' @param overwrite Overwrite \code{save_path} if it already exists?
#'   Default \code{FALSE} -- errors instead, so a typo'd path can't silently
#'   clobber a previously-built annotation you meant to keep.
#' @param verbose Message progress. Default \code{TRUE}.
#'
#' @return The annotation \code{GRanges}, invisibly (it's also written to
#'   \code{save_path}, which is the point -- most callers won't need the
#'   return value in the same session, but it's there so you can chain
#'   straight into
#'   \code{CreateATACObjects(..., annotation = BuildATACAnnotation(...))}
#'   without a round-trip through disk if you want to).
#'
#' @seealso \code{\link{CreateATACObjects}}, whose \code{annotation}
#'   argument accepts either the \code{GRanges} this returns directly, or
#'   the \code{save_path} you saved it to.
#' @export
BuildATACAnnotation <- function(gtf_path, fai_path, organism, genome_label,
                                save_path, overwrite = FALSE, verbose = TRUE) {
  if (!requireNamespace("ensembldb", quietly = TRUE)) {
    stop("'ensembldb' is required for BuildATACAnnotation(). Install with: ",
        "BiocManager::install('ensembldb')")
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("'rtracklayer' is required for BuildATACAnnotation(). Install with: ",
        "BiocManager::install('rtracklayer')")
  }
  if (!file.exists(gtf_path)) {
    stop("`gtf_path` not found: '", gtf_path, "'")
  }
  if (!file.exists(fai_path)) {
    stop("`fai_path` not found: '", fai_path, "'")
  }
  if (file.exists(save_path) && !isTRUE(overwrite)) {
    stop("'", save_path, "' already exists. Pass `overwrite = TRUE` to ",
        "replace it, or choose a different `save_path`.")
  }

  if (verbose) message(sprintf('--- Importing %s ---', gtf_path))
  gtf_gr <- rtracklayer::import(gtf_path, format = "gtf")

  fai <- data.table::fread(fai_path, header = FALSE, data.table = FALSE,
                           col.names = c("name", "length", "offset",
                                        "linebases", "linewidth"))
  fai_lengths <- setNames(fai$length, fai$name)
  missing_lengths <- setdiff(GenomeInfoDb::seqlevels(gtf_gr), names(fai_lengths))
  if (length(missing_lengths) > 0) {
    stop("'", gtf_path, "' references contig(s) not found in '", fai_path,
        "': ", paste(missing_lengths, collapse = ", "))
  }
  GenomeInfoDb::seqlengths(gtf_gr) <- fai_lengths[GenomeInfoDb::seqlevels(gtf_gr)]

  # --- Sanity check: every chromosome with transcript/exon records also has
  # at least one `gene`-type record. This is exactly the failure mode
  # described in Details above -- gene rows on only a couple of chromosomes
  # despite full transcript/exon data on all of them -- which
  # ensembldb::ensDbFromGRanges() reflects faithfully (silently building a
  # database with genes on only those chromosomes) rather than erroring, so
  # it would otherwise surface three steps later as an apparently-unrelated
  # dropped-seqlevel mismatch on the final Seurat object instead of here,
  # where the actual cause -- and which chromosomes are affected -- is
  # obvious.
  gene_rows <- gtf_gr[gtf_gr$type == "gene"]
  chroms_with_genes <- unique(as.character(GenomicRanges::seqnames(gene_rows)))
  chroms_with_transcripts <- unique(as.character(GenomicRanges::seqnames(
    gtf_gr[gtf_gr$type %in% c("transcript", "exon")])))
  missing_gene_chroms <- setdiff(chroms_with_transcripts, chroms_with_genes)
  if (length(missing_gene_chroms) > 0) {
    warning(sprintf(paste0(
      "'%s' has transcript/exon records on %d chromosome(s) with NO ",
      "`gene`-type record at all: %s. ensembldb::ensDbFromGRanges() will ",
      "silently build genes for every OTHER chromosome only -- if that's ",
      "not what you expect (e.g. if this is a cellranger reference's ",
      "bundled genes.gtf.gz), use a different/more complete GTF instead."),
      gtf_path, length(missing_gene_chroms),
      paste(missing_gene_chroms, collapse = ", ")))
  }

  # exon_number gapfill -- see CreateATACObjects.R's cellranger_ref branch
  # for the full explanation of why ensembldb::ensDbFromGRanges() needs
  # this: when exon_number is missing for any exon, its own internal
  # exon-ordering validity check produces NA instead of TRUE/FALSE and
  # crashes with a confusing "missing value where TRUE/FALSE needed".
  # 1-based rank by genomic position within each transcript, ascending for
  # +/* strand and descending for - strand (the standard 5'->3' convention).
  is_exon <- gtf_gr$type == "exon"
  if (is.null(gtf_gr$exon_number)) {
    gtf_gr$exon_number <- rep(NA_character_, length(gtf_gr))
  }
  needs_exon_number <- is_exon & (is.na(gtf_gr$exon_number) |
                                  !nzchar(gtf_gr$exon_number))
  if (any(needs_exon_number)) {
    affected_tx <- unique(gtf_gr$transcript_id[needs_exon_number])
    if (verbose) {
      message(sprintf(
        paste('--- %d exon(s) across %d transcript(s) had no exon_number',
             'in the GTF -- computing one from genomic position ---'),
        sum(needs_exon_number), length(affected_tx)))
    }
    fix_idx  <- which(is_exon & gtf_gr$transcript_id %in% affected_tx)
    tx_id    <- gtf_gr$transcript_id[fix_idx]
    is_minus <- as.character(GenomicRanges::strand(gtf_gr))[fix_idx] == "-"
    pos_key  <- ifelse(is_minus, -GenomicRanges::start(gtf_gr)[fix_idx],
                       GenomicRanges::start(gtf_gr)[fix_idx])
    ordered_idx <- fix_idx[order(tx_id, pos_key)]
    gtf_gr$exon_number[ordered_idx] <- as.character(ave(
      seq_along(ordered_idx), gtf_gr$transcript_id[ordered_idx],
      FUN = seq_along
    ))
  }

  if (verbose) {
    message(sprintf('--- Building EnsDb (organism = %s, genomeVersion = %s) ---',
                    organism, genome_label))
  }
  ensdb_path <- ensembldb::ensDbFromGRanges(
    x             = gtf_gr,
    outfile       = tempfile(fileext = ".sqlite"),
    organism      = organism,
    genomeVersion = genome_label,
    version       = 1L
  )
  ensdb_obj <- ensembldb::EnsDb(ensdb_path)

  # standard.chromosomes = FALSE -- the default (TRUE) applies
  # GenomeInfoDb::keepStandardChromosomes() internally, which doesn't
  # recognize accession-style contig names (e.g. NCBI RefSeq
  # "NC_133406.1") and would filter out every one of them, returning NULL
  # instead of a GRanges.
  annotations <- Signac::GetGRangesFromEnsDb(ensdb = ensdb_obj,
                                             standard.chromosomes = FALSE)
  GenomeInfoDb::genome(annotations) <- genome_label

  if (verbose) {
    message(sprintf('--- Built annotation: %d seqlevel(s), %d feature(s) ---',
                    length(GenomeInfoDb::seqlevels(annotations)),
                    length(annotations)))
  }

  saveRDS(annotations, save_path)
  if (verbose) message(sprintf('--- Saved to %s ---', save_path))

  invisible(annotations)
}
