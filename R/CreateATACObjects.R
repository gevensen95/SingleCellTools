#' Create Seurat ATAC Objects
#'
#' This function creates multiple Seurat objects. It takes a list of directories
#' as input. In each directory, there should be at least the following files:
#' peaks.bed, singlecell.csv and fragments.tsv.gz. You can directly give the
#' output folder from cellranger into this function. While reading in the data,
#' the function creates a common peak list from the samples. Following this,
#' peaks on scaffolds and not the major chromosomes are removed.
#'
#'
#' @param data_dirs Path to the directories with ATAC data
#' @param add_treatment Whether to add a treatment column to metadata
#' @param treatment Treatment metadata value (vector)
#' @param genome One of \code{"mm10"} (default), \code{"hg38"}, or
#'   \code{"custom"}. \code{"mm10"}/\code{"hg38"} select the gene annotation
#'   (\code{EnsDb.Mmusculus.v79} / \code{EnsDb.Hsapiens.v86}) and genome
#'   sequence (\code{BSgenome.Mmusculus.UCSC.mm10} / \code{BSgenome.Hsapiens.UCSC.hg38})
#'   used for scaffold-chromosome filtering and the \code{ChromatinAssay} annotation.
#'   Both are optional (Suggests) packages -- only the one matching your chosen
#'   \code{genome} needs to be installed, and you'll get an actionable error naming
#'   it if it's missing. For any other species (e.g. rhesus macaque), pass
#'   \code{genome = "custom"} -- there's no way to auto-resolve an
#'   annotation/genome-sequence package for an arbitrary species the way
#'   there is for mouse/human, so this needs either \code{cellranger_ref}
#'   (recommended -- derives everything from your own cellranger reference
#'   package) or \code{annotation}/\code{main_chroms} supplied directly.
#'   See those parameters' docs below.
#' @param annotation Only used when \code{genome = "custom"} and
#'   \code{cellranger_ref} is NOT supplied. A \code{GRanges} of gene
#'   annotations, e.g. from \code{Signac::GetGRangesFromEnsDb()} on your own
#'   \code{EnsDb} (for a species without a dedicated \code{EnsDb.*} package,
#'   build one from an Ensembl GTF via \code{ensembldb::ensDbFromGtf()}).
#'   Must already use the same chromosome-naming style (e.g. \code{"chr1"}
#'   vs \code{"1"}) as your \code{peaks.bed}/\code{fragments.tsv.gz} files
#'   -- unlike the \code{"mm10"}/\code{"hg38"} paths, this is NOT
#'   auto-converted via \code{seqlevelsStyle()}, since there's no reliable
#'   way to know which style a custom genome's fragments are in.
#' @param main_chroms Only used when \code{genome = "custom"}. Character
#'   vector of the standard/main chromosome names to keep when filtering out
#'   scaffolds (e.g. \code{paste0("chr", c(1:20, "X", "Y"))} for rhesus
#'   rheMac10, matching whatever naming style your peaks use). For
#'   \code{"mm10"}/\code{"hg38"} this is derived automatically from the
#'   corresponding \code{BSgenome} package via
#'   \code{GenomeInfoDb::standardChromosomes()}. Without \code{cellranger_ref},
#'   this has to be passed explicitly -- there's no BSgenome dependency for a
#'   custom genome to derive it from. \strong{With \code{cellranger_ref}},
#'   this is normally derived automatically the same way (from the
#'   reference's own \code{fasta/*.fa} contig names, via
#'   \code{GenomeInfoDb::standardChromosomes()}) and does NOT need to be
#'   passed -- but MAY still be supplied explicitly as an override, and then
#'   takes precedence over that auto-detection. This override path exists
#'   because \code{standardChromosomes()}'s \code{chr1}/\code{1},
#'   \code{chrX}/\code{X}, \code{chrM}/\code{MT} heuristic doesn't recognize
#'   accession-style contig names (e.g. NCBI RefSeq \code{"NC_133406.1"}),
#'   which some genome assemblies -- including custom T2T builds -- use
#'   throughout, with no scaffolds to exclude at all; auto-detection then
#'   matches zero "standard" names and errors rather than guessing. Every
#'   name passed here is validated against the reference's own \code{.fai}
#'   index, so a typo or a name from the wrong assembly fails immediately
#'   with a clear error instead of silently keeping zero peaks on that
#'   "chromosome" later. (\code{annotation}, unlike \code{main_chroms}, stays
#'   mutually exclusive with \code{cellranger_ref} -- see that param's doc.)
#' @param genome_label Required when \code{genome = "custom"}. A short
#'   string identifying the genome build (e.g. \code{"rheMac10"}), used to
#'   tag the \code{ChromatinAssay}'s own genome slot the same way
#'   \code{"mm10"}/\code{"hg38"} do, and as the \code{genomeVersion} when
#'   \code{cellranger_ref} builds an \code{EnsDb} internally.
#' @param cellranger_ref Only used when \code{genome = "custom"}. Path to a
#'   reference package directory built with \code{cellranger-atac mkref} /
#'   \code{cellranger-arc mkref} (i.e. the same reference your
#'   \code{data_dirs} were actually aligned/called against). When supplied,
#'   \code{annotation} and \code{main_chroms} are derived automatically from
#'   files already inside it instead of being passed by hand: the gene
#'   annotation comes from \code{<cellranger_ref>/genes/genes.gtf.gz}, imported
#'   via \code{rtracklayer::import()} with seqlengths filled in from this
#'   genome's own \code{.fai} index (deliberately \emph{not}
#'   \code{ensembldb::ensDbFromGtf()} directly -- that leaves seqlengths NA
#'   and then tries to fetch them from Ensembl over FTP, which fails for a
#'   non-Ensembl custom assembly and/or without outbound internet), then
#'   \code{ensembldb::ensDbFromGRanges()} + \code{Signac::GetGRangesFromEnsDb()},
#'   and \code{main_chroms} comes from the contig names in
#'   \code{<cellranger_ref>/fasta/*.fa} (indexed with
#'   \code{Rsamtools::indexFa()} if a \code{.fai} isn't already there),
#'   filtered down via \code{GenomeInfoDb::standardChromosomes()}. This is
#'   the recommended way to use a non-mouse/human species like rhesus
#'   macaque: since the annotation and chromosome list both come from the
#'   exact same reference files cellranger used, there's no risk of a
#'   genome-build/assembly mismatch or a chromosome-naming-style mismatch
#'   against your actual \code{peaks.bed}/\code{fragments.tsv.gz} -- unlike
#'   downloading a GTF separately and hoping it matches. Requires
#'   \code{organism} (see below) and \code{genome_label} to also be
#'   supplied, since a cellranger reference's \code{genes.gtf.gz} isn't
#'   named the way \code{ensembldb} expects to auto-extract those from a
#'   filename. Mutually exclusive with passing \code{annotation} directly
#'   (ambiguous which gene annotation should win) -- but \code{main_chroms}
#'   MAY still be passed alongside \code{cellranger_ref}, as an override; see
#'   its own doc above.
#' @param organism Required when \code{genome = "custom"} and
#'   \code{cellranger_ref} is supplied. A short organism name/label (e.g.
#'   \code{"Macaca_mulatta"}) recorded in the \code{EnsDb} built from
#'   \code{cellranger_ref}'s GTF -- metadata only, not used for matching.
#' @param object_names Optional character vector of names for the returned list,
#'   the same length as \code{data_dirs}. \code{NULL} (default) uses
#'   \code{basename(data_dirs)}.
#' @param peakwidths_max Max peak width for finding combined peaks
#' @param peakwidths_min Min peak width for finding combined peaks
#' @param passed_filters_value Min value for filtering cells based on
#' passed_filters column
#' @param workers Number of parallel workers to use (via \code{future.apply})
#'   for building each sample's Seurat object -- reading singlecell.csv,
#'   the fragment file, computing the peak x cell FeatureMatrix, and the
#'   ChromatinAssay QC metrics, all fully independent across samples once
#'   the combined peak set is built. Defaults to \code{length(data_dirs)}
#'   (one worker per sample); errors up front if that (or an explicit
#'   value) exceeds \code{parallel::detectCores()}, naming the number of
#'   cores actually available. Pass \code{workers = 1} to run sequentially
#'   instead. \code{workers > 1} spins up that many parallel workers via
#'   \code{future::plan()} -- forked processes (\code{future::multicore})
#'   on Unix-likes outside RStudio, or background R sessions
#'   (\code{future::multisession}) on Windows / in RStudio, where forking
#'   isn't available -- restored on exit. Forked workers share memory with
#'   the main process via copy-on-write, but a \code{multisession} fallback
#'   holds its own copy of each sample's fragments/counts, so peak memory
#'   scales with \code{workers} in that case.
#' @return A list of Seurat objects
#' @export

CreateATACObjects <-
  function(data_dirs, add_treatment = FALSE, treatment = NULL,
           genome = c("mm10", "hg38", "custom"),
           annotation = NULL, main_chroms = NULL, genome_label = NULL,
           cellranger_ref = NULL, organism = NULL,
           object_names = NULL,
           peakwidths_max = 10000, peakwidths_min = 20,
           passed_filters_value = 500, workers = length(data_dirs)) {
    workers <- .resolve_workers(workers, n_samples = length(data_dirs),
                                was_default = missing(workers))
    genome <- match.arg(genome)

    if (workers > 1) {
      # Same BLAS/LAPACK thread-clamp as CreateRNAObjects() -- FeatureMatrix()
      # and the fragment-counting step underneath it can multithread
      # internally, so `workers` background sessions doing that at once would
      # otherwise oversubscribe the CPU the same way unclamped PCA/scaling
      # did there. unset = NA distinguishes "never set" from "set to empty
      # string" so restoring on exit uses Sys.unsetenv() rather than
      # Sys.setenv(VAR = ""), which would leave a literal empty-string value
      # behind (see CreateRNAObjects.R for the OMP_NUM_THREADS warning that
      # causes).
      old_blas_env <- Sys.getenv(c("VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS",
                                   "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"),
                                 unset = NA)
      Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1",
                OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
      on.exit({
        was_set <- !is.na(old_blas_env)
        if (any(was_set)) {
          do.call(Sys.setenv, as.list(old_blas_env[was_set]))
        }
        if (any(!was_set)) {
          Sys.unsetenv(names(old_blas_env)[!was_set])
        }
      }, add = TRUE)

      # See workers_utils.R -- shared by all six workers-taking loaders.
      cleanup <- .setup_future_plan(workers)
      on.exit(cleanup(), add = TRUE)
    }

    if(add_treatment == FALSE & is.null(treatment) == FALSE) {
      stop('\n\n  Error: Treatment vector was added, but add_treatment set to FALSE.\nSet add_treatment to TRUE before proceeding.')
    }
    if (isTRUE(add_treatment) && is.null(treatment)) {
      warning('add_treatment = TRUE but `treatment` is NULL -- no Treatment column will be added.')
    }
    if (!is.null(treatment) && length(treatment) != length(data_dirs)) {
      stop(sprintf(
        "`treatment` has length %d but there are %d `data_dirs` -- these must match one-to-one, or samples would silently get wrong/NA treatment labels via recycling.",
        length(treatment), length(data_dirs)))
    }
    if (!is.null(object_names) && length(object_names) != length(data_dirs)) {
      stop('`object_names` must be the same length as `data_dirs` (', length(data_dirs),
          '), got ', length(object_names), '.')
    }

    # ---- Resolve genome-specific annotation/sequence packages -----------------
    # Both mouse and human packages are optional (Suggests): only the one
    # matching the requested `genome` needs to actually be installed.
    #
    # For genome = "custom" there's no package to resolve at all -- an
    # arbitrary species (e.g. rhesus macaque) generally doesn't have a
    # dedicated EnsDb.*/BSgenome.* CRAN/Bioconductor package the way mouse
    # and human do (rhesus's own gene annotation, for instance, has to be
    # built from an Ensembl GTF via ensembldb::ensDbFromGtf(), then
    # Signac::GetGRangesFromEnsDb() on that). So this branch just validates
    # and uses whatever the caller supplied instead of trying to auto-resolve
    # anything -- see the annotation/main_chroms/genome_label docs above.
    if (identical(genome, "custom")) {
      if (is.null(genome_label) || !nzchar(genome_label)) {
        stop("genome = 'custom' requires `genome_label`, a short string ",
            "identifying the genome build (e.g. 'rheMac10') -- used to tag ",
            "the ChromatinAssay's own genome slot the same way 'mm10'/",
            "'hg38' do.")
      }

      if (!is.null(cellranger_ref)) {
        # ---- Derive annotation/main_chroms from a cellranger mkref reference ----
        # Recommended path: since both come from the exact reference files
        # cellranger used to align/call peaks, there's no genome-build or
        # chromosome-naming-style mismatch risk the way there is when
        # supplying `annotation`/`main_chroms` from a separately-downloaded
        # GTF -- see the cellranger_ref doc above.
        if (!is.null(annotation)) {
          stop("Pass either `cellranger_ref` OR `annotation` directly, not ",
              "both -- it's ambiguous which gene annotation should win. ",
              "`main_chroms` MAY still be supplied alongside `cellranger_ref` ",
              "-- see its doc -- as an override for when this genome's ",
              "contig names don't follow the chr1/1, chrX/X, chrM/MT ",
              "convention GenomeInfoDb::standardChromosomes() expects (e.g. ",
              "NCBI RefSeq accession-style names like 'NC_133406.1').")
        }
        if (is.null(organism) || !nzchar(organism)) {
          stop("genome = 'custom' with `cellranger_ref` requires `organism` ",
              "(e.g. 'Macaca_mulatta') -- a cellranger reference's ",
              "genes.gtf.gz isn't named the way ensembldb::ensDbFromGtf() ",
              "expects to auto-extract organism/genome/version from a ",
              "filename, so these have to be supplied explicitly.")
        }
        if (!requireNamespace("ensembldb", quietly = TRUE)) {
          stop("'ensembldb' is required for `cellranger_ref`. Install with: ",
              "BiocManager::install('ensembldb')")
        }
        if (!requireNamespace("Rsamtools", quietly = TRUE)) {
          stop("'Rsamtools' is required for `cellranger_ref`. Install with: ",
              "BiocManager::install('Rsamtools')")
        }
        if (!dir.exists(cellranger_ref)) {
          stop("`cellranger_ref` directory not found: '", cellranger_ref, "'")
        }

        genes_gtf <- file.path(cellranger_ref, "genes", "genes.gtf.gz")
        if (!file.exists(genes_gtf)) {
          stop("Expected a cellranger mkref reference at '", cellranger_ref,
              "' but didn't find '", genes_gtf, "'. Is this the top-level ",
              "reference package directory (the one containing 'fasta/' ",
              "and 'genes/' subdirectories)?")
        }
        fasta_candidates <- Sys.glob(file.path(cellranger_ref, "fasta", "*.fa"))
        if (length(fasta_candidates) != 1) {
          stop("Expected exactly one *.fa file under '",
              file.path(cellranger_ref, "fasta"), "', found ",
              length(fasta_candidates), ".")
        }
        fasta_path <- fasta_candidates[[1]]
        fai_path   <- paste0(fasta_path, ".fai")
        if (!file.exists(fai_path)) {
          message('--- Indexing reference FASTA (', basename(fasta_path), ') ---')
          Rsamtools::indexFa(fasta_path)
        }
        fai <- data.table::fread(fai_path, header = FALSE, data.table = FALSE,
                                 col.names = c("name", "length", "offset",
                                              "linebases", "linewidth"))

        if (!is.null(main_chroms)) {
          # Explicit override -- skip the naming-convention heuristic below
          # entirely. Needed for assemblies whose contigs use accession-style
          # names (e.g. NCBI RefSeq "NC_133406.1") that
          # GenomeInfoDb::standardChromosomes() doesn't recognize as
          # "standard" even when, as in that case, every contig in the file
          # genuinely is a real chromosome (or the mitochondrial genome) with
          # no scaffolds to exclude -- confirmed via the fasta index itself
          # (`cut -f1 <fai> | ...`), not assumed from the accession format
          # alone. Validated against the fasta index so a typo'd or
          # copy-pasted-from-elsewhere name fails loudly here rather than
          # silently keeping zero peaks on that "chromosome" downstream.
          unknown_chroms <- setdiff(main_chroms, fai$name)
          if (length(unknown_chroms) > 0) {
            stop("`main_chroms` includes name(s) not found in '", fasta_path,
                "': ", paste(unknown_chroms, collapse = ", "))
          }
          main.chroms <- main_chroms
          message(sprintf('--- Using user-supplied main_chroms (%d contigs) ---',
                          length(main.chroms)))
        } else {
          ref_seqinfo <- GenomeInfoDb::Seqinfo(seqnames   = fai$name,
                                               seqlengths = fai$length,
                                               genome     = genome_label)
          main.chroms <- GenomeInfoDb::standardChromosomes(ref_seqinfo)
          if (length(main.chroms) == 0) {
            stop(
              "GenomeInfoDb::standardChromosomes() matched none of the ",
              "contig names in '", fasta_path, "' -- e.g. accession-style ",
              "names like 'NC_133406.1' aren't recognized by its chr1/1, ",
              "chrX/X, chrM/MT heuristic. Pass `main_chroms` explicitly ",
              "instead (this IS supported alongside `cellranger_ref` -- only ",
              "`annotation` is mutually exclusive with it). If this assembly ",
              "has no scaffolds to exclude, use every contig in the index: ",
              "main_chroms = data.table::fread('", fai_path,
              "', header = FALSE)$V1"
            )
          }
        }

        if (!requireNamespace("rtracklayer", quietly = TRUE)) {
          stop("'rtracklayer' is required for `cellranger_ref`. Install with: ",
              "BiocManager::install('rtracklayer')")
        }

        message(sprintf(
          '--- Building EnsDb from %s (organism = %s, genomeVersion = %s) ---',
          genes_gtf, organism, genome_label))

        # ensembldb::ensDbFromGtf() imports the GTF with rtracklayer::import(),
        # which leaves every seqlength NA (a GTF has no sequence-length
        # header) -- ensDbFromGRanges() then notices the NAs and tries to
        # fill them in by fetching seqlengths from Ensembl over FTP. That
        # makes no sense for a non-Ensembl custom assembly like this one
        # (there's no real "release" to fetch -- the FTP paths it tries are
        # built from the placeholder `version = 1L` below) and fails outright
        # on a compute node without outbound internet, which then cascades
        # into a fatal, confusingly-unrelated "missing value where TRUE/FALSE
        # needed" error from ensembldb's own internal exon-ordering check
        # once it's left holding NA seqlengths. Importing the GTF ourselves
        # and setting real seqlengths -- from this genome's own .fai index,
        # the same one main.chroms was already validated against above --
        # before handing it to ensDbFromGRanges() avoids the fetch entirely:
        # it only fires when a seqlength is still NA.
        gtf_gr <- rtracklayer::import(genes_gtf, format = "gtf")
        fai_lengths <- setNames(fai$length, fai$name)
        missing_lengths <- setdiff(GenomeInfoDb::seqlevels(gtf_gr), names(fai_lengths))
        if (length(missing_lengths) > 0) {
          stop("'", genes_gtf, "' references contig(s) not found in this ",
              "genome's .fai index (", fai_path, "): ",
              paste(missing_lengths, collapse = ", "))
        }
        GenomeInfoDb::seqlengths(gtf_gr) <- fai_lengths[GenomeInfoDb::seqlevels(gtf_gr)]

        # ensembldb::ensDbFromGRanges()'s own internal validity check compares
        # each exon's GTF-stated `exon_number` against the order it computes
        # itself from genomic position; when `exon_number` is missing for any
        # exon, that comparison produces NA instead of TRUE/FALSE, and its own
        # `if (any(Different))` check crashes with a confusing "missing value
        # where TRUE/FALSE needed" -- e.g. every entry on a manually
        # converted/appended contig (mtDNA annotations merged in from a
        # different source than the rest of the GTF are a common way this
        # happens, since mitochondrial gene models are usually single-exon
        # and easy to lose an `exon_number` attribute for during conversion).
        # Rather than require every custom genome's GTF to already carry a
        # complete `exon_number` for every exon, compute it ourselves for
        # whatever's missing -- 1-based rank by genomic position within each
        # transcript, ascending for +/* strand and descending for - strand
        # (the standard 5'->3' exon-numbering convention) -- before handing
        # the GRanges to ensDbFromGRanges().
        is_exon <- gtf_gr$type == "exon"
        if (is.null(gtf_gr$exon_number)) {
          gtf_gr$exon_number <- rep(NA_character_, length(gtf_gr))
        }
        needs_exon_number <- is_exon & (is.na(gtf_gr$exon_number) |
                                        !nzchar(gtf_gr$exon_number))
        if (any(needs_exon_number)) {
          affected_tx <- unique(gtf_gr$transcript_id[needs_exon_number])
          message(sprintf(
            paste('--- %d exon(s) across %d transcript(s) had no exon_number',
                 'in the GTF -- computing one from genomic position ---'),
            sum(needs_exon_number), length(affected_tx)))
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

        ensdb_path <- ensembldb::ensDbFromGRanges(
          x             = gtf_gr,
          outfile       = tempfile(fileext = ".sqlite"),
          organism      = organism,
          genomeVersion = genome_label,
          version       = 1L
        )
        ensdb_obj   <- ensembldb::EnsDb(ensdb_path)
        # standard.chromosomes = FALSE -- GetGRangesFromEnsDb()'s own default
        # (TRUE) applies GenomeInfoDb::keepStandardChromosomes() internally,
        # which uses the same chr1/1, chrX/X, chrM/MT heuristic that doesn't
        # recognize this genome's NCBI RefSeq accession-style contig names
        # (e.g. "NC_133406.1") -- with every contig unrecognized, it filters
        # out all of them, and GetGRangesFromEnsDb() returns NULL instead of
        # a GRanges. main.chroms was already derived/validated from this same
        # reference's own .fai index above, so this redundant filtering is
        # both unnecessary and actively destructive here -- skip it.
        annotations <- Signac::GetGRangesFromEnsDb(ensdb = ensdb_obj,
                                                    standard.chromosomes = FALSE)
        # No seqlevelsStyle conversion here, deliberately -- annotations and
        # main.chroms both trace back to this same reference's own files, so
        # whatever naming convention it used, they already agree with each
        # other and with the peaks/fragments cellranger produced from it.
        genome_tag <- genome_label
      } else {
        # ---- Manual custom path: annotation/main_chroms supplied directly ----
        if (is.null(annotation) || !methods::is(annotation, "GRanges")) {
          stop("genome = 'custom' requires either `cellranger_ref`, or ",
              "`annotation` (a GRanges of gene annotations, e.g. from ",
              "Signac::GetGRangesFromEnsDb() on your own EnsDb) already in ",
              "the same chromosome-naming style as your peaks.bed/",
              "fragments.tsv.gz files.")
        }
        if (is.null(main_chroms) || length(main_chroms) == 0) {
          stop("genome = 'custom' requires either `cellranger_ref`, or ",
              "`main_chroms` (a character vector of the standard ",
              "chromosome names to keep, e.g. paste0('chr', c(1:20, 'X', ",
              "'Y')) for rhesus rheMac10) in the same naming style as your ",
              "peaks/fragments files -- there's no BSgenome package to ",
              "derive this from automatically for a custom genome.")
        }
        annotations <- annotation
        main.chroms <- main_chroms
        genome_tag  <- genome_label
      }
      message(sprintf('--- Using custom annotation/genome (%s, %d chromosomes) ---',
                      genome_tag, length(main.chroms)))
    } else {
      genome_pkgs <- switch(
        genome,
        mm10 = list(ensdb = 'EnsDb.Mmusculus.v79', bsgenome = 'BSgenome.Mmusculus.UCSC.mm10'),
        hg38 = list(ensdb = 'EnsDb.Hsapiens.v86',   bsgenome = 'BSgenome.Hsapiens.UCSC.hg38')
      )
      if (!requireNamespace(genome_pkgs$ensdb, quietly = TRUE)) {
        stop("'", genome_pkgs$ensdb, "' is required for genome = '", genome, "'. ",
            "Install with: BiocManager::install('", genome_pkgs$ensdb, "')")
      }
      if (!requireNamespace(genome_pkgs$bsgenome, quietly = TRUE)) {
        stop("'", genome_pkgs$bsgenome, "' is required for genome = '", genome, "'. ",
            "Install with: BiocManager::install('", genome_pkgs$bsgenome, "')")
      }
      # Both packages export an object of the same name as the package itself
      # (the standard Bioconductor EnsDb/BSgenome annotation-package convention).
      ensdb_obj    <- getExportedValue(genome_pkgs$ensdb, genome_pkgs$ensdb)
      bsgenome_obj <- getExportedValue(genome_pkgs$bsgenome, genome_pkgs$bsgenome)

      message(sprintf('--- Loading gene annotations from %s ---', genome_pkgs$ensdb))
      annotations <- Signac::GetGRangesFromEnsDb(ensdb = ensdb_obj)

      # EnsDb annotations come out in Ensembl seqname style ("1", "2", ...,
      # "MT"); ChromatinAssay/fragment data here is UCSC style ("chr1",
      # "chr2", ..., "chrM"). seqlevelsStyle<- handles the full mapping
      # (including the MT -> chrM special case), unlike a manual
      # paste0("chr", ...), which would produce the wrong "chrMT" and leave
      # every other non-numbered contig unmapped. Only done for the known
      # mm10/hg38 packages -- see the `annotation` doc above for why this
      # is NOT applied to a custom annotation.
      GenomeInfoDb::seqlevelsStyle(annotations) <- "UCSC"

      main.chroms <- GenomeInfoDb::standardChromosomes(bsgenome_obj)
      genome_tag  <- genome
    }
    GenomeInfoDb::genome(annotations) <- genome_tag

    message(sprintf('--- Reading peak sets (%d directories) ---', length(data_dirs)))
    # Use lapply to read the peak sets for each sample
    peak_data_list <- lapply(data_dirs, function(dir) {
      # Read peaks
      # data.table::fread() instead of read.table() -- peaks.bed can be a
      # large genome-wide peak set; fread's parser is substantially faster
      # for this than read.table()'s.
      peak_data <- data.table::fread(file = paste(dir, '/outs/peaks.bed', sep = ''),
                                     header = FALSE, col.names = c("chr", "start", "end"),
                                     data.table = FALSE)
      # Make GRanges objects
      gr <- GenomicRanges::makeGRangesFromDataFrame(peak_data)
    })

    message('--- Building combined peak set ---')
    # Create combined peak set. reduce() over the concatenation of every
    # sample's peaks at once -- NOT a for loop reassigning combined.peaks
    # from just peak_data_list[[1]] and peak_data_list[[i]] each iteration,
    # which discarded every previous iteration's result and silently kept
    # only sample 1 + the last sample's peaks for any run with >2 samples
    # (and errored outright with exactly 1 sample, since 2:length(x) counts
    # backward to 2:1 when length(x) == 1).
    combined.peaks <- GenomicRanges::reduce(do.call(c, peak_data_list))
    peakwidths <- GenomicRanges::width(combined.peaks)
    combined.peaks <- combined.peaks[peakwidths < peakwidths_max &
                                       peakwidths > peakwidths_min]
    message(sprintf('  Peaks within width range [%d, %d]: %d',
                    peakwidths_min, peakwidths_max, length(combined.peaks)))

    message('--- Removing peaks on scaffolds (keeping main chromosomes) ---')
    # remove scaffolds not in genome -- main.chroms/annotations/genome_tag
    # were already resolved above (mm10/hg38 auto-detected, or supplied
    # directly for genome = "custom").
    keep.peaks <- as.logical(GenomeInfoDb::seqnames(GenomicRanges::granges(combined.peaks)) %in% main.chroms)
    combined.peaks <- combined.peaks[keep.peaks, ]
    message(sprintf('  Peaks after scaffold removal: %d', length(combined.peaks)))

    message(sprintf('--- Building Seurat ATAC objects per sample%s ---',
                    if (workers > 1) sprintf(' (%d parallel workers)', workers) else ''))
    # Create Seurat objects. Fully independent per sample given the shared
    # combined.peaks/annotations built above, so this parallelizes cleanly
    # when workers > 1.
    .build_one <- function(idx) {
      dir <- data_dirs[[idx]]
      if (workers == 1) {
        message(sprintf('  Building object %d of %d: %s',
                        idx, length(data_dirs), basename(dir)))
      }

      # Load metadata for each sample. data.table::fread() instead of
      # read.table() -- singlecell.csv is one row per barcode (often
      # hundreds of thousands for ATAC) and read.table() is slow at that
      # size. fread() has no row.names arg (data.table has no rownames), so
      # column 1 is moved to rownames manually to match read.table(row.names
      # = 1)'s behavior.
      md <- data.table::fread(file = paste(dir, "/outs/singlecell.csv", sep = ''),
                              sep = ",", header = TRUE, data.table = FALSE)
      rownames(md) <- md[[1]]
      md[[1]] <- NULL
      md <- md[-1, ] # remove the first row
      md <- md[md$passed_filters > passed_filters_value, ]

      # Create fragment objects
      frag.obj <- Signac::CreateFragmentObject(path = paste(dir, '/outs/fragments.tsv.gz', sep = ''),
                                               cells = rownames(md))

      # Create Feature matrix objects
      counts <- Signac::FeatureMatrix(
        fragments = frag.obj,
        features = combined.peaks,
        cells = rownames(md)
      )

      # Create chromatin assay and final object with QC metrics. `genome_tag`
      # tags the assay's own seqinfo (previously only `annotations` got a
      # genome tag, leaving the assay itself untagged) -- "mm10"/"hg38" for
      # the built-in paths, or the caller's own `genome_label` for custom.
      assay <- Signac::CreateChromatinAssay(counts, fragments = frag.obj, genome = genome_tag)
      seurat.obj <- Seurat::CreateSeuratObject(assay, assay = "ATAC", meta.data = md,
                                               project = basename(dir))

      # add the gene information to the object
      Signac::Annotation(seurat.obj) <- annotations

      seurat.obj <- Signac::NucleosomeSignal(seurat.obj)
      seurat.obj$nucleosome_group <- ifelse(seurat.obj$nucleosome_signal > 4,
                                            'NS > 4', 'NS < 4')
      seurat.obj <- Signac::TSSEnrichment(seurat.obj)
      seurat.obj$pct_reads_in_peaks <- seurat.obj$peak_region_fragments /
        seurat.obj$passed_filters * 100
      seurat.obj$blacklist_ratio <- seurat.obj$blacklist_region_fragments /
        seurat.obj$peak_region_fragments

      return(seurat.obj)
    }

    seurat_objects <- if (workers > 1) {
      future.apply::future_lapply(seq_along(data_dirs), .build_one, future.seed = TRUE)
    } else {
      lapply(seq_along(data_dirs), .build_one)
    }

    names(seurat_objects) <- if (!is.null(object_names)) object_names else basename(data_dirs)

    message('--- Generating ATAC QC plots ---')
    # Row-bind just the metadata rather than merge()-ing the objects --
    # these QC plots never touch the fragment/counts data merge() would also
    # combine, and merging ChromatinAssay objects has its own fragment-path
    # gotchas that this sidesteps entirely (matching the fix already applied
    # to CreateRNAObjects.R).
    meta <- dplyr::bind_rows(lapply(seurat_objects, function(x) x@meta.data))
    orig.ident <- pct_reads_in_peaks <- peak_region_fragments <- NULL
    TSS.enrichment <- blacklist_ratio <- nucleosome_signal <- NULL  # silence R CMD check NSE notes

    pct_reads_in_peaks.plot <- ggplot2::ggplot(meta,
                                      ggplot2::aes(orig.ident, pct_reads_in_peaks)) + ggplot2::geom_boxplot() + Ol_Reliable()
    peak_region_fragments.plot <- ggplot2::ggplot(meta,
                                         ggplot2::aes(orig.ident, peak_region_fragments)) + ggplot2::geom_boxplot() + Ol_Reliable()
    TSS.enrichment.plot <- ggplot2::ggplot(meta,
                                  ggplot2::aes(orig.ident, TSS.enrichment)) + ggplot2::geom_boxplot() + Ol_Reliable()
    blacklist_ratio.plot <- ggplot2::ggplot(meta,
                                   ggplot2::aes(orig.ident, blacklist_ratio)) + ggplot2::geom_boxplot() + Ol_Reliable()
    nucleosome_signal.plot <- ggplot2::ggplot(meta,
                                     ggplot2::aes(orig.ident, nucleosome_signal)) + ggplot2::geom_boxplot() + Ol_Reliable()

    print(pct_reads_in_peaks.plot + peak_region_fragments.plot +
            TSS.enrichment.plot +
            blacklist_ratio.plot + nucleosome_signal.plot +
            patchwork::plot_layout(ncol = 3))

    # Add a column to metadata to specify treatment
    if (is.null(treatment) == FALSE ){
      message('--- Adding Treatment metadata column ---')
      seurat_objects <- setNames(lapply(seq_along(seurat_objects), function(i) {
        seurat_obj <- seurat_objects[[i]]
        seurat_obj[["Treatment"]] <- treatment[i]
        return(seurat_obj)
      }), names(seurat_objects))
    }

    return(seurat_objects)

  }
