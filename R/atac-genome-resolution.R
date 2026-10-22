# ============================================================================
# .resolve_atac_genome() -- shared genome-resolution logic for every
# function that builds a ChromatinAssay from a `genome` argument of
# "mm10" / "hg38" / "custom" (currently CreateATACObjects() and
# CreateATACObjectsFilter()). Originally written once, inline, inside
# CreateATACObjects() -- extracted here verbatim (no logic changes) so
# CreateATACObjectsFilter() can share it instead of re-implementing (and
# inevitably drifting from) the custom-genome path, which is the
# considerably more involved of the two: deriving gene annotation and
# standard-chromosome names either from a cellranger mkref reference
# directory (preferred -- see `cellranger_ref`'s doc on the two exported
# functions) or from caller-supplied `annotation`/`main_chroms`, with no
# BSgenome/EnsDb package to fall back on for an arbitrary non-model
# species. See CreateATACObjects.R's own parameter docs (`genome`,
# `annotation`, `main_chroms`, `genome_label`, `cellranger_ref`,
# `organism`) for the full rationale behind each piece of this -- they
# apply identically here.
#
# Returns list(annotations = <GRanges>, main.chroms = <character>,
# assay_genome = <Seqinfo or character scalar>) -- `assay_genome` is a
# full Seqinfo object for genome = "custom" (registering a genome label
# that genuinely isn't in UCSC/NCBI's registry) and the plain "mm10"/
# "hg38" string otherwise (which Signac::CreateChromatinAssay()'s
# `genome` argument resolves via that registry directly) -- see the
# comment right before `assay_genome <- GenomeInfoDb::Seqinfo(...)` below
# for why those two cases can't be handled the same way.
# ============================================================================

#' @keywords internal
#' @noRd
.resolve_atac_genome <- function(genome, annotation, main_chroms,
                                  genome_label, cellranger_ref, organism) {
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
        if (is.null(annotation)) {
          stop("genome = 'custom' requires either `cellranger_ref`, or ",
              "`annotation` (a GRanges of gene annotations -- e.g. from ",
              "Signac::GetGRangesFromEnsDb() on your own EnsDb, or ",
              "BuildATACAnnotation()'s return value -- or a path to an ",
              ".rds file containing one, e.g. what BuildATACAnnotation() ",
              "saves via its `save_path` argument) already in the same ",
              "chromosome-naming style as your peaks.bed/fragments.tsv.gz ",
              "files.")
        }
        # `annotation` may be given as a path (character) instead of a
        # GRanges directly -- e.g. the save_path BuildATACAnnotation() wrote
        # to. Load it here so the rest of this branch can treat `annotation`
        # uniformly as a GRanges from this point on. Keep the original path
        # string around (annotation_path) purely for error messages, since
        # `annotation` itself gets overwritten by readRDS()'s result below.
        if (is.character(annotation)) {
          if (length(annotation) != 1) {
            stop("`annotation` was given as a character vector of length ",
                length(annotation), " -- it must be a single file path ",
                "(or a GRanges object directly), not a vector of paths.")
          }
          annotation_path <- annotation
          if (!file.exists(annotation_path)) {
            stop("`annotation` was given as a path, but no file exists at '",
                annotation_path, "'.")
          }
          annotation <- readRDS(annotation_path)
          if (!methods::is(annotation, "GRanges")) {
            stop("The file at '", annotation_path, "' didn't contain a ",
                "GRanges object (got a `", class(annotation)[1], "`) -- ",
                "`annotation` must be a GRanges, or a path to an .rds file ",
                "(e.g. saved by BuildATACAnnotation()) containing one.")
          }
        } else if (!methods::is(annotation, "GRanges")) {
          stop("`annotation` must be a GRanges object, or a path to an ",
              ".rds file containing one (e.g. saved by ",
              "BuildATACAnnotation()) -- got a `", class(annotation)[1],
              "`.")
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

      # Signac::CreateChromatinAssay()'s `genome` argument, when given a
      # plain string, calls GenomeInfoDb::Seqinfo(genome = <string>) to
      # build a Seqinfo automatically -- which only works for a genome ID
      # that's registered with UCSC/NCBI ("mm10", "hg38", ...). genome_tag
      # here is this build's own genome_label, which isn't registered
      # anywhere, so that call fails downstream (per-sample, inside the
      # parallel workers, where the error is much harder to trace back to
      # this line) with '"<genome_label>" is not a registered NCBI
      # assembly or UCSC genome'. GenomeInfoDb::Seqinfo() only makes that
      # registry lookup when `seqnames` is NULL -- passing `seqnames`
      # explicitly (main.chroms, already validated above) builds the
      # Seqinfo directly instead, tagged with genome_tag as a plain label
      # rather than something it tries to resolve. Real seqlengths are
      # included when available (from this genome's own .fai index, via
      # `cellranger_ref`); NA seqlengths are fine otherwise -- Seqinfo
      # supports them and CreateChromatinAssay doesn't require real ones.
      genome_lengths <- if (exists("fai_lengths", inherits = FALSE)) {
        unname(fai_lengths[main.chroms])
      } else {
        unname(GenomeInfoDb::seqlengths(annotations)[main.chroms])
      }
      assay_genome <- GenomeInfoDb::Seqinfo(seqnames   = main.chroms,
                                            seqlengths = genome_lengths,
                                            genome     = genome_tag)
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
      # "mm10"/"hg38" ARE registered UCSC genome IDs, so the plain string
      # works fine as CreateChromatinAssay()'s `genome` argument here --
      # see the custom-genome branch above for why that's not true there.
      assay_genome <- genome_tag
    }
    GenomeInfoDb::genome(annotations) <- genome_tag

  list(annotations = annotations, main.chroms = main.chroms,
       assay_genome = assay_genome)
}
