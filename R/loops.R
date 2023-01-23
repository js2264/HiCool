#' @title Finding loops in contact map
#' @name getLoops
#' @rdname getLoops
#' @description
#' 
#' Find loops using chromosight
#'
#' @param x A `HiCExperiment` object
#' @param resolution Which resolution to use to search loops
#' @param output_prefix Prefix to chromosight output (default: "chromosight/chromo")
#' @param pearson Minimum Pearson correlation score to use to filter for 
#' significant loops 
#' @param min.dist,max.dist Min and max distance to use to filter for significant 
#' loops 
#' @param min.separation Minimum separation between anchors of potential loops
#' @param n.mads Number of MADs to use to filter relevant bins to search for 
#' loops 
#' @param nreads Number of reads to subsample to before searching for loops
#' @param ncores Number of cores for chromosight
#' @param norm Normalization parameter for chromosight
#' @return A `HiCExperiment` object with a new "loops" topologicalFeatures 
#' storing significant interactions identified by chromosight, and an additional 
#' `chromosight_args` metadata entry.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom vroom vroom
#' @importFrom InteractionSet GInteractions
#' @importFrom IRanges IRanges
#' @export
#' @examples 
#' contacts_yeast <- contacts_yeast()
#' contacts_yeast <- getLoops(contacts_yeast)
#' S4Vectors::metadata(contacts_yeast)$chromosight_args
#' topologicalFeatures(contacts_yeast, 'loops')

getLoops <- function(
    x, 
    resolution = NULL, 
    output_prefix = file.path('chromosight', 'chromo'), 
    norm = 'auto', 
    max.dist = 'auto', 
    min.dist = 'auto', 
    min.separation = 'auto', 
    n.mads = 5L, 
    pearson = 'auto', 
    nreads = 'no', 
    ncores = 1L  
) {
    proc <- basilisk::basiliskStart(env_HiCool)
    on.exit(basilisk::basiliskStop(proc))
    if (is.null(resolution)) resolution <- resolution(x)
    path <- paste0(fileName(x), '::/resolutions/', as.integer(resolution))
    dir.create(dirname(output_prefix), showWarnings = FALSE)
    tmpdir <- tempdir()
    args <- list(
        # non-settable parameters
        "--pattern" = "loops",
        "--dump" = tmpdir, 
        "--inter" = FALSE, 
        "--iterations" = 'auto', 
        "--kernel-config" = NULL, 
        "--perc-zero" = 'auto', 
        "--perc-undetected" = 'auto', 
        "--tsvd" = FALSE, 
        "--win-fmt" = 'json', 
        "--win-size" = "auto", 
        "--no-plotting" = TRUE, 
        "--smooth-trend" = FALSE,
        # settable parameters
        "--norm" = 'auto', 
        "<contact_map>" = path, 
        "--max-dist" = max.dist, 
        "--min-dist" = min.dist, 
        "--min-separation" = min.separation, 
        "--n-mads" = n.mads, 
        "<prefix>" = output_prefix,
        "--pearson" = pearson, 
        "--subsample" = nreads, 
        "--threads" = ncores 
    )
    loops <- basilisk::basiliskRun(
        env = env_HiCool, 
        fun = .getLoops,
        chromosight_args = args
    )
    topologicalFeatures(x, 'loops') <- loops
    metadata(x)[['chromosight_args']] <- args
    return(x)
}

.getLoops <- function(chromosight_args) {
    cs <- reticulate::import("chromosight")
    cs$cli$chromosight$cmd_detect(chromosight_args)
    df <- vroom::vroom(paste0(chromosight_args[['<prefix>']], '.tsv'), show_col_types = FALSE)
    loops <- InteractionSet::GInteractions(
        anchor1 = GenomicRanges::GRanges(
            df$chrom1, IRanges::IRanges(df$start1+1, df$end1)
        ),
        anchor2 = GenomicRanges::GRanges(
            df$chrom2, IRanges::IRanges(df$start2+1, df$end2)
        ),
        bin_id1 = df$bin1, 
        bin_id2 = df$bin2, 
        score = df$score, 
        pvalue = df$pvalue, 
        qvalue = df$qvalue
    )
    return(loops)
}