#' @rdname HiCool
#' @export

importHiCoolFolder <- function(output, hash, resolution = NULL) {
    files <- list.files(output, pattern = hash, full.names = TRUE, recursive = TRUE)
    log_file <- grep('\\.log', files, value = TRUE)
    mcool_file <- grep('mcool', files, value = TRUE)
    pairs_file <- grep('pairs', files, value = TRUE)

    ## --- Check that all required fields exist
    if (!length(mcool_file))
        stop("No cool file was found.")
    if (!length(pairs_file))
        message("Warning: No pairs file was found.")
    if (!file.exists(mcool_file))
        stop("The associated cool file is missing.")
    if (!file.exists(pairs_file))
        message("Warning: The associated pairs file is missing.")
    if (is.null(log_file))
        stop("Missing log file.")
    if (!file.exists(log_file))
        stop("Log file does not exist.")

    x <- HiCExperiment::CoolFile(
        path = mcool_file, 
        resolution = resolution,
        pairs = pairs_file, 
        metadata = list(
            log = log_file, 
            args = getHiCoolArgs(log_file),
            stats = getHicStats(log_file)
        )
    )
    return(x)
}

#' @rdname HiCool
#' @export

getHiCoolArgs <- function(log) {
    args <- list()
    lines <- readLines(log)
    fromIdx <- which(grepl('HiCool working directory :::', lines)) + 1
    toIdx <- which(grepl('^----------------$', lines)) - 1
    if (length(fromIdx)) {
        wd <- gsub('.*::: ', '', lines[fromIdx-1])
        argsl <- strsplit(gsub('.*::: ', '', lines[fromIdx:toIdx]), ': ')
        args <- lapply(argsl, '[', 2)
        names(args) <- lapply(argsl, '[', 1)
        args <- lapply(args, function(x) {
            ifelse(x == "TRUE" | x == "FALSE", as.logical((x)), x)
        })
        args$threads <- as.numeric(args$threads)
        args$wd <- wd
    }
    else {
        message("Warning: HiCool arguments could not be retrieved from the log file.")
        args <- list()
    }
    return(args)
}

#' @rdname HiCool
#' @export

getHicStats <- function(log) {
    stats <- list()
    lines <- readLines(log)
    filtered <- any(grepl('INFO :: Filtering with thresholds', lines))

    ## -- Parse log file. Numbers are formatted differently in logs from
    ## -- hicstuff >= 3.2.5, e.g. `(613/53553 pairs)` instead of
    ## -- `(613 / 53553 pairs) `: both are supported.
    mapped <- .logNumbers(lines, "mapped with Q >= [0-9]+ \\(([0-9]+) ?/ ?([0-9]+)\\)")
    pcr <- .logNumbers(lines, "PCR duplicates have been filtered out \\(([0-9]+) ?/ ?([0-9]+) pairs\\)")
    nFragments <- mapped[2] / 2
    nDups <- pcr[1]
    if (!filtered) {
        nDangling <- 0
        nSelf <- 0
        nDumped <- 0
        nFiltered <- pcr[2]
        threshold_uncut <- NA
        threshold_self <- NA
    }
    else {
        discarded <- .logNumbers(lines, "pairs discarded: Loops: ([0-9]+), Uncuts: ([0-9]+), Weirds: ([0-9]+)")
        nSelf <- discarded[1] # "loop"
        nDangling <- discarded[2] # "Uncut"
        nDumped <- discarded[3] # "Weird"
        nFiltered <- .logNumbers(lines, "INFO :: ([0-9]+) pairs kept")[1]
        thresholds <- .logNumbers(lines, "Filtering with thresholds: uncuts=([0-9]+) loops=([0-9]+)")
        threshold_uncut <- thresholds[1]
        threshold_self <- thresholds[2]
    }

    stats[["nFragments"]] <- nFragments
    stats[["nPairs"]] <- nFiltered + nDangling + nSelf + nDumped
    stats[["nDangling"]] <- nDangling # "Uncut"
    stats[["nSelf"]] <- nSelf # "loop"
    stats[["nDumped"]] <- nDumped # "Weird"
    stats[["nFiltered"]] <- nFiltered
    stats[["nDups"]] <- nDups
    stats[["nUnique"]] <- nFiltered - nDups
    stats[["threshold_uncut"]] <- threshold_uncut
    stats[["threshold_self"]] <- threshold_self

    return(stats)
}

## Numbers captured by the groups of `pattern` in the first matching log line
.logNumbers <- function(lines, pattern) {
    hit <- grep(pattern, lines, value = TRUE)
    if (!length(hit)) return(NA_real_)
    as.numeric(regmatches(hit[1], regexec(pattern, hit[1]))[[1]][-1])
}

## Run a command-line tool of the HiCool conda environment in a separate
## process. Its python modules are thus never loaded in the R session, where
## they can clash with libraries that R has already loaded.
.runInEnv <- function(env_dir, command, args, log) {
    bin <- file.path(env_dir, 'bin')
    vars <- c(
        PATH = paste(bin, Sys.getenv('PATH'), sep = .Platform$path.sep),
        PYTHONNOUSERSITE = '1',
        PYTHONPATH = '',
        PYTHONHOME = '',
        MPLBACKEND = 'Agg'
    )
    status <- system2(
        file.path(bin, command), shQuote(as.character(args)),
        stdout = log, stderr = log,
        env = paste0(names(vars), '=', shQuote(vars))
    )
    if (status != 0) {
        out <- if (file.exists(log)) utils::tail(readLines(log), 20)
        stop(
            "HiCool :: `", command, "` failed (exit status ", status, "). ",
            "Last lines of its output (", log, "):\n",
            paste(out, collapse = '\n'),
            call. = FALSE
        )
    }
    invisible(log)
}

## Command-line arguments from a list of docopt arguments (as used by
## chromosight): options set to TRUE are flags, options set to FALSE or NULL
## are dropped, and positional arguments (`<...>`) come last, in order.
.docoptArgs <- function(args) {
    positional <- grepl('^<.*>$', names(args))
    opts <- args[!positional]
    opts <- opts[!vapply(opts, function(x) is.null(x) || isFALSE(x), logical(1))]
    flags <- vapply(opts, isTRUE, logical(1))
    values <- vapply(opts[!flags], function(x) {
        if (is.numeric(x)) format(x, scientific = FALSE, trim = TRUE)
        else as.character(x)
    }, character(1))
    c(
        names(opts)[flags],
        paste0(names(values), '=', values, recycle0 = TRUE),
        unlist(args[positional], use.names = FALSE)
    )
}

.dhms <- function(t) {
    paste(
        t %/% (60*60*24), ' days, ',
        paste(
            formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0"), 'h ',
            formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"), 'min ',
            formatC(t %% 60, width = 2, format = "d", flag = "0"), 's ',
            sep = ""
        )
    )
}

.checkGenome <- function(genome) {
    
    ## -- Fetch genome bowtie2 index from refgenie
    if (genome %in% c('mm10', 'hg38', 'dm6')) {
        bfc <- BiocFileCache::BiocFileCache()
        rid_index <- BiocFileCache::bfcquery(bfc, query = paste0('HiCool_', genome))$rid
        if (!length(rid_index)) {
            message( "HiCool :: Fetching bowtie genome index archive from regenie..." )
            archive <- dplyr::case_when(
                genome == 'hg38' ~ "http://refgenomes.databio.org/v3/assets/archive/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/bowtie2_index?tag=default", 
                genome == 'mm10' ~ "http://refgenomes.databio.org/v3/assets/archive/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1/bowtie2_index?tag=default",
                genome == 'dm6' ~ "http://refgenomes.databio.org/v3/assets/archive/8baf9d24ad8f5678f0fe1f5b21a812d410755d49e3123158/bowtie2_index?tag=default"
            )
            bfcentry <- BiocFileCache::bfcadd( 
                bfc, 
                rname = paste0('HiCool_', genome), 
                fpath = archive 
            )
            rid_index <- names(bfcentry)
        }
        tmp_dir <- tempdir()
        message( "HiCool :: Unzipping bowtie2 genome index from refgenie..." )
        utils::untar(BiocFileCache::bfcrpath(bfc, rids = rid_index), exdir = tmp_dir)
        idx.files <- list.files(file.path(tmp_dir, 'default'), full.names = TRUE)
        idx.base <- gsub('\\..*', '', basename(idx.files)[1])
        for (file in idx.files) {
            file.rename(file, gsub(idx.base, genome, file))
        }
        genome <- file.path(dirname(idx.files[1]), genome)
        .checkGenome(genome)
    }

    ## -- Fetch genome bowtie2 index from AWS S3 iGenomes
    if (genome %in% c('R64-1-1', 'WBcel235', 'GRCz10', 'Galgal4')) {
        bfc <- BiocFileCache::BiocFileCache()
        rid_indices <- BiocFileCache::bfcquery(bfc, query = paste0('HiCool_', genome))$rid
        if (!length(rid_indices)) {
            message( "HiCool :: Fetching bowtie genome index files from AWS iGenomes S3 bucket..." )
            s3url <- "https://ngi-igenomes.S3.amazonaws.com/"
            S3basepath <- dplyr::case_when(
                genome == 'R64-1-1' ~ "igenomes/Saccharomyces_cerevisiae/Ensembl/R64-1-1/Sequence/Bowtie2Index/", 
                genome == 'WBcel235' ~ "igenomes/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/Bowtie2Index/", 
                genome == 'GRCz10' ~ "igenomes/Danio_rerio/Ensembl/GRCz10/Sequence/Bowtie2Index/", 
                genome == 'Galgal4' ~ "igenomes/Gallus_gallus/Ensembl/Galgal4/Sequence/Bowtie2Index/"
            )
            for (idx in c('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2')) {
                path <- paste0(s3url, S3basepath, 'genome', idx)
                BiocFileCache::bfcadd( 
                    bfc, 
                    rname = paste0('HiCool_', genome, idx), 
                    fpath = path 
                )
            }
            rid_indices <- BiocFileCache::bfcquery(bfc, query = paste0('HiCool_', genome))$rid
        }
        tmp_dir <- tempdir()
        message( "HiCool :: Recovering bowtie2 genome index from AWS iGenomes..." )
        idx.files <- BiocFileCache::bfcpath(bfc, rid_indices)
        for (file in idx.files) {
            idx.base <- gsub('\\..*', '', basename(file))
            file.copy(file, file.path(tmp_dir, basename(gsub(idx.base, genome, file))))
        }
        genome <- file.path(tmp_dir, genome)
        genome <- .checkGenome(genome)
        genome
    }

    ## -- Fetch local fasta file
    if (grepl('.fa$|.fasta$', genome)) {
        if (!file.exists(genome)) {
            stop("Genome fasta file not found.")
        }
        else {
            return(normalizePath(genome))
        }
    }

    ## -- Fetch local bowtie2 index files
    idx1 <- paste0(genome, '.1.bt2')
    idx2 <- paste0(genome, '.2.bt2')
    idx3 <- paste0(genome, '.3.bt2')
    idx4 <- paste0(genome, '.4.bt2')
    idx5 <- paste0(genome, '.rev.1.bt2')
    idx6 <- paste0(genome, '.rev.2.bt2')
    if (all(file.exists(idx1, idx2, idx3, idx4, idx5, idx6))) {
        genome <- gsub('.1.bt2', '', (normalizePath(idx1)))
        return(genome)
    }
    stop("Genome bowtie2 index detected, but some index files are missing.")
}
