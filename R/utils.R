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

#' @export

getHicStats <- function(log) {
    stats <- list()
    lines <- readLines(log)
    filtered <- any(grepl('INFO :: Filtering with thresholds', lines))

    ## -- Parse log file
    if (!filtered) {
        nFragments <- grep("INFO :: .* mapped with Q", lines, value = TRUE) |>
            stringr::str_replace_all(".*/|\\)", "") |>
            as.numeric() |>
            (`/`)(2)
        nDangling <- 0
        nSelf <- 0
        nDumped <- 0
        nPairs <- grep("INFO :: .* PCR duplicates have", lines, value = TRUE) |>
            stringr::str_replace_all(".* / | pairs\\) $", "") |>
            as.numeric()
        nFiltered <- nPairs
        nDups <- grep("INFO :: .* PCR duplicates have", lines, value = TRUE) |>
            stringr::str_replace_all(".*out \\(| \\/ .*", "") |>
            as.numeric()
        nUnique <- nFiltered - nDups
        threshold_uncut <- NA
        threshold_self <- NA
    }
    else {
        nFragments <- grep("INFO :: .* mapped with Q", lines, value = TRUE) |>
            stringr::str_replace_all(".*/|\\)", "") |>
            as.numeric() |>
            (`/`)(2)
        nDangling <- grep("INFO :: .* pairs discarded", lines, value = TRUE) |>
            stringr::str_replace_all(".* Uncuts: |, Weirds.*", "") |>
            as.numeric() # "Uncut"
        nSelf <- grep("INFO :: .* pairs discarded", lines, value = TRUE) |>
            stringr::str_replace_all(".* Loops: |, Uncuts.*", "") |>
            as.numeric() # "loop"
        nDumped <- grep("INFO :: .* pairs discarded", lines, value = TRUE) |>
            stringr::str_replace_all(".* Weirds:", "") |>
            as.numeric() # "Weird"
        nFiltered <- grep("INFO :: .* pairs kept", lines, value = TRUE) |>
            stringr::str_replace_all(".*INFO ::| pairs.*", "") |>
            as.numeric()
        nPairs <- nFiltered + nDangling + nSelf + nDumped
        nDups <- grep("INFO :: .* PCR duplicates have", lines, value = TRUE) |>
            stringr::str_replace_all(".*out \\(| \\/ .*", "") |>
            as.numeric()
        nUnique <- nFiltered - grep("INFO :: .* PCR duplicates have", lines, value = TRUE) |>
            stringr::str_replace_all(".*out \\(| \\/ .*", "") |>
            as.numeric()
        threshold_uncut <- grep("INFO :: Filtering with thresholds", lines, value = TRUE) |>
            stringr::str_replace_all(".*thresholds: | loops=.*", "") |>
            stringr::str_replace(".*=", "") |>
            as.numeric()
        threshold_self <- grep("INFO :: Filtering with thresholds", lines, value = TRUE) |>
            stringr::str_replace_all(".*=", "") |>
            as.numeric()
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
