#' @title HiCool utils
#' 
#' @name utils
#' @rdname HiCool
#' 
#' @description 
#' 
#' Utility functions provided by HiCool package 
#' 
#' @param hash Unique 6-letter ID used to identify files from a specific 
#'   HiCool processing run.
#' @param resolution Resolution used to import the mcool file
#' 
#' @importFrom HiCExperiment CoolFile
#' @importFrom stringr str_replace_all
#' @export
#' 
#' @examples 
#' importHiCoolFolder(hcf)

importHiCoolFolder <- function(output, hash, resolution = NULL) {
    files <- list.files(output, pattern = hash, full.names = TRUE, recursive = TRUE)
    HiCExperiment::CoolFile(
        path = grep('mcool', files, value = TRUE), 
        resolution = resolution,
        pairs = grep('pairs', files, value = TRUE), 
        metadata = list(
            hash = hash,
            log = grep('\\.log', files, value = TRUE)
        )
    )
}

.fixHOME <- function(x) {
    gsub('~', Sys.getenv('HOME'), x)
}

.getHicStats <- function(log, filtered, iterative) {
    stats <- list()
    lines <- readLines(log)
    
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
    if (grepl('.fa$|.fasta$', genome)) {
        if (!file.exists(genome)) {
            stop("Genome fasta file not found.")
        }
    }
    else {
        idx1 <- paste0(genome, '.1.bt2')
        idx2 <- paste0(genome, '.2.bt2')
        idx3 <- paste0(genome, '.3.bt2')
        idx4 <- paste0(genome, '.4.bt2')
        idx5 <- paste0(genome, '.rev.1.bt2')
        idx6 <- paste0(genome, '.rev.2.bt2')
        if (any(!file.exists(idx1, idx2, idx3, idx4, idx5, idx6))) {
            stop("Genome bowtie2 index detected, but some index files are missing.")
        }
    }
}