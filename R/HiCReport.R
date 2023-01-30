#' @title HiC processing report
#' 
#' @name HiCReport
#' @rdname HiCReport
#' 
#' @param x an CoolFile object, generated from `HiCool::HiCool()` or 
#'   `HiCool::importHiCoolFolder()`, or directly from calling 
#'   `HiCExperiment::CoolFile()`.
#' @param output Path to save output HTML file.
#'
#' @importFrom rmarkdown render
#' @importFrom rmdformats readthedown
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom sessioninfo os_name
#' @importFrom sessioninfo package_info
#' @importFrom sessioninfo python_info
#' @importFrom BiocIO resource
#' @importFrom S4Vectors metadata
#' @return String to the generated HTML report file
#' @export
#' @examples 
#' mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#' log_path <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'HiCool_log')
#' cf <- CoolFile(mcool_path, pairs = pairs_path, metadata = list(log = log_path))
#' HiCReport(cf)

HiCReport <- function(x, output = NULL) {

    ## --- Parse CoolFile for info
    mcool <- normalizePath(BiocIO::resource(x))
    pairs <- normalizePath(pairsFile(x))
    log_file <- normalizePath(S4Vectors::metadata(x)$log)

    ## --- Check that all required fields exist
    if (!length(mcool))
        stop("No cool file was found.")
    if (!length(pairs))
        message("Warning: No pairs file was found.")
    if (!file.exists(mcool))
        stop("The associated cool file is missing.")
    if (!file.exists(pairs))
        message("Warning: The associated pairs file is missing.")
    if (is.null(log_file))
        stop("Missing log file.")
    if (!file.exists(log_file))
        stop("Log file does not exist.")
    
    ## --- Prepare output files
    if (is.null(output)) {
        output <- file.path(
            dirname(dirname(mcool)), 
            gsub('\\..?cool', '.html', basename(mcool))
        )
        if (!grepl('.html', output))
            output <- paste0(output, '.html')
    }
    if (file.exists(output)) {
        if (interactive()) {
            qtxt <- sprintf(
                "Output report %s already exists. \nOverwrite it? [Y/n]: ",
                output
            )
            repeat {
                cat(qtxt)
                answer <- readLines(n = 1)
                if (answer %in% c("y", "Y", "n", "N")) break
            }
            tolower(answer)
            if ("n" == answer)
                stop("Aborting report computing now.")
        }
        file.remove(output)
    }

    ## --- Generate report
    tmpdir <- tempdir()
    tmpRmd <- file.path(tmpdir, basename(gsub('.html$', '.Rmd', output)))
    tmpHtml <- file.path(tmpdir, basename(output))
    if (file.exists(tmpRmd)) file.remove(tmpRmd)
    if (file.exists(tmpHtml)) file.remove(tmpHtml)
    file.copy(
        system.file('templates', 'HiCoolFile_report_template.Rmd', package = 'HiCool'), 
        tmpRmd
    )
    writeLines(gsub("%COOL%", mcool, readLines(tmpRmd)), tmpRmd)
    writeLines(gsub("%PAIRS%", pairs, readLines(tmpRmd)), tmpRmd)
    writeLines(gsub("%LOG%", log_file, readLines(tmpRmd)), tmpRmd)
    rmarkdown::render(tmpRmd, quiet = TRUE)
    writeLines(gsub("&amp;quot;", '"', readLines(tmpHtml)), tmpHtml)
    file.copy(tmpHtml, output)
    if (file.exists(tmpRmd)) file.remove(tmpRmd)
    if (file.exists(tmpHtml)) file.remove(tmpHtml)
    message("HiCool :: Report generated and available @ ", output)
    return(output)
}
