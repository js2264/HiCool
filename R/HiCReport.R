#' @title HiC processing report
#' 
#' @name HiCReport
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
#' @export
#' @examples 
#' HiCReport(hcf)

HiCReport <- function(x, output = NULL) {

    ## --- Parse CoolFile for info
    mcool <- BiocIO::resource(x)
    pairs <- pairsFile(x)
    log_file <- NULL
    args <- NULL
    log_file <- S4Vectors::metadata(x)$log
    
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
    }

    ## --- Generate report
    tmpRmd <- basename(gsub('.html$', '.Rmd', output))
    unlink(tmpRmd)
    unlink(output)
    file.copy(
        fs::path_package('templates', 'HiCoolFile_report_template.Rmd', package = 'HiCool'), 
        tmpRmd
    )
    writeLines(gsub("%COOL%", mcool, readLines(tmpRmd)), tmpRmd)
    writeLines(gsub("%PAIRS%", pairs, readLines(tmpRmd)), tmpRmd)
    writeLines(gsub("%LOG%", log_file, readLines(tmpRmd)), tmpRmd)
    rmarkdown::render(tmpRmd)
    writeLines(gsub("&amp;quot;", '"', readLines(basename(output))), basename(output))
    file.copy(basename(output), output)
    unlink(tmpRmd)
    unlink(basename(output))
    message("Report generated and available @ ", output)
}

HiCReports <- function(x, output = NULL) {

    ## --- Parse CoolFiles 
    log_files <- lapply(x, function(y) metadata(y)$log) |> 
        paste(collapse = ',')
    
    ## --- Prepare output files
    if (is.null(output)) {
        output <- 'HiCReports.html'
    }

    ## --- Generate report
    tmpRmd <- basename(gsub('.html$', '.Rmd', output))
    unlink(tmpRmd)
    unlink(output)
    file.copy(
        fs::path_package('templates', 'HiCoolFile_reports_template.Rmd', package = 'HiCool'), 
        tmpRmd
    )
    writeLines(gsub("%LOGS%", log_files, readLines(tmpRmd)), tmpRmd)
    rmarkdown::render(tmpRmd)
    writeLines(gsub("&amp;quot;", '"', readLines(basename(output))), basename(output))
    file.copy(basename(output), output)
    unlink(tmpRmd)
    message("Report generated and available @ ", output)
}
