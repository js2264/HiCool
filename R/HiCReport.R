#' @title HiC processing report
#' 
#' @name HiCReport
#' 
#' @param x an CoolFile object, generated from `HiCool::HiCool()` or 
#'   `HiCool::importHiCoolFolder()`, or directly from calling 
#'   `HiCExperiment::CoolFile()`.
#' @param output Path to save output HTML file.
#'
#' @importFrom usethis use_template
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
    hash <- NULL
    args <- NULL
    log_file <- S4Vectors::metadata(x)$log
    hash <- S4Vectors::metadata(x)$hash
    args <- S4Vectors::metadata(x)$args
    
    ## --- Check that all required fields exist
    if (!file.exists(mcool))
        stop("The associated cool file is missing.")
    if (!file.exists(pairs))
        stop("The associated pairs file is missing.")
    if (is.null(log_file))
        stop("Missing log file.")
    if (!file.exists(log_file))
        stop("Log file does not exist.")
    if (is.null(hash))
        stop("Missing `args` metadata.")
    if (is.null(args))
        stop("Missing `args` metadata.")
    
    ## --- Prepare output files
    if (is.null(output)) {
        output <- file.path(
            dirname(dirname(mcool)), 
            gsub('\\..?cool', '.html', basename(mcool))
        )
    }
    tmpRmd <- gsub('.html$', '.Rmd', output)

    ## --- Prepare data
    data <- list(
        x = x, 
        wd = getwd()
    )

    ## --- Generate report
    unlink(tmpRmd)
    usethis::use_template('HiCoolFile_report_template.Rmd', package = 'HiCool', save_as = tmpRmd, data = data)
    rmarkdown::render(tmpRmd, quiet = TRUE)
    writeLines(gsub("&amp;quot;", '"', readLines(output)), output)
    unlink(tmpRmd)
}
