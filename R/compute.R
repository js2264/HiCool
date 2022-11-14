#' @title Process Hi-C raw paired-end fastq files into a (m)cool processed data file.
#' @rdname HiCool
#'
#' @importFrom basilisk basiliskStart
#' @importFrom basilisk basiliskStop
#' @importFrom basilisk basiliskRun
#' @export

HiCool <- function(
    r1, 
    r2, 
    genome, 
    resolutions = '1000,2000,4000,8000,16000', 
    restriction = 'DpnII,HinfI', 
    iterative = TRUE, 
    filter = TRUE, 
    threads = 16L, 
    output_folder = 'HiCool', 
    exclude_chr = 'Mito|chrM|MT', 
    scratch = 'tmp'  
)
{
    proc <- basilisk::basiliskStart(env_HiCool)
    on.exit(basilisk::basiliskStop(proc))
    hash <- basilisk::basiliskRun(
        env = env_HiCool, 
        fun = .HiCool,
        r1 = .fix_HOME(r1), 
        r2 = .fix_HOME(r2), 
        genome = .fix_HOME(genome), 
        resolutions = resolutions, 
        restriction = restriction, 
        iterative = iterative, 
        filter = filter, 
        threads = threads, 
        output_folder = .fix_HOME(output_folder), 
        exclude_chr = exclude_chr, 
        scratch = scratch  
    )
    files <- list.files(output_folder, pattern = hash, full.names = TRUE, recursive = TRUE)
    log_file <- grep('log', files, value = TRUE)
    res <- HiCExperiment::HiCoolFile(
        mcool = grep('mcool', files, value = TRUE), 
        resolution = as.integer(gsub(',.*', '', resolutions)),
        pairs = grep('pairs', files, value = TRUE), 
        log = list(

        )
    )
    return(res)
}

#' @importFrom reticulate import
#' 

.HiCool <- function(
    r1, 
    r2, 
    genome, 
    resolutions, 
    restriction, 
    iterative, 
    filter, 
    threads, 
    output_folder, 
    exclude_chr, 
    scratch  
) {

    ## -------- Import python libraries
    hs <- reticulate::import("hicstuff")
    cooler <- reticulate::import("cooler")

    ## -------- Define variables
    hash <- paste0(sample(c(LETTERS, 0:9), 6, replace = TRUE), collapse = '')
    tmp_folder <- file.path(scratch, hash)
    prefix <- paste0(
        gsub('_[rR].*', '', basename(r1)), 
        '^', basename(genome), 
        '^', hash
    )
    first_res <- as.integer(gsub(',.*', '', resolutions))
    frags <- file.path(tmp_folder, paste0(prefix, '.frags.tsv'))
    chroms <- file.path(tmp_folder, paste0(prefix, '.chr.tsv'))
    filtered_chroms <- file.path(tmp_folder, paste0(prefix, '.chr_filtered.tsv'))
    contact_map <- file.path(tmp_folder, paste0(prefix, '.cool'))
    rebinned_prefix <- file.path(tmp_folder, paste0(prefix, '_', first_res))
    contact_map_rebinned <- file.path(tmp_folder, paste0(prefix, '_', first_res, '.cool'))
    contact_map_filtered <- file.path(tmp_folder, paste0(prefix, '_', first_res, '_filtered.cool'))
    contact_map_mcool <- file.path(tmp_folder, paste0(prefix, '_', first_res, '.mcool'))
    sinked_log <- file.path(tmp_folder, paste0(prefix, '.Rlog'))

    ## -------- Save R console output to sinked_log
    dir.create(tmp_folder, showWarnings = FALSE, recursive = TRUE)

    ## -------- Map reads with hicstuff
    hs$pipeline$full_pipeline(
        input1 = r1, 
        input2 = r2, 
        genome = genome, 
        enzyme = restriction, 
        filter_events = filter, 
        force = TRUE, 
        mapping = ifelse(iterative, "iterative", "normal"),
        mat_fmt = "cool",
        # mat_fmt = "graal",
        no_cleanup = TRUE,
        out_dir = tmp_folder, 
        pcr_duplicates = TRUE, 
        plot = TRUE, 
        prefix = prefix, 
        threads = threads,
        distance_law = TRUE
    ) |> reticulate::py_capture_output() |> write(sinked_log)

    ## -------- Rebin with hicstuff rebin 
    hs$commands$Rebin(
        command_args = paste0(
            " --binning ", paste0(first_res/1000, 'kb'), 
            " --frags ", frags, 
            " --chroms ", chroms,
            " --force ", 
            contact_map, ' ', rebinned_prefix 
        ),
        global_args = ""
    )$execute() |> reticulate::py_capture_output() |> write(sinked_log, append = TRUE)
    
    ## -------- Exclude unwanted chr. from cool file
    chr <- readLines(file.path(tmp_folder, paste0(prefix, '.chr.tsv'))) 
    chr <- grep(exclude_chr, chr, invert = TRUE, value = TRUE)
    chr <- grep('contig', chr, invert = TRUE, value = TRUE)
    writeLines(chr, filtered_chroms)
    cooler$cli$dump$dump$callback(
        contact_map_rebinned, 
        table = "pixels", 
        columns = NULL, 
        header = FALSE, 
        na_rep = "", 
        float_format = 'g', 
        range = NULL, 
        range2 = NULL, 
        matrix = FALSE, 
        balanced = FALSE, 
        join = TRUE, 
        annotate = NULL, 
        one_based_ids = FALSE, 
        one_based_starts = FALSE, 
        chunksize = NULL, 
        out = paste0(contact_map_filtered, '_tmp')
    ) |> reticulate::py_capture_output() |> write(sinked_log, append = TRUE)
    cooler$cli$load$load$callback(
        bins_path = paste0(
            filtered_chroms, ":", first_res
        ), 
        pixels_path = paste0(contact_map_filtered, '_tmp'), 
        cool_path = contact_map_filtered, 
        format = 'bg2', 
        metadata = NULL, 
        assembly = NULL, 
        chunksize = 20e6L, 
        field = "", 
        count_as_float = FALSE, 
        one_based = FALSE, 
        comment_char = "#", 
        input_copy_status = NULL, 
        no_symmetric_upper = FALSE, 
        storage_options = NULL
    ) |> reticulate::py_capture_output() |> write(sinked_log, append = TRUE)

    ## -------- Generate a multi-resolution mcool file
    cooler$zoomify_cooler(
        base_uris = contact_map_filtered, 
        outfile = contact_map_mcool, 
        resolutions = strsplit(resolutions, ',')[[1]] |> as.integer(), 
        chunksize = 10000000L, 
        nproc = threads, 
        columns = NULL, 
        dtypes = NULL, 
        agg = NULL
    ) |> reticulate::py_capture_output() |> write(sinked_log, append = TRUE)
    cooler$cli$zoomify$invoke_balance(
        args = paste0("--nproc ", threads, " --cis-only --min-nnz 3 --mad-max 7"), 
        resolutions = strsplit(resolutions, ',')[[1]] |> as.integer(), 
        outfile = contact_map_mcool  
    ) |> reticulate::py_capture_output() |> write(sinked_log, append = TRUE)

    ## -------- Tidy-up everything
    # Mcool
    dir.create(file.path(output_folder, 'cool'), showWarnings = FALSE, recursive = TRUE)
    file.copy(
        contact_map_mcool, 
        file.path(output_folder, 'cool', paste0(prefix, '.mcool'))
    )

    # Bam
    dir.create(file.path(output_folder, 'bam'), showWarnings = FALSE, recursive = TRUE)
    file.copy(
        file.path(tmp_folder, 'tmp', paste0(prefix, '.for.bam')), 
        file.path(output_folder, 'bam', paste0(prefix, '.fwd.bam'))
    )
    file.copy(
        file.path(tmp_folder, 'tmp', paste0(prefix, '.rev.bam')), 
        file.path(output_folder, 'bam', paste0(prefix, '.rev.bam'))
    )

    # Pairs
    dir.create(file.path(output_folder, 'pairs'), showWarnings = FALSE, recursive = TRUE)
    file.copy(
        list.files(tmp_folder, pattern = paste0(hash, '.valid_idx_pcrfree.pairs'), recursive = TRUE, full.names = TRUE), 
        file.path(output_folder, 'pairs', paste0(prefix, '.pairs'))
    )

    # Pairs
    dir.create(file.path(output_folder, 'pairs'), showWarnings = FALSE, recursive = TRUE)
    file.copy(
        list.files(tmp_folder, pattern = paste0(hash, '.valid_idx_pcrfree.pairs'), recursive = TRUE, full.names = TRUE), 
        file.path(output_folder, 'pairs', paste0(prefix, '.pairs'))
    )

    # Plots
    dir.create(file.path(output_folder, 'plots'), showWarnings = FALSE, recursive = TRUE)
    file.copy(
        list.files(tmp_folder, pattern = paste0(hash, '_event_distance.pdf'), recursive = TRUE, full.names = TRUE), 
        file.path(output_folder, 'plots', paste0(prefix, '_event_distance.pdf'))
    )
    file.copy(
        list.files(tmp_folder, pattern = paste0(hash, '_event_distribution.pdf'), recursive = TRUE, full.names = TRUE), 
        file.path(output_folder, 'plots', paste0(prefix, '_event_distribution.pdf'))
    )

    # Log
    dir.create(file.path(output_folder, 'logs'), showWarnings = FALSE, recursive = TRUE)
    file.copy(
        list.files(tmp_folder, pattern = paste0(hash, '.hicstuff_'), full.names = TRUE), 
        file.path(output_folder, 'logs', paste0(prefix, '.log'))
    )

    return(hash)
}
