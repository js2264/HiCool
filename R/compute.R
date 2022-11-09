#' @title Process Hi-C raw paired-end fastq files into a (m)cool processed data file.
#' @rdname HiCompute
#'
#' @importFrom basilisk basiliskStart
#' @importFrom basilisk basiliskStop
#' @importFrom basilisk basiliskRun
#' @export

HiCompute <- function(
    r1, 
    r2, 
    genome, 
    resolutions = '1000,2000,4000,8000,16000', 
    restriction = 'DpnII,HinfI', 
    iterative = TRUE, 
    filter = TRUE, 
    threads = 16L, 
    output_folder = 'HiCompute/', 
    exclude_chr = 'Mito|chrM|MT'  
)
{
    proc <- basilisk::basiliskStart(env_HiCompute)
    on.exit(basilisk::basiliskStop(proc))
    tmp_folder <- paste0('/tmp/Rtmp', paste0(sample(c(letters, 0:9), 6, replace = TRUE), collapse = ''), '/')
    basilisk::basiliskRun(
        env = env_HiCompute, 
        fun = .HiCompute,
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
        tmp_folder = tmp_folder  
    )
    return(TRUE)
}

#' @importFrom reticulate import
#' 

.HiCompute <- function(
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
    tmp_folder  
) {

    ## -------- Import python libraries
    hs <- reticulate::import("hicstuff")
    cooler <- reticulate::import("cooler")

    ## -------- Define variables
    hash <- paste0(sample(c(LETTERS, 0:9), 6, replace = TRUE), collapse = '')
    prefix <- paste0(gsub('_[rR].*', '', basename(r1)), '^', hash)
    first_res <- as.integer(gsub(',.*', '', resolutions))
    frags <- paste0(tmp_folder, prefix, '.frags.tsv')
    chroms <- paste0(tmp_folder, prefix, '.chr.tsv')
    filtered_chroms <- paste0(tmp_folder, prefix, '.chr_filtered.tsv')
    contact_map <- paste0(tmp_folder, prefix, '.cool')
    out_prefix <- paste0(tmp_folder, prefix, '_', first_res)
    contact_map_rebinned <- paste0(out_prefix, '.cool')
    contact_map_filtered <- paste0(out_prefix, '_filtered.cool')
    mcool <- paste0(out_prefix, '.mcool')

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
    )

    ## -------- Rebin with hicstuff rebin 
    hs$commands$Rebin(
        command_args = paste0(
            " --binning ", paste0(first_res/1000, 'kb'), 
            " --frags ", frags, 
            " --chroms ", chroms,
            " --force ", 
            contact_map, ' ', out_prefix 
        ),
        global_args = ""
    )$execute()
    
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
    )
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
    )

    ## -------- Generate a multi-resolution mcool file
    cooler$zoomify_cooler(
        base_uris = contact_map_filtered, 
        outfile = mcool, 
        resolutions = strsplit(resolutions, ',')[[1]] |> as.integer(), 
        chunksize = 10000000L, 
        nproc = threads, 
        columns = NULL, 
        dtypes = NULL, 
        agg = NULL
    )
    cooler$cli$zoomify$invoke_balance(
        args = "--cis-only --min-nnz 3 --mad-max 7", 
        resolutions = strsplit(resolutions, ',')[[1]] |> as.integer(), 
        outfile = mcool  
    )

    ## -------- Save commands to R script 
    dir.create(file.path(output_folder, 'logs'), showWarnings = FALSE, recursive = TRUE)

    ## -------- Tidy-up everything
    # Mcool
    dir.create(file.path(output_folder, 'cool', basename(genome)), showWarnings = FALSE, recursive = TRUE)
    file.copy(
        mcool, 
        file.path(output_folder, 'cool', basename(genome), paste0(prefix, '.mcool'))
    )

    # Bam
    dir.create(file.path(output_folder, 'bam', basename(genome)), showWarnings = FALSE, recursive = TRUE)
    file.copy(
        file.path(tmp_folder, 'tmp', paste0(prefix, '.for.bam')), 
        file.path(output_folder, 'bam', basename(genome), paste0(prefix, '.fwd', '.bam'))
    )
    file.copy(
        file.path(tmp_folder, 'tmp', paste0(prefix, '.rev.bam')), 
        file.path(output_folder, 'bam', basename(genome), paste0(prefix, '.rev', '.bam'))
    )

    # Pairs
    dir.create(file.path(output_folder, 'pairs', basename(genome)), showWarnings = FALSE, recursive = TRUE)
    file.copy(
        list.files(tmp_folder, pattern = paste0(hash, '.valid_idx_pcrfree.pairs'), recursive = TRUE, full.names = TRUE), 
        file.path(output_folder, 'pairs', basename(genome), paste0(prefix, '.rev', '.bam'))
    )

    # Log
    dir.create(file.path(output_folder, 'logs'), showWarnings = FALSE, recursive = TRUE)
    file.copy(
        list.files(tmp_folder, pattern = paste0(hash, '.hicstuff_'), full.names = TRUE), 
        file.path(output_folder, 'logs', paste0(prefix, '.log'))
    )

    return(TRUE)
}
