#' @title Processing Hi-C paired-end fastq files in R
#'
#' @name HiCool
#' 
#' @description 
#' 
#' `HiCool::HiCool()` automatically processes paired-end HiC sequencing files 
#' by performing the following steps: 
#' 
#' 1. Automatically setting up an appropriate conda environment using basilisk;  
#' 2. Mapping the reads to the provided genome reference using hicstuff;  
#' 3. Binning the mapped fragments into a cool file at a chosen resolution;  
#' 4. Filtering the resulting cool file to remove unwanted chromosomes (e.g. chrM);  
#' 5. Generating a multi-resolution mcool file;  
#' 6. Normalizing matrices at each resolution by iterative corretion using cooler.
#' 
#' @param r1 Path to fastq file (R1 read)
#' @param r2 Path to fastq file (R2 read)
#' @param genome Genome used to map the reads on, provided either 
#'   as a fasta file (in which case the bowtie2 index will be automatically 
#'   generated), or as a prefix to a bowtie2 index (e.g. `mm10` for 
#'   `mm10.*.bt2` files)
#' @param resolutions Resolutions used to bin the final mcool file 
#'   (Default: "1000,2000,4000,8000,16000")
#' @param restriction Restriction enzyme(s) used in HiC (Default: "DpnII,HinfI")
#' @param iterative Should the read mapping be performed iteratively? 
#'   (Default: TRUE)
#' @param filter Should read pairs be filtered (using filtering approach 
#'   described in Cournac et al., BMC Genomics 2012)? (Default: TRUE)
#' @param threads Number of CPUs used for parallelization. (Default: 16)
#' @param exclude_chr Chromosomes excluded from the final .mcool file. This will 
#'   not affect the pairs file. (Default: "Mito|chrM|MT")
#' @param output Path to output directory where processed files will 
#'   be created. (Default: `./HiCool`)
#' @param keep_bam Should the bam files be kept? (Default: FALSE)
#' @param build_report Should an automated report be computed? (Default: FALSE)
#' @param scratch Path to temporary directory where processing will take place. 
#'   (Default: `tempdir()`)
#' 
#' @return A `CoolFile` object with prefilled `pairsFile` and `metadata` slots.
#' 
#' @import HiCExperiment
#' @importClassesFrom HiCExperiment CoolFile
#' @importFrom basilisk basiliskStart
#' @importFrom basilisk basiliskStop
#' @importFrom basilisk basiliskRun
#' @export
#' 
#' @examples 
#' ## -------- Get fasta file from AnnotationHub
#' ah <- AnnotationHub::AnnotationHub()
#' twobit <- ah[['AH106801']]
#' on.exit(unlink('seq.fa'))
#' BiocIO::export(BiocIO::import(twobit), con = 'seq.fa', format = 'fasta')
#' 
#' ## -------- Get fastq files from HiContactsData
#' #r1 <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'r1_fastq')
#' #r2 <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'r2_fastq')
#' 
#' ## -------- Run HiCool::HiCool()
#' on.exit(unlink('./HiCool/'), add = TRUE)
#' #hcf <- HiCool(r1, r2, 'seq.fa', output = './HiCool/')
#' #hcf

#r1 = '/home/rsg/repos/HiContactsData/data/HiC_wt_yeast.R1.fq.gz'
#r2 = '/home/rsg/repos/HiContactsData/data/HiC_wt_yeast.R2.fq.gz'
#hcf <- HiCool(r1, r2, '~/genomes/S288c/S288c', output = './HiCool/')

HiCool <- function(
    r1, 
    r2, 
    genome, 
    resolutions = NULL, 
    restriction = 'DpnII,HinfI', 
    iterative = TRUE, 
    filter = TRUE, 
    threads = 1L, 
    output = 'HiCool', 
    exclude_chr = 'Mito|chrM|MT', 
    keep_bam = FALSE, 
    build_report = FALSE, 
    scratch = tempdir()  
)
{
    r1 <- .fixHOME(r1)
    r2 <- .fixHOME(r2)
    genome <- .fixHOME(genome)
    output <- .fixHOME(output)
    .checkGenome(genome)

    proc <- basilisk::basiliskStart(env_HiCool)
    on.exit(basilisk::basiliskStop(proc))
    on.exit(unlink(scratch), add = TRUE)
    hash <- basilisk::basiliskRun(
        env = env_HiCool, 
        fun = .processFastq,
        r1 = r1, 
        r2 = r2, 
        genome = genome, 
        resolutions = resolutions, 
        restriction = restriction, 
        iterative = iterative, 
        filter = filter, 
        threads = as.integer(threads), 
        output = output, 
        exclude_chr = exclude_chr, 
        keep_bam = keep_bam, 
        scratch = scratch  
    )
    hcf <- importHiCoolFolder(output, hash)
    metadata(hcf)$stats <- .getHicStats(gsub('.log$', '.stats', metadata(hcf)$log), filtered = filter, iterative = iterative)
    metadata(hcf)$args <- list(
        r1 = r1,
        r2 = r2,
        genome = genome,
        resolutions = resolutions,
        restriction = restriction,
        iterative = iterative,
        filter = filter,
        threads = threads,
        output = output,
        exclude_chr = exclude_chr,
        keep_bam = keep_bam,
        scratch = scratch
    )
    message("HiCool :: .fastq to .mcool processing done!")
    message("HiCool :: Check ", output, "folder to find the generated files")

    ## -- Create HiCool report 
    if (build_report) {
        message("HiCool :: Creating HiCool report...")
        HiCReport(hcf)
        message("HiCool :: All processing successfully achieved. Congrats!")
    }
    else {
        message("HiCool :: Run `HiCool::HiCReport(x)` to generate a report for this CoolFile.")
    }

    return(hcf)
}

.processFastq <- function(
    r1, 
    r2, 
    genome, 
    resolutions, 
    restriction, 
    iterative, 
    filter, 
    threads, 
    output, 
    exclude_chr, 
    keep_bam, 
    scratch  
) {

    ## -------- Import python libraries
    hs <- reticulate::import("hicstuff")
    cooler <- reticulate::import("cooler")

    ## -------- Define variables
    hash <- paste0(sample(c(LETTERS, 0:9), 6, replace = TRUE), collapse = '')
    tmp_folder <- file.path(scratch, hash)
    message("HiCool :: Initiating processing of fastq files [tmp folder: ", tmp_folder, "]...")
    prefix <- paste0(
        gsub('[._][rR][12].*', '', basename(r1)), 
        '^mapped-', gsub('.fa$', '', basename(genome)), 
        '^', hash
    )
    frags <- file.path(tmp_folder, paste0(prefix, '.frags.tsv'))
    chroms <- file.path(tmp_folder, paste0(prefix, '.chr.tsv'))
    filtered_chroms <- file.path(tmp_folder, paste0(prefix, '.chr_filtered.tsv'))
    contact_map <- file.path(tmp_folder, paste0(prefix, '.cool'))
    rebinned_prefix <- file.path(tmp_folder, paste0(prefix, '_res0'))
    contact_map_rebinned <- file.path(tmp_folder, paste0(prefix, '_res0.cool'))
    contact_map_filtered <- file.path(tmp_folder, paste0(prefix, '_res0_filtered.cool'))
    contact_map_mcool <- file.path(tmp_folder, paste0(prefix, '_res0.mcool'))
    sinked_log <- file.path(tmp_folder, paste0(prefix, '.Rlog'))

    ## -------- Save R console output to sinked_log
    dir.create(tmp_folder, showWarnings = FALSE, recursive = TRUE)
    
    ## -------- Map reads with hicstuff
    message("HiCool :: Mapping fastq files...")
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
    log_file <- list.files(tmp_folder, pattern = paste0(hash, '.hicstuff_'), full.names = TRUE)
    log <- readLines(log_file)
    log[length(log)+1] <- paste0('#### ::: WD ::: ', getwd()) 

    ## -------- Automatically deduce appropriate resolutions if unspecified
    chrs <- read.delim(file.path(tmp_folder, paste0(prefix, '.chr.tsv')), sep = '\t') 
    if (is.null(resolutions)) { 
        tot_length <- sum(chrs$length)
        resolutions_idx <- dplyr::case_when(
            tot_length < 16000 ~ 1,
            tot_length >= 16000 & tot_length < 100000000 ~ 2,
            tot_length >= 100000000 & tot_length < 1000000000 ~ 3,
            tot_length >= 1000000000 ~ 4
        )
        list_resolutions <- list(
            c(100, 200, 400, 800, 1600),
            c(1000, 2000, 4000, 8000, 16000),
            c(4000, 8000, 16000, 32000, 64000, 128000, 256000, 512000),
            c(10000, 20000, 40000, 80000, 160000, 320000, 640000, 1280000, 2560000)
        )
        first_res <- list_resolutions[[resolutions_idx]][1]
        message("HiCool :: Best-suited resolution automatically inferred: ", first_res)
        resolutions <- paste(list_resolutions[[resolutions_idx]], collapse = ',')
    }
    else {
        first_res <- as.integer(gsub(',.*', '', resolutions))
    }

    ## -------- Rebin with hicstuff rebin 
    message("HiCool :: Binning chimeric fragments...")
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
    excludable_chrs <- grep(exclude_chr, chrs$contig, value = TRUE)
    if (length(excludable_chrs)) {
        message("HiCool :: Remove unwanted chromosomes...")
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
    }
    else {
        file.copy(contact_map_rebinned, contact_map_filtered)
    }

    ## -------- Generate a multi-resolution mcool file
    message("HiCool :: Generating multi-resolution .mcool file...")
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
    message("HiCool :: Balancing .mcool file...")
    cooler$cli$zoomify$invoke_balance(
        args = paste0("--nproc ", threads, " --cis-only --min-nnz 3 --mad-max 7"), 
        resolutions = strsplit(resolutions, ',')[[1]] |> as.integer(), 
        outfile = contact_map_mcool  
    ) |> reticulate::py_capture_output() |> write(sinked_log, append = TRUE)

    ## -------- Tidy-up everything
    message("HiCool :: Tidying up everything for you :)...")
    # Matrices
    dir.create(file.path(output, 'matrices'), showWarnings = FALSE, recursive = TRUE)
    file.copy(
        contact_map_mcool, 
        file.path(output, 'matrices', paste0(prefix, '.mcool'))
    )

    # Bam
    if (keep_bam) {
        dir.create(file.path(output, 'bam'), showWarnings = FALSE, recursive = TRUE)
        file.copy(
            file.path(tmp_folder, 'tmp', paste0(prefix, '.for.bam')), 
            file.path(output, 'bam', paste0(prefix, '.fwd.bam'))
        )
        file.copy(
            file.path(tmp_folder, 'tmp', paste0(prefix, '.rev.bam')), 
            file.path(output, 'bam', paste0(prefix, '.rev.bam'))
        )
    }

    # Pairs
    dir.create(file.path(output, 'pairs'), showWarnings = FALSE, recursive = TRUE)
    file.copy(
        list.files(tmp_folder, pattern = paste0(hash, '.valid_idx_pcrfree.pairs'), recursive = TRUE, full.names = TRUE), 
        file.path(output, 'pairs', paste0(prefix, '.pairs'))
    )

    # Plots
    dir.create(file.path(output, 'plots'), showWarnings = FALSE, recursive = TRUE)
    file.copy(
        list.files(tmp_folder, pattern = paste0(hash, '_event_distance.pdf'), recursive = TRUE, full.names = TRUE), 
        file.path(output, 'plots', paste0(prefix, '_event_distance.pdf'))
    )
    file.copy(
        list.files(tmp_folder, pattern = paste0(hash, '_event_distribution.pdf'), recursive = TRUE, full.names = TRUE), 
        file.path(output, 'plots', paste0(prefix, '_event_distribution.pdf'))
    )

    # Log
    dir.create(file.path(output, 'logs'), showWarnings = FALSE, recursive = TRUE)
    file.copy(
        log_file, 
        file.path(output, 'logs', paste0(prefix, '.log'))
    )
    writeLines(
        log,
        file.path(output, 'logs', paste0(prefix, '.stats'))
    )
    return(hash)
}
