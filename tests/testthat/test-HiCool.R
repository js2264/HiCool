test_that("Check that packages are installed", {
    ## -------- Get fasta file from AnnotationHub
    # ah <- AnnotationHub::AnnotationHub()
    # twobit <- ah[['AH106801']]
    # on.exit(unlink('seq.fa'))
    # BiocIO::export(BiocIO::import(twobit), con = 'seq.fa', format = 'fasta')
    #' 
    ## -------- Get fastq files from HiContactsData
    # r1 <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'r1_fastq')
    # r2 <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'r2_fastq')
    #' 
    ## -------- Run HiCool::HiCool()
    # on.exit(unlink('./HiCool/'), add = TRUE)
    # hcf <- HiCool(r1, r2, 'seq.fa', output_folder = './HiCool/')
    # expect_s4_class(cf, 'CoolFile')

    cool_path <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'cool')
    cf <- HiCExperiment::CoolFile(cool_path)
    expect_s4_class(cf, 'CoolFile')
})
