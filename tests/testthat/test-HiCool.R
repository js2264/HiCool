test_that("Check that HiCool works with fasta or AWS S3 iGenomes", {
    
    ## -------- Get fastq/fasta file from AnnotationHub
    r1 <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'fastq_R1')
    r2 <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'fastq_R2')
    # ah <- AnnotationHub::AnnotationHub()
    # twobit <- ah[['AH106801']]
    # tmpfile <- tempfile()
    # BiocIO::export(BiocIO::import(twobit), con = tmpfile, format = 'fasta')

    ## -------- Run HiCool::HiCool() with fasta
    # hcf <- HiCool(r1, r2, tmpfile, output = './HiCool/')
    # expect_s4_class(hcf, 'CoolFile')

    ## -------- Run HiCool::HiCool() with AWS S3 iGenomes
    hcf <- HiCool(r1, r2, 'R64-1-1', output = './HiCool/')
    expect_s4_class(hcf, 'CoolFile')

})
