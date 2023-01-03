test_that("Check that HiCool works with fasta", {
    
    ## -------- Get fasta file from AnnotationHub
    ah <- AnnotationHub::AnnotationHub()
    twobit <- ah[['AH106801']]
    on.exit(unlink('seq.fa'))
    BiocIO::export(BiocIO::import(twobit), con = 'seq.fa', format = 'fasta')

    ## -------- Get fastq files from HiContactsData
    r1 <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'fastq_R1')
    r2 <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'fastq_R2')

    ## -------- Run HiCool::HiCool()
    on.exit(unlink('./HiCool/'), add = TRUE)
    hcf <- HiCool(r1, r2, 'seq.fa', output = './HiCool/')
    expect_s4_class(hcf, 'CoolFile')

})

test_that("Check that HiCool works with AWS S3 iGenomes", {
    
    r1 <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'fastq_R1')
    r2 <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'fastq_R2')

    on.exit(unlink('./HiCool/'), add = TRUE)
    hcf <- HiCool(r1, r2, 'R64-1-1', output = './HiCool/')
    expect_s4_class(hcf, 'CoolFile')

})
