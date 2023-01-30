test_that("Check that HiCool works with AWS S3 iGenomes", {
    
    ## -------- Get fastq file from AnnotationHub
    r1 <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'fastq_R1')
    r2 <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'fastq_R2')

    ## -------- Run HiCool::HiCool() with AWS S3 iGenomes
    ## HiCReport() is tested here as well, since it is executed by HiCool() as long as `build_report == TRUE`
    hcf <- HiCool(r1, r2, 'R64-1-1', output = './HiCool/')
    expect_s4_class(hcf, 'CoolFile')

})
