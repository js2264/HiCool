test_that("Check that reports work", {

    mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
    pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
    log_path <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'HiCool_log')
    cf <- CoolFile(mcool_path, pairs = pairs_path, metadata = list(log = log_path))
    expect_no_error(HiCReport(cf))
    
})
