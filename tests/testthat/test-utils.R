test_that("Check that utils work", {

    log_path <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'HiCool_log')
    expect_type(getHicStats(log_path), 'list')
    expect_type(getHiCoolArgs(log_path), 'list')
    
})
