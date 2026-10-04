test_that("Check that utils work", {

    log_path <- HiContactsData::HiContactsData(sample = 'yeast_wt', format = 'HiCool_log')
    expect_type(getHicStats(log_path), 'list')
    expect_type(getHiCoolArgs(log_path), 'list')

})

test_that("getHicStats() parses logs of hicstuff < 3.2.5 and >= 3.2.5", {

    filtering <- c(
        "2026-10-04,13:53:02 :: INFO :: Filtering with thresholds: uncuts=7 loops=7",
        "2026-10-04,13:53:02 :: INFO :: 11208 pairs discarded: Loops: 1910, Uncuts: 9266, Weirds: 32",
        "2026-10-04,13:53:02 :: INFO :: 53553 pairs kept (82.69%)"
    )
    log_3.2.4 <- c(
        "2026-10-04,13:53:00 :: INFO :: 76% reads (single ends) mapped with Q >= 30 (151804/200000)",
        filtering,
        "2026-10-04,13:53:05 :: INFO :: 1% PCR duplicates have been filtered out (613 / 53553 pairs) "
    )
    log_3.2.5 <- c(
        "2026-10-04,13:53:00 :: INFO :: 76.00% reads (single ends) mapped with Q >= 30 (151804/200000)",
        filtering,
        "2026-10-04,13:53:05 :: INFO :: 1.1% PCR duplicates have been filtered out (613/53553 pairs)"
    )
    expected <- list(
        nFragments = 100000, nPairs = 64761, nDangling = 9266, nSelf = 1910,
        nDumped = 32, nFiltered = 53553, nDups = 613, nUnique = 52940,
        threshold_uncut = 7, threshold_self = 7
    )
    for (lines in list(log_3.2.4, log_3.2.5)) {
        log <- tempfile(fileext = '.log')
        writeLines(lines, log)
        expect_equal(getHicStats(log), expected)
    }

    ## -- Without filtering of events
    log <- tempfile(fileext = '.log')
    writeLines(log_3.2.5[-(2:4)], log)
    stats <- getHicStats(log)
    expect_equal(stats$nPairs, 53553)
    expect_equal(stats$nUnique, 52940)
    expect_true(is.na(stats$threshold_uncut))

})

test_that("chromosight arguments are turned into command-line arguments", {

    args <- list(
        "--pattern" = "loops",
        "--inter" = FALSE,
        "--kernel-config" = NULL,
        "--no-plotting" = TRUE,
        "<contact_map>" = "x.mcool::/resolutions/1000",
        "--max-dist" = 2e6,
        "--n-mads" = 5L,
        "<prefix>" = "out/chromo"
    )
    expect_identical(.docoptArgs(args), c(
        "--no-plotting", "--pattern=loops", "--max-dist=2000000", "--n-mads=5",
        "x.mcool::/resolutions/1000", "out/chromo"
    ))

})
