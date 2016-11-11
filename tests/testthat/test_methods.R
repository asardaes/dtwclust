context("Test methods")

test_that("Methods for dtwclust objects are dispatched correctly.", {
    pc <- dtwclust(data_subset[-1L], type = "partitional", k = 4L,
                   distance = "sbd", preproc = zscore)

    expect_error(plot(pc, type = "dendrogram"), "dendrogram.*hierarchical", ignore.case = TRUE)

    ## extra argument for preproc so that it is used in predict
    hc <- dtwclust(data_subset[-1L], type = "hierarchical", k = 4L,
                   distance = "sbd", preproc = zscore, center = FALSE)

    mute <- capture.output(show_return <- show(hc))
    expect_null(show_return, info = "show")

    expect_null(plot(hc), info = "plot dendrogram")
    expect_true(inherits(plot(hc, type = "sc", plot = FALSE), "ggplot"), info = "plot sc")
    expect_true(inherits(plot(hc, type = "sc", data = data_subset[-1L], plot = FALSE), "ggplot"),
                info = "plot sc with data")

    expect_true(predict(hc)[1L] == predict(hc, newdata = data_subset[1L]), info = "predict")

    suppressMessages(method_update <- update(hc))
    expect_identical(hc, method_update)
    method_update <- update(hc, method = "complete", distmat = hc@distmat)
    expect_true(inherits(method_update, "dtwclust") & validObject(method_update), info = "update")

    require(flexclust)
    clusterSim(hc)
    clusterSim(hc, symmetric = TRUE)
    clusterSim(hc, method = "centers")

    skip_on_cran()

    method_update <- reset_nondeterministic(method_update)
    expect_equal_to_reference(method_update, file_name(method_update, x32 = TRUE))
})
