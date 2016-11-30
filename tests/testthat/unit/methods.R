context("\tMethods")

# =================================================================================================
# setup
# =================================================================================================

## Original objects in env
ols <- ls()

# =================================================================================================
# generics
# =================================================================================================

test_that("Methods for dtwclust objects are dispatched correctly.", {
    mute <- capture.output(pc <- dtwclust(data_subset[-1L], type = "partitional", k = 4L,
                                          distance = "sbd", preproc = zscore,
                                          control = list(trace = TRUE)))

    expect_error(plot(pc, type = "dendrogram"), "dendrogram.*hierarchical", ignore.case = TRUE)

    ## fuzzy
    mute <- capture.output(fc <- dtwclust(data_reinterpolated_subset[-1L], type = "fuzzy",
                                          k = 4L, distance = "sbd", preproc = zscore,
                                          control = list(trace = TRUE)))

    expect_error(plot(fc, type = "dendrogram"), "dendrogram.*hierarchical", ignore.case = TRUE)

    ## extra argument for preproc so that it is used in predict
    mute <- capture.output(hc <- dtwclust(data_subset[-1L], type = "hierarchical", k = 4L,
                                          distance = "sbd", preproc = zscore, center = FALSE,
                                          control = list(trace = TRUE)))

    mute <- capture.output(show_return <- show(pc))
    expect_null(show_return, info = "show")
    mute <- capture.output(show_return <- show(hc))
    expect_null(show_return, info = "show")
    mute <- capture.output(show_return <- show(fc))
    expect_null(show_return, info = "show")

    expect_null(plot(hc), info = "plot dendrogram")
    expect_true(inherits(plot(hc, type = "sc", plot = FALSE), "ggplot"), info = "plot sc")
    expect_true(inherits(plot(hc, type = "sc", data = data_subset[-1L], plot = FALSE), "ggplot"),
                info = "plot sc with data")

    expect_true(predict(hc)[1L] == predict(hc, newdata = data_subset[1L]), info = "predict")

    suppressMessages(method_update <- update(hc))
    expect_identical(hc, method_update)
    mute <- capture.output(method_update <- update(hc, method = "complete", distmat = hc@distmat))
    expect_true(inherits(method_update, "dtwclust") & validObject(method_update), info = "update")

    require(flexclust)
    clusterSim(hc)
    clusterSim(hc, symmetric = TRUE)
    clusterSim(hc, method = "centers")

    ## refs
    method_update <- reset_nondeterministic(method_update)
    assign("method_update", method_update, persistent)
})

# =================================================================================================
# clean
# =================================================================================================
rm(list = setdiff(ls(), ols))
