context("Test methods")

test_that("Methods for dtwclust objects are dispatched correctly.", {
     ## extra argument for preproc so that it is used in predict
     hc <- dtwclust(data_subset[-1L], type = "hierarchical", k = 4L,
                    distance = "sbd", preproc = zscore, center = FALSE)

     mute <- capture.output(show_return <- show(hc))
     expect_null(show_return, info = "show")
     expect_null(plot(hc), info = "plot dendrogram")

     expect_true(inherits(plot(hc, type = "sc", plot = FALSE), "ggplot"), info = "plot sc")
     expect_true(predict(hc)[1L] == predict(hc, newdata = data_subset[1L]), info = "predict")

     method_update <- update(hc, method = "complete", distmat = hc@distmat)

     expect_true(inherits(method_update, "dtwclust") & validObject(method_update), info = "update")

     skip_on_cran()

     method_update <- reset_nondeterministic(method_update)
     expect_equal_to_reference(method_update, file_name(method_update))
})
