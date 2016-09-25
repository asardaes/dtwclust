context("Test dtw_basic consistency")

test_that("dtw_basic gives the same results as dtw/dtw2", {
     for (rep in 1L:100L) {
          i <- sample(length(data), 1L)
          j <- sample(length(data), 1L)
          norm <- sample(c("L1", "L2"), 1L)

          if (i <= 50L)
               window.size <- 20L
          else
               window.size <- NULL

          if (j <= 50L)
               step.pattern <- symmetric1
          else
               step.pattern <- symmetric2

          if (is.null(window.size))
               window.type <- "none"
          else
               window.type <- "slantedband"

          x <- data[[i]]
          y <- data[[j]]

          if (norm == "L1")
               d1 <- dtw(x, y, step.pattern = step.pattern,
                         window.type = window.type, window.size = window.size,
                         keep.internals = TRUE)
          else
               d1 <- dtw2(x, y, step.pattern = step.pattern,
                          window.type = window.type, window.size = window.size,
                          keep.internals = TRUE)

          storage.mode(d1$index1) <- "integer"
          storage.mode(d1$index2) <- "integer"

          d2 <- dtw_basic(x, y, window.size = window.size, norm = norm,
                          step.pattern = step.pattern, backtrack = TRUE)

          expect_identical(d1$distance, d2$distance, info = paste("Indices:", i, j))
          expect_identical(d1$index1, d2$index1, info = paste("Indices:", i, j))
          expect_identical(d1$index2, d2$index2, info = paste("Indices:", i, j))

          if (identical(step.pattern, symmetric2)) {
               d3 <- dtw_basic(x, y, window.size = window.size, norm = norm,
                               step.pattern = step.pattern, normalize = TRUE)

               expect_identical(d1$normalizedDistance, d3, info = paste("Indices:", i, j))
          }
     }

     D1 <- proxy::dist(data[1L:16L], data[1L:16L], method = "dtw")
     D2 <- proxy::dist(data[1L:16L], data[1L:16L], method = "dtw_basic")

     expect_equal(D1, D2, tolerance = 0, check.attributes = FALSE,
                  info = "dtw vs dtw_basic")

     D1 <- proxy::dist(data[1L:16L], data[1L:16L], method = "dtw2")
     D2 <- proxy::dist(data[1L:16L], data[1L:16L], method = "dtw_basic", norm = "L2")

     expect_equal(D1, D2, tolerance = 0, check.attributes = FALSE,
                  info = "dtw2 vs dtw_basic")
})
