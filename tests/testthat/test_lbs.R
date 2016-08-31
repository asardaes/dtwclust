context("Test lower bounds consistency")

# =================================================================================================
# L1 norm
# =================================================================================================

lbk <- lb_keogh(data_matrix[2L, ], data_matrix[1L, ], window.size = 15L)$d
lbi <- lb_improved(data_matrix[2L, ], data_matrix[1L, ], window.size = 15L)
dtwd <- dtw(data_matrix[2L, ], data_matrix[1L, ],
            window.type = "sakoechiba", window.size = 15L,
            distance.only = TRUE)$distance

lbks <- proxy::dist(data_list[-1L], data_list[1L],
                    method = "lbk", window.size = 15L)
lbis <- proxy::dist(data_list[-1L], data_list[1L],
                    method = "lbi", window.size = 15L)
dtwds <- proxy::dist(data_list[-1L], data_list[1L],
                     method = "dtw", window.type = "sakoechiba",
                     window.size = 15L)

test_that("LBK is leq LBI with L1 norm",
          expect_true(lbk <= lbi))

test_that("LBK is leq DTW with L1 norm",
          expect_true(lbk <= dtwd))

test_that("LBI is leq DTW with L1 norm",
          expect_true(lbi <= dtwd))

test_that("multiple LBK are leq LBI with L1 norm",
          expect_true(all(lbks <= lbis)))

test_that("multiple LBK are leq DTW with L1 norm",
          expect_true(all(lbks <= dtwds)))

test_that("multiple LBI are leq DTW with L1 norm",
          expect_true(all(lbis <= dtwds)))

test_that("single LBK_1 matches proxy version",
          expect_true(lbk == lbks[1L, 1L]))

test_that("single LBI_1 matches proxy version",
          expect_true(lbi == lbis[1L, 1L]))

test_that("single DTW_1 matches proxy version",
          expect_true(dtwd == dtwds[1L, 1L]))

# =================================================================================================
# L2 norm
# =================================================================================================

lbk <- lb_keogh(data_matrix[2L, ], data_matrix[1L, ],
                window.size = 15L, norm = "L2")$d
lbi <- lb_improved(data_matrix[2L, ], data_matrix[1L, ],
                   window.size = 15L, norm = "L2")
dtwd <- proxy::dist(data_matrix[2L, ], data_matrix[1L, ], method = "L1")
dtwd <- sqrt(dtw(dtwd^2,
                 window.type = "sakoechiba", window.size = 15L,
                 distance.only = TRUE)$distance)

lbks <- proxy::dist(data_list[-1L], data_list[1L],
                    method = "lbk", window.size = 15L,
                    norm = "L2")
lbis <- proxy::dist(data_list[-1L], data_list[1L],
                    method = "lbi", window.size = 15L,
                    norm = "L2")
dtwds <- proxy::dist(data_list[-1L], data_list[1L],
                     method = "dtw2", window.type = "sakoechiba",
                     window.size = 15L)

test_that("LBK is leq LBI with L2 norm",
          expect_true(lbk <= lbi))

test_that("LBK is leq DTW with L2 norm",
          expect_true(lbk <= dtwd))

test_that("LBI is leq DTW with L2 norm",
          expect_true(lbi <= dtwd))

test_that("multiple LBK are leq LBI with L2 norm",
          expect_true(all(lbks <= lbis)))

test_that("multiple LBK are leq DTW with L2 norm",
          expect_true(all(lbks <= dtwds)))

test_that("multiple LBI are leq DTW with L1 norm",
          expect_true(all(lbis <= dtwds)))

test_that("single LBK_2 matches proxy version",
          expect_true(lbk == lbks[1L, 1L]))

test_that("single LBI_2 matches proxy version",
          expect_true(lbi == lbis[1L, 1L]))

test_that("single DTW_2 matches proxy version",
          expect_true(dtwd == dtwds[1L, 1L]))
