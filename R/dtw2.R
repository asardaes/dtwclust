dtw2 <- function(x, y, ...) {
     lcm <- proxy::dist(x, y, method = "L1")

     sqrt(dtw::dtw(lcm^2, distance.only = TRUE, ...)$distance)
}
