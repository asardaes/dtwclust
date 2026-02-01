# all.equal doesn't work with dist objects anymore because they're considered numeric vectors but their names are lists because they're conceptually 2-dimensional
all.equal.dist <- function(target, current, ...) {
    NextMethod(check.names = FALSE)
}

registerS3method("all.equal", "dist", all.equal.dist, .GlobalEnv)
