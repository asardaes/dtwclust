# Explore tab, characteristics
characteristics_table <- quote({
    data.frame(
        "Amount of series" = length(.series_),
        "Multivariate" = is_multivariate(.series_),
        "Same lengths" = all(diff(lengths(.series_)) == 0L),
        "Centered" = all(sapply(.series_, function(series) {
            if (is.null(dim(series)))
                mu <- mean(series)
            else
                mu <- colMeans(series)
            all(sapply(mu, function(m) { isTRUE(all.equal(0, m)) }))
        })),
        "Scaled" = all(sapply(.series_, function(series) {
            if (is.null(dim(series)))
                sigma <- sd(series)
            else
                sigma <- apply(series, 2L, sd)
            all(sapply(sigma, function(s) { isTRUE(all.equal(1, s)) }))
        })),
        check.names = FALSE
    )
})

# Evaluate tab, cvis
cvis_table <- quote({
    as.data.frame(rbind(cvi(result(), b = b, type = c(internal, external))))
})

# Evaluate tab, fcluster
fcluster_table <- quote({
    out@fcluster
})

# Evaluate tab, clusinfo
clusinfo_table <- quote({
    df <- t(out@clusinfo)
    rownames(df) <- c("Cluster size", "Average distance")
    colnames(df) <- paste("Cluster", 1L:out@k)
    df
})

# Evaluate tab, cl
cl_table <- quote({
    df <- cbind(out@cluster, out@cldist)
    colnames(df) <- c("Cluster ID", "Distance to centroid")
    rownames(df) <- names(.series_)
    df
})

# Evaluate tab, distmat
distmat_table <- quote({
    as.data.frame(dm)
})
