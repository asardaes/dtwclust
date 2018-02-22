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

# Evaluate tab, raw
raw_table <- quote({
    df <- this_result$results[[1L]]
    df <- df[setdiff(colnames(df), colnames(this_result$scores[[1L]]))]
    df
})

