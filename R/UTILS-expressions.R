prepare_expr <- quote({
    pairwise <- isTRUE(pairwise)
    dim_names <- list(names(x), names(y))
    if (pairwise && length(x) != length(y))
        stop("Pairwise distances require the same amount of series in 'x' and 'y'.")
    # Get appropriate matrix (big.matrix was needed before)
    D <- allocate_distmat(length(x), length(y), pairwise, symmetric, lower_triangular_only, diagonal) # UTILS-utils.R

    fill_type <- if (pairwise) {
        "PAIRWISE"
    } else if (lower_triangular_only) {
        if (diagonal) "LOWER_TRIANGULAR_DIAGONAL" else "LOWER_TRIANGULAR"
    } else if (symmetric) {
        "SYMMETRIC"
    } else {
        "PRIMARY"
    }

    mat_type <- "R_MATRIX"
})
