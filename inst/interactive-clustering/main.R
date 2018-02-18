# main clustering routine
main <- quote({
    shinyjs::disable("cluster__cluster")
    this_result <- tryCatch({
        type <- input$cluster__clus_type
        k <- as.integer(input$cluster__k)
        distance <- input$cluster__dist
        seed <- as.integer(input$cluster__seed)
        if (seed == 0L) seed <- NULL
        trace <- input$cluster__trace
        if (input$cluster__cent_custom)
            centroid <- match.fun(input$cluster__cent_func)
        else
            centroid <- input$cluster__cent
        # controls
        control <- switch(
            input$cluster__clus_type,
            "p" = {
                partitional_control(
                    iter.max = as.integer(input$cluster__part_iter),
                    pam.precompute = input$cluster__part_pam,
                    symmetric = input$cluster__part_symmetric,
                    distmat = distmat()
                )
            },
            "h" = {
                if (input$cluster__hier_method_custom)
                    method <- match.fun(input$cluster__hier_method_func)
                else
                    method <- input$cluster__hier_method

                hierarchical_control(
                    method = method,
                    symmetric = input$cluster__hier_symmetric,
                    distmat = distmat()
                )
            },
            "f" = {
                fuzzy_control(
                    iter.max = as.integer(input$cluster__fuzz_iter),
                    fuzziness = input$cluster__fuzz_m,
                    delta = input$cluster__fuzz_delta,
                    symmetric = input$cluster__fuzz_symmetric,
                    distmat = distmat()
                )
            },
            "t" = {
                tadpole_control(
                    dc = input$cluster__tadp_dc,
                    window.size = as.integer(input$cluster__tadp_window),
                    lb = input$cluster__tadp_lb
                )
            }
        )
        # dist args
        dist_args <- input$cluster__dist_args
        dist_args <- if (nzchar(dist_args)) parse_input(dist_args) else list()
        # cent args
        cent_args <- input$cluster__cent_args
        cent_args <- if (nzchar(cent_args)) parse_input(cent_args) else list()
        # dots
        dots <- input$cluster__dots
        dots <- if (nzchar(dots)) parse_input(dots) else list()
        # return
        args <- enlist(
            series = .series_,
            type = type,
            k = k,
            distance = distance,
            centroid = centroid,
            control = control,
            seed = seed,
            trace = trace,
            error.check = FALSE,
            args = tsclust_args(dist = dist_args, cent = cent_args),
            dots = dots
        )
        if (is.character(centroid) && centroid == "default")
            args$centroid <- NULL
        if (type == "t")
            args$distance <- NULL
        do.call(tsclust, args, TRUE)
    },
    error = function(e) {
        e
    })
    shinyjs::enable("cluster__cluster")
    if (inherits(this_result, "error")) {
        shinyjs::alert(this_result$message)
    }
    else {
        result(this_result)
        distmat(this_result@distmat)
    }
})
