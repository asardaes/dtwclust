# main clustering routine
main <- quote({
    shinyjs::disable("cluster__cluster")
    disable_buttons()
    best_window(NA_integer_)
    constraints(data.frame())
    window_flags(data.frame())
    this_result <- tryCatch({
        type <- input$cluster__clus_type
        k <- as.integer(input$cluster__k)
        distance <- "dtw_basic"
        seed <- as.integer(input$cluster__seed)
        if (seed == 0L) seed <- NULL
        trace <- input$cluster__trace
        if (input$cluster__cent_custom)
            centroid <- input$cluster__cent_func
        else
            centroid <- input$cluster__cent
        # windows
        window_sizes <- input$cluster__windows
        window_sizes <- seq(from = window_sizes[1L],
                            to = window_sizes[2L],
                            by = as.integer(input$cluster__windows_step))
        window_sizes <- as.integer(unique(round(min(lengths(.series_)) * window_sizes / 100)))
        if (window_sizes[1L] == 0L && !all(diff(lengths(.series_)) == 0L))
            stop("A window of size 0 should not be used with series of different length.")
        # controls
        control <- switch(
            type,
            "p" = {
                list(partitional = partitional_control(
                    iter.max = as.integer(input$cluster__part_iter),
                    pam.precompute = input$cluster__part_pam,
                    nrep = as.integer(input$cluster__part_nrep)
                ))
            },
            "h" = {
                if (input$cluster__hier_method_custom)
                    method <- input$cluster__hier_method_func
                else
                    method <- input$cluster__hier_method

                list(hierarchical = hierarchical_control(
                    method = method
                ))
            },
            "t" = {
                list(tadpole = tadpole_control(
                    dc = input$cluster__tadp_dc,
                    window.size = as.integer(window_sizes),
                    lb = input$cluster__tadp_lb
                ))
            }
        )
        # dist
        dist_args <- input$cluster__dist_args
        dist_args <- if (nzchar(dist_args)) parse_input(dist_args) else list()
        dist_args$window.size <- as.integer(window_sizes)
        # cent
        cent_args <- input$cluster__cent_args
        cent_args <- if (nzchar(cent_args)) parse_input(cent_args) else list()
        if (tolower(centroid) == "dba") cent_args$window.size <- as.integer(window_sizes)
        cent_cfg <- as.data.frame(c(list(centroid = centroid), cent_args), stringsAsFactors = FALSE)
        cent_cfg <- list(cent_cfg)
        names(cent_cfg) <- match.arg(type, c("partitional", "hierarchical", "tadpole"))
        # dots
        dots <- input$cluster__dots
        dots <- if (nzchar(dots)) parse_input(dots) else list()
        # configs
        cfgs <- compare_clusterings_configs(
            types = type,
            k = k,
            no.expand = "window.size",
            controls = control,
            preprocs = pdc_configs(
                "preproc",
                none = list(),
                share.config = type
            ),
            distances = suppressWarnings(pdc_configs(
                "distance",
                dtw_basic = dist_args,
                share.config = type
            )),
            centroids = cent_cfg
        )
        # return
        args <- enlist(
            series = .series_,
            type = type,
            configs = cfgs,
            seed = seed,
            trace = trace,
            score.clus = score_fun,
            dots = dots
        )
        do_call("compare_clusterings", args)
    },
    error = function(e) {
        e
    })
    shinyjs::enable("cluster__cluster")
    if (inherits(this_result, "error")) {
        shinyjs::alert(this_result$message)
    }
    else {
        # post-processing --------------------------------------------------------------------------
        result(this_result)
        output$evaluate__raw <- renderTable(raw_table, quoted = TRUE)

        agg_ids <- this_result$scores[[1L]]
        if (
            # partitional
            (input$cluster__clus_type == "p" && input$cluster__part_nrep > 1L) ||
            # hierarchical
            (input$cluster__clus_type == "h" && input$cluster__hier_method == "all"))
        {
            cfgs <- this_result$results[[1L]]$config_id
            agg_cfgs <- sapply(strsplit(cfgs, "_"), "[", 1L)
            split_ids <- split(agg_ids, agg_cfgs)
            method <- switch(input$cluster__clus_type,
                             "p" = input$cluster__part_agg,
                             "h" = input$cluster__hier_agg)
            agg_ids <- dplyr::bind_rows(lapply(split_ids, function(cfgs_ids) {
                partitions <- apply(cfgs_ids, 1L, function(ids) {
                    clue::as.cl_hard_partition(as.integer(ids))
                })
                ensemble <- clue::cl_ensemble(list = partitions)
                as.data.frame(rbind(unclass(clue::cl_medoid(ensemble, method)$.Data)))
            }))
            agg_ids <- data.frame(config_id = unique(agg_cfgs),
                                  window_size = window_sizes,
                                  agg_ids)
        }
        else {
            agg_ids <- data.frame(this_result$results[[1L]]["config_id"],
                                  window_size = window_sizes,
                                  agg_ids)
        }
        cluster_ids <<- agg_ids

        pair_tracker <<- PairTracker$new(length(.series_)) # S4-PairTracker.R
        enable_buttons()
        pair_ids(pair_tracker$get_unseen_pair())
    }
})
