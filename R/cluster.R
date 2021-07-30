#' @export
setup_cluster <- function (cluster, default_resources) {
    if (cluster == 'local') {
        default_resources$mem <- default_resources$njobs <- default_resources$walltime <- NULL
        default_resources$ncpus <- parallel::detectCores() - 1
        default_resources$cl <- parallel::makeCluster(default_resourecs$ncpus)

    } else if (cluster == "slurm") {
        plan(batchtools_slurm, template = template, resources = default_resources)
    }
}

#' @export
run_jobs <- function(inputs, FUN, njobs, ...) {
    # setting up all required parameters
    do_pb <- dopb()
    res <- listenv()
    
    if (do_pb) {
        pb <- startpb(0, njobs)
        on.exit(closepb(pb), add = TRUE)
    }
    
    # submit futures to cluster or cpu
    for (input in inputs) {
        res[[input]] %<-% { 
            FUN(input, ...)
        }
    }
    
    # monitor progression in computation and return list
    monitor_progress(res, njobs, do_pb, pb)
    return (as.list(res))
}

monitor_progress <- function(res, njobs, do_pb, pb) {
    resolved <- c()
    
    while (length(resolved) < njobs) {
        new_resolved <- c()

        for (i in seq_len(njobs)) {
            f <- futureOf(res[[i]])

            if (resolved(f) & !i %in% resolved) {
                resolved <- c(resolved, i)
                new_resolved <- c(new_resolved, i)
            }
        } 

        if (do_pb & length(new_resolved) > 0) {
            setpb(pb, length(resolved))
        }
    }
}