unregister_parallel <- function() {
    doParallel::stopImplicitCluster()
}

register_parallel <- function(cores=NULL) {
    if (is.null(cores)) {
        cores <- parallel::detectCores()
    }
    cores <- max(1, cores)
    if (cores != 1) {
        doParallel::registerDoParallel(cores = cores)
    }
}
