.onUnload <- function(libpath) {
  .stop_cluster()
  library.dynam.unload("Rcausality", libpath)
}

.onLoad <- function(libname, pkgname) {
  # the PSOCK cluster cached on Windows is closed when the session ends
  reg.finalizer(.par_env, function(e) {
    cl <- e$cluster
    if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
  }, onexit = TRUE)
  invisible(NULL)
}
