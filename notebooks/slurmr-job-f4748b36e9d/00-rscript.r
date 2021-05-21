.libPaths(c("/home/baker02/miniconda3/envs/tidyscreen/lib/R/library"))
message("[slurmR info] Loading variables and functions... ", appendLF = FALSE)
Slurm_env <- function (x = "SLURM_ARRAY_TASK_ID") 
{
    y <- Sys.getenv(x)
    if ((x == "SLURM_ARRAY_TASK_ID") && y == "") {
        return(1)
    }
    y
}
ARRAY_ID  <- as.integer(Slurm_env("SLURM_ARRAY_TASK_ID"))

# The -snames- function creates the write names for I/O of files as a 
# function of the ARRAY_ID
snames    <- function (type, array_id = NULL, tmp_path = NULL, job_name = NULL) 
{
    if (length(array_id) && length(array_id) > 1) 
        return(sapply(array_id, snames, type = type, tmp_path = tmp_path, 
            job_name = job_name))
    type <- switch(type, r = "00-rscript.r", sh = "01-bash.sh", 
        out = "02-output-%A-%a.out", rds = if (missing(array_id)) "03-answer-%03i.rds" else sprintf("03-answer-%03i.rds", 
            array_id), job = "job.rds", stop("Invalid type, the only valid types are `r`, `sh`, `out`, and `rds`.", 
            call. = FALSE))
    sprintf("%s/%s/%s", tmp_path, job_name, type)
}
TMP_PATH  <- "/mnt/scratcha/fmlab/baker02/tidyscreen/notebooks"
JOB_NAME  <- "slurmr-job-f4748b36e9d"

# The -tcq- function is a wrapper of tryCatch that on error tries to recover
# the message and saves the outcome so that slurmR can return OK.
tcq <- function (...) 
{
    ans <- tryCatch(..., error = function(e) e)
    if (inherits(ans, "error")) {
        ARRAY_ID. <- get("ARRAY_ID", envir = .GlobalEnv)
        msg <- paste("An error has ocurred while evualting the expression:\n", 
            paste(deparse(match.call()[[2]]), collapse = "\n"), 
            "\n in ", "ARRAY_ID # ", ARRAY_ID.)
        warning(msg, immediate. = TRUE, call. = FALSE)
        ans$message <- paste(ans$message, msg)
        saveRDS(ans, snames("rds", tmp_path = get("TMP_PATH", 
            envir = .GlobalEnv), job_name = get("JOB_NAME", envir = .GlobalEnv), 
            array_id = ARRAY_ID.))
        q("no")
    }
    invisible(ans)
}
message("done loading variables and functions.")
message("[slurmR info] Loading packages ... ")
tcq({
  library(dynfeature, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(dynguidelines, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(dynmethods, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(dynplot, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(dynwrap, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(dyno, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(purrr, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(matrixStats, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(MatrixGenerics, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(BiocGenerics, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(S4Vectors, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(IRanges, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(GenomeInfoDb, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(GenomicRanges, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(Biobase, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(SummarizedExperiment, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(SingleCellExperiment, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(scran, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(ggplot2, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(scater, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(dyngen, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(igraph, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(dynutils, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(tidyverse, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(tibble, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(tidyr, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(readr, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(dplyr, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(stringr, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(forcats, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(assertthat, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
  library(data.table, lib.loc = "/home/baker02/miniconda3/envs/tidyscreen/lib/R/library")
})
message("[slurmR info] done loading packages.")
tcq({
  INDICES <- readRDS("/mnt/scratcha/fmlab/baker02/tidyscreen/notebooks/slurmr-job-f4748b36e9d/INDICES.rds")
})
tcq({
  X <- readRDS(sprintf("/mnt/scratcha/fmlab/baker02/tidyscreen/notebooks/slurmr-job-f4748b36e9d/X_%04d.rds", ARRAY_ID))
})
tcq({
  FUN <- readRDS("/mnt/scratcha/fmlab/baker02/tidyscreen/notebooks/slurmr-job-f4748b36e9d/FUN.rds")
})
tcq({
  mc.cores <- readRDS("/mnt/scratcha/fmlab/baker02/tidyscreen/notebooks/slurmr-job-f4748b36e9d/mc.cores.rds")
})
tcq({
  seeds <- readRDS("/mnt/scratcha/fmlab/baker02/tidyscreen/notebooks/slurmr-job-f4748b36e9d/seeds.rds")
})
set.seed(seeds[ARRAY_ID], kind = NULL, normal.kind = NULL)
tcq({
  ans <- parallel::mclapply(
    X                = X,
    FUN              = FUN,
    mc.cores         = mc.cores
)
})
saveRDS(ans, sprintf("/mnt/scratcha/fmlab/baker02/tidyscreen/notebooks/slurmr-job-f4748b36e9d/03-answer-%03i.rds", ARRAY_ID), compress = TRUE)
