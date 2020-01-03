library(tidyverse)
library(here)

### Read main child & adult results
#for(phmrc.type in c("child", 'adult')) {
phmrc.type <- 'child'
files <- list.files(here("CSMFPredictionsWeightedResampling",
                         paste0("cluster_output_", phmrc.type)),
                    pattern = "*.rds", full.names = TRUE)
results.list <- lapply(seq_along(files), function(i){
    x <- readRDS(files[i])
    return(cbind(x$csmf.acc, x$setting))
})
results.df <- do.call(rbind, results.list)
data.dir <- here("CSMFPredictionsWeightedResampling", paste0(phmrc.type, "_data"))
saveRDS(results.df, file.path(data.dir, "csmfa_results.rds"))
#}

### Get CSMF estimates


### Read main child & adult results
files <- list.files(here("CSMFPredictionsWeightedResampling", "cluster_output_child"),
                    full.names = TRUE)
results.list <- lapply(seq_along(files), function(i){
    x <- readRDS(files[i])
    return(cbind(x$csmf, x$setting))
})
results.df <- do.call(rbind, results.list)
data.dir <- here("CSMFPredictionsWeightedResampling", "child_data")
saveRDS(results.df, file.path(data.dir, "csmf_results.rds"))


### Read child data with expanded causes

files <- list.files(here("CSMFPredictionsWeightedResampling",
                         "cluster_output_child_expanded_causes"), full.names = TRUE)
results.list <- lapply(seq_along(files), function(i){
    x <- readRDS(files[i])
    return(cbind(x$csmf.acc, x$setting))
})
results.df <- do.call(rbind, results.list)
data.dir <- here("CSMFPredictionsWeightedResampling", "child_data")
saveRDS(results.df, file.path(data.dir, "csmfa_expanded_causes_results.rds"))

