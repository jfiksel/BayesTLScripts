library(here)
library(tidyverse)
child_results <- readRDS(here("CSMFPredictionsWeightedResampling", "child_data", "csmfa_results.rds"))

### Run for InSilico in India to see why low does worse than representative
insilico_india_setting_low <-
    child_results %>%
    filter(method == "insilico_calibva", country == "India", calib.size == 200,
           cod == 5, alpha == 5, calib.csfma == "low") %>%
    mutate(median_csmfa = median(csmfa),
           diff = abs(csmfa - median_csmfa)) %>%
    slice(which.min(diff))

insilico_india_setting_rep <-
    child_results %>%
    filter(method == "insilico_calibva", country == "India", calib.size == 200,
           cod == 5, alpha == 5, calib.csfma == "representative") %>%
    mutate(median_csmfa = median(csmfa),
           diff = abs(csmfa - median_csmfa)) %>%
    slice(which.min(diff))

insilico_india_setting <- bind_rows(insilico_india_setting_low, insilico_india_setting_rep)


### Runs for InSilico in Tanzania, comparing cod = 5 to cod = 7 
insilico_tanzania_setting <-
    child_results %>%
    filter(method == "insilico_calibva", country == "Tanzania", calib.size %in% c(50, 200),
           cod == 5, alpha == 5, calib.csfma == "low") %>%
    group_by(calib.size) %>%
    mutate(median_csmfa = median(csmfa),
           diff = abs(csmfa - median_csmfa)) %>%
    slice(which.min(diff))
insilico_tanzania_setting <- insilico_tanzania_setting[c(1:2, 1:2),]
insilico_tanzania_setting$cod[3:4] <- c(7, 7)

setting.df <- bind_rows(insilico_india_setting, insilico_tanzania_setting)

calibration_list <- lapply(1:nrow(setting.df), function(i) {
    setting <- setting.df[i,]
    calib.dir <- here("CSMFPredictionsWeightedResampling", "cluster_output_child", "calib_indices")
    calib.file <- file.path(calib.dir,
                            paste0(setting$country,
                                   "run", setting$run,
                                   "calibsize", setting$calib.size,
                                   "csfma", setting$calib.csfma, ".rds"))
    calib_ind <- readRDS(calib.file)
    return(list(setting = setting, calib_ind = calib_ind))
})
saveRDS(calibration_list, here("CSMFPredictionsWeightedResampling", "child_data",
                               "calib_indices_for_eda.rds"))
