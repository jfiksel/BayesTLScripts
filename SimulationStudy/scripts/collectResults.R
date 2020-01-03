library(dplyr)
library(here)
setting.df <- expand.grid(run = 1:100, Mtype=1:3, calib.size = c(50, 100, 200, 400),
                          fake=c("tariff","insilico"), csmfa.calib=c("low","medium","high"))
results.list <- lapply(1:nrow(setting.df), function(i){
    #print(i)
    setting <- setting.df[i,]
    sim.dir <- here("SimulationStudy", "cluster_output", paste0("M",setting$Mtype),
                    setting$fake,
                    paste0("size_", setting$calib.size),
                    paste0("calib_",setting$csmfa.calib),
                    paste0("run_", setting$run))
    results.file <- file.path(sim.dir, "results.rds")
    results <- readRDS(results.file)
    csmf.acc.df <-
        results$csmf.acc  %>%
        rename(accuracy = csmf.acc)
    avg.ccc.df <-
        results$avg.ccc %>%
        rename(accuracy = ccc.mean)
    combined.results <- rbind(csmf.acc.df, avg.ccc.df)
    
    causes=c("Pneumonia", "Diarrhea/Dysentery", "Sepsis", "Other")
    
    for(c in causes) {
        cause.acc.df <- results$csmf %>% 
            filter(cause.text==c) %>% 
            transmute(accuracy = (csmf.est - true.csmf), method) 
            combined.results=rbind(combined.results,cause.acc.df)
    }
    
    combined.results$measure <- c(rep(c('csmf', 'ccc'),c(nrow(csmf.acc.df),nrow(avg.ccc.df))),
     unlist(sapply(causes,rep,nrow(cause.acc.df))))
    combined.results <- cbind(combined.results, setting)    
    combined.results$csmfa.test.train=results$csmf.test.train
    combined.results$csmfa.test.calib=results$csmf.test.calib
    #combined.results[,colnames(setting.df)]=setting
    return(combined.results)
})
results.df <- do.call(rbind, results.list)
data.dir <-  here("SimulationStudy", "data")
if(!dir.exists(data.dir)){
    dir.create(data.dir, recursive = TRUE)
}
saveRDS(results.df, file.path(data.dir, "simulation_results.rds"))
