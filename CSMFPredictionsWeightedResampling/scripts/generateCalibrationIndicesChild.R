library(here)
### i is row of setting.df
i <- as.numeric(commandArgs(trailingOnly = TRUE))
setting.df <- expand.grid(run = 1:500,
                          country = c("India", "Tanzania"),
                          calib.size = c(50, 100, 200),
                          calib.csfma = c("low", "representative"))

setting <- setting.df[i,]
### Quit if results file already exists
sim.dir <- here("CSMFPredictionsWeightedResampling", "cluster_output_child", "calib_indices")
if(!dir.exists(sim.dir)) {
    dir.create(sim.dir, recursive = TRUE) 
}

results.file <- file.path(sim.dir,
                          paste0(setting$country,
                                 "run", setting$run,
                                 "calibsize", setting$calib.size,
                                 "csfma", setting$calib.csfma, ".rds"))
if(file.exists(results.file)){
    quit('no')
}
### Make data directory for this run

library(openVA)
library(CalibratedVA)
library(tidyverse)
library(coda)
library(ggmcmc)
library(gtools)
set.seed(123)
### Create list of seeds
all.seeds <- sample(1e6, size = 500, replace = F)
### Read in child data
child.raw <- read.csv(getPHMRC_url("child"))
child.raw$gs_text34 <- as.character(child.raw$gs_text34)
### Get data frame matching va34 code to actual COD
cause.df <- unique(child.raw[,c("gs_text34", "va34")])
cause.df$va34 <- as.character(as.numeric(cause.df$va34))
### Clean data into output usable by Tariff & InsilicoVA
child.clean <- ConvertData.phmrc(child.raw, phmrc.type = "child")$output
child.clean$Cause <- cause.df$gs_text34[match(child.clean$Cause, cause.df$va34)]
### Assign countries
countries <- rep(NA, nrow(child.raw))
countries[child.raw$site %in% c("AP", "UP")] <- "India"
countries[child.raw$site %in% c("Mexico")] <- "Mexico"
countries[child.raw$site %in% c("Dar", "Pemba")] <- "Tanzania"
countries[child.raw$site %in% c("Bohol")] <- "Philippines"
### Create setting data frame with number of runs (1-500) and each split type on it
### ncod will be the number of total causes: top (ncod - 1) causes + other

ncod <- setting$ncod
set.seed(all.seeds[setting$run])
### split data into country we are using and other
country.data <- child.clean[countries == setting$country,]
train.final <- child.clean[countries != setting$country,]
### Create cause map
external <- c("Road Traffic",
              "Falls",
              "Fires",
              "Drowning",
              "Poisonings",
              "Bite of Venomous Animal",
              "Violent Death")
### Going to keep Pneumonia & Diarrhea/Dysentery
other_infectious <- c("AIDS",
                      "Malaria",
                      "Other Infectious Diseases",
                      "Encephalitis",
                      "Hemorrhagic fever",
                      "Measles",
                      "Meningitis",
                      "Sepsis")
other <- c("Other Cardiovascular Diseases",
           "Other Cancers",
           "Other Defined Causes of Child Deaths",
           "Other Digestive Diseases")

cause.df$broad <- ifelse(cause.df$gs_text34 %in% external, "External",
                         ifelse(cause.df$gs_text34 %in% other, "Other",
                                ifelse(cause.df$gs_text34 %in% other_infectious, "Other Infectious",
                                       cause.df$gs_text34)))


### get causes
top.cause.df <- data.frame(cause = unique(cause.df$broad))


### Get distribution of these causes in the train set
### first, function to change all guesses not equal to top 3 causes to other
changeTopCOD <- function(topcod) {
    topcod <- as.character(topcod)
    m <- match(topcod, cause.df$gs_text34)
    return(cause.df$broad[m])
}

top.cause.df$ptrain <- sapply(top.cause.df$cause, function(c) mean(changeTopCOD(train.final$Cause) == c))
top.cause.df$ptest <- sapply(top.cause.df$cause, function(c) mean(changeTopCOD(country.data$Cause) == c))
top.cause.df$ntest <- sapply(top.cause.df$cause, function(c) sum(changeTopCOD(country.data$Cause) == c))


### Get pcalib
if(setting$calib.csfma == "low") {
    rdirichlet_mod=function(v) {
        p=rdirichlet(1,v)
        while(min(p) < 0.05 |  max(p) > .5) {
            p=rdirichlet(1,v)
        } 
        return(p)
    }
    
    flag=0
    while(flag==0){
        print("Generating CSMFs with desired CSMFA")
        pcalib.mat=sapply(1:100000,function(i) rdirichlet_mod(rep(1, nrow(top.cause.df))))
        csmfa.train=apply(pcalib.mat,2,getCSMF_accuracy, top.cause.df$ptest)
        ind.calib=which(csmfa.train > 0 & 
                            csmfa.train < .4)
        if(length(ind.calib) > 0) flag=1
    }
    
    pcalib=pcalib.mat[,sample(ind.calib,1)]
    top.cause.df$pcalib <- pcalib
    
    ### Sample according to probabilities
    country_causes <- changeTopCOD(country.data$Cause)
    ind_probs <- top.cause.df$pcalib[match(country_causes, top.cause.df$cause)]
    ind_n <- top.cause.df$ntest[match(country_causes, top.cause.df$cause)]
    ind_weights <- ind_probs / ind_n
    
    calib_ind <- sample(seq_len(nrow(country.data)), size = setting$calib.size,
                        replace = F, prob = ind_weights)
} else {
    top.cause.df$pcalib <- top.cause.df$ptest
    calib_ind <- sample(seq_len(nrow(country.data)), size = setting$calib.size,
                        replace = F)
}
saveRDS(calib_ind, results.file)
quit('no')
