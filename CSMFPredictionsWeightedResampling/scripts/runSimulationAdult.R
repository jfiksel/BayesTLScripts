library(here)
id <- as.numeric(commandArgs(trailingOnly = TRUE))
sim.dir <- here("CSMFPredictionsWeightedResampling", "cluster_output_adult")
if(!dir.exists(sim.dir)) {
    dir.create(sim.dir, recursive = TRUE) 
}
results.file <- file.path(sim.dir, paste0("run-", id, ".rds"))
if(file.exists(results.file)){
    quit('no')
}
library(tidyverse)
library(ggmcmc)
library(CalibratedVA)
library(coda)
#library(openVA)
library(gtools)

countries <- c("India", "Mexico", "Philippines", "Tanzania")
setting.df <- expand.grid(run = 1:500,
                          country = c("India", "Tanzania"),
                          calib.size = c(50, 100, 200, 400),
                          calib.csfma = c("low", "representative"))
setting <- setting.df[id,]
seed.index <- setting$run
c <- as.character(setting$country)
calib.size <- setting$calib.size

phmrc <- readRDS(here("CSMFPredictionsWeightedResampling","adult_data", "phmrc_adult.rds"))
country.df <- data.frame(site = c("AP", "Bohol", "Dar", "Mexico", "Pemba", "UP"),
                         country = c("India", 
                                     "Philippines",
                                     "Tanzania",
                                     "Mexico",
                                     "Tanzania",
                                     "India"),
                          stringsAsFactors=FALSE)
country <- country.df$country[match(phmrc$site, country.df$site)]
country.data <- phmrc[country == setting$country,]

causes <- c("external", "circulatory","non_communicable", "infectious", "maternal")
C <- length(causes)

### Cause map
external <- c("Road Traffic",
              "Falls",
              "Homicide",
              "Suicide",
              "Fires",
              "Drowning",
              "Other Injuries",
              "Poisonings",
              "Bite of Venomous Animal")
circulatory <- c("Stroke",
                 "Other Cardiovascular Diseases",
                 "Acute Myocardial Infarction")
non_communicable <- c("Other Non-communicable Diseases",
                      "Colorectal Cancer",
                      "Breast Cancer",
                      "Leukemia/Lymphomas",
                      "Prostate Cancer",
                      "Esophageal Cancer",
                      "Stomach Cancer",
                      "Lung Cancer",
                      "Cervical Cancer",
                      "Renal Failure",
                      "Epilepsy",
                      "Cirrhosis",
                      "COPD",
                      "Diabetes",
                      "Asthma")
infectious <- c("Malaria",
                "Pneumonia",
                "Diarrhea/Dysentery",
                "AIDS",
                "TB",
                "Other Infectious Diseases")
maternal <- c("Maternal")

cause.df <- data.frame(causes = c(external, circulatory, non_communicable,
                                  infectious, maternal),
                       broad = rep(c("external", "circulatory",
                                     "non_communicable", "infectious", "maternal"),
                                   c(length(external), length(circulatory),
                                     length(non_communicable), length(infectious),
                                     length(maternal))),
                       stringsAsFactors=FALSE)


top.cause.df <- data.frame(cause = unique(cause.df$broad))


### Get distribution of these causes in the train set
### first, function to change all guesses not equal to top 3 causes to other
changeTopCOD <- function(topcod) {
    topcod <- as.character(topcod)
    m <- match(topcod, cause.df$causes)
    return(cause.df$broad[m])
}


top.cause.df$ptest <- sapply(top.cause.df$cause, function(c) mean(changeTopCOD(country.data$gs_text34) == c))
top.cause.df$ntest <- sapply(top.cause.df$cause, function(c) sum(changeTopCOD(country.data$gs_text34) == c))

#########################################
### Get calibration and test indices
getCSMF_accuracy <- function(csmf, truth){
    acc <- 1 - sum(abs(truth - csmf)) / (2  * (1 - min(truth)))
    return(acc)
}
set.seed(123)
seeds <- sample(-1e6:1e6, size = 500, replace = F)
init.seed = seeds[seed.index]
set.seed(init.seed)

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
        pcalib.mat=sapply(1:100000,function(i) rdirichlet_mod(rep(1,nrow(top.cause.df))))
        csmfa.train=apply(pcalib.mat,2,getCSMF_accuracy, top.cause.df$ptest)
        ind.calib=which(csmfa.train > 0 & 
                            csmfa.train < .4)
        if(length(ind.calib) > 0) flag=1
    }
    
    pcalib=pcalib.mat[,sample(ind.calib,1)]
    top.cause.df$pcalib <- pcalib
    
    ### Sample according to probabilities
    country_causes <- changeTopCOD(country.data$gs_text34)
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

### read in models
models_dir <- here("CSMFPredictionsWeightedResampling","adult_data")
insilico_cod <- as.character(readRDS(file.path(models_dir, paste0("insilico_model_", c, "_cod.rds"))))
tariff_cod <- as.character(readRDS(file.path(models_dir, paste0("tariff_model_", c, "_cod.rds"))))
nbc_cod <- as.character(readRDS(file.path(models_dir, paste0("nbc_model_", c, "_cod.rds"))))

### run calibration

test.final <- country.data
calib.final <- country.data[calib_ind,]


### Get causes (these will be the ordering we will use)
causes <- as.character(top.cause.df$cause)
################################

#Calibration truth
calib.truth <- changeTopCOD(calib.final$gs_text34)

alpha <- 5
beta <- .5
tau <- .5
epsilon <- .001
tau.vec <- rep(tau, length(causes))
delta <- 1
gamma.init <- 1
ndraws <- 50E3
nchains <- 3

### calibva with InSilicoVA
set.seed(123)
calibva.seeds <- sample(1e6, nchains, replace = F)
insilico.calibva <-  calibva.sampler(test.cod = insilico_cod,
                                     calib.cod = insilico_cod[calib_ind],
                                     calib.truth = calib.truth, causes = causes,
                                     epsilon = epsilon, alpha=alpha, beta=beta,
                                     tau.vec=tau.vec, delta=delta,
                                     gamma.init=gamma.init, ndraws = ndraws,
                                     nchains = nchains,
                                     init.seed = calibva.seeds)

tariff.calibva <-  calibva.sampler(test.cod = tariff_cod,
                                   calib.cod = tariff_cod[calib_ind],
                                   calib.truth = calib.truth, causes = causes,
                                   epsilon = epsilon, alpha=alpha, beta=beta,
                                   tau.vec=tau.vec, delta=delta,
                                   gamma.init=gamma.init, ndraws = ndraws,
                                   nchains = nchains,
                                   init.seed = calibva.seeds)

nbc.calibva <-  calibva.sampler(test.cod = nbc_cod,
                                calib.cod = nbc_cod[calib_ind],
                                calib.truth = calib.truth, causes = causes,
                                epsilon = epsilon, alpha=alpha, beta=beta,
                                tau.vec=tau.vec, delta=delta,
                                gamma.init=gamma.init, ndraws = ndraws,
                                nchains = nchains,
                                init.seed = calibva.seeds)

### ensemble
test.cod.mat <- matrix(c(insilico_cod,
                         tariff_cod,
                         nbc_cod),
                       ncol = 3)
calib.cod.mat <- matrix(c(insilico_cod[calib_ind],
                          tariff_cod[calib_ind],
                          nbc_cod[calib_ind]),
                        ncol = 3)
ensemble.calibva <- calibva.ensemble.lite.sampler(test.cod.mat = test.cod.mat,
                                                  calib.cod.mat = calib.cod.mat,
                                                  calib.truth = calib.truth, causes = causes,
                                                  epsilon = epsilon, alpha=alpha, beta=beta,
                                                  tau.vec=tau.vec, delta=delta,
                                                  gamma.init=gamma.init, ndraws = ndraws,
                                                  nchains = nchains,
                                                  init.seed = calibva.seeds)

burnin <- 10E3
thin <- 10
tariff.calibva <- window(tariff.calibva, start = burnin, thin = thin)
insilico.calibva <- window(insilico.calibva, start = burnin, thin = thin)
nbc.calibva <- window(nbc.calibva, start = burnin, thin = thin)
ensemble.calibva <- window(ensemble.calibva, start = burnin, thin = thin)


tariff.calibva.csmf.df <- calibvaCSMFPosteriorSamples(tariff.calibva, causes = causes)
insilico.calibva.csmf.df <- calibvaCSMFPosteriorSamples(insilico.calibva, causes = causes)
nbc.calibva.csmf.df <- calibvaCSMFPosteriorSamples(nbc.calibva, causes = causes)
ensemble.calibva.csmf.df <- calibvaCSMFPosteriorSamples(ensemble.calibva, causes = causes)

### Function to get CSMF estimate from CalibratedVA
calibratedVACSMF <- function(calib_df, causes) {
    csmf_df <- calibvaCSMFPosteriorSummary(calib_df)
    m <- match(causes, csmf_df$cause)
    return(csmf_df$mean[m])
}

### Function to get CSMF estimate from openVA
openVACSMF <- function(topcod, causes) {
    csmf <- sapply(causes, function(c) mean(topcod == c))
    return(unname(csmf))
}

### Get CSMF estimates 
tariff.calibva.csmf <- calibratedVACSMF(tariff.calibva.csmf.df, causes)
insilico.calibva.csmf <- calibratedVACSMF(insilico.calibva.csmf.df, causes)
nbc.calibva.csmf <- calibratedVACSMF(nbc.calibva.csmf.df, causes)
ensemble.calibva.csmf <- calibratedVACSMF(ensemble.calibva.csmf.df, causes)


tariff.train.csmf <- openVACSMF(tariff_cod, causes)
insilico.train.csmf <- openVACSMF(insilico_cod, causes)
nbc.train.csmf <- openVACSMF(nbc_cod, causes)
calib.csmf <- openVACSMF(calib.truth, causes)

### CSMF  data frame
methods <- c("tariff_calibva",
             "insilico_calibva",
             "nbc_calibva",
             "ensemble_calibva",
             "tariff_train",
             "insilico_train",
             "nbc_train",
             "calib")
csmf.df <- data.frame(csmf.est = c(tariff.calibva.csmf,
                                   insilico.calibva.csmf,
                                   nbc.calibva.csmf,
                                   ensemble.calibva.csmf,
                                   tariff.train.csmf,
                                   insilico.train.csmf,
                                   nbc.train.csmf,
                                   calib.csmf),
                      cause = rep(causes, length(methods)),
                      true.csmf = rep(top.cause.df$ptest, length(methods)),
                      method = rep(methods, each = nrow(top.cause.df)))

### Now CSMF Accuracy
csmf.acc.df <-
    csmf.df %>%
    group_by(method) %>%
    summarize(csmfa = getCSMF_accuracy(csmf.est, true.csmf))

results <- list(topcause.df = top.cause.df,
                csmf = csmf.df,
                csmf.acc = csmf.acc.df,
                setting = setting) 
saveRDS(results, results.file)
quit('no')


