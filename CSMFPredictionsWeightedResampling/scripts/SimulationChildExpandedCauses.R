library(here)
### i is row of setting.df
i <- as.numeric(commandArgs(trailingOnly = TRUE))
### Make data directory for this run
sim.dir <- here("CSMFPredictionsWeightedResampling", "cluster_output_child_expanded_causes")
if(!dir.exists(sim.dir)) {
    dir.create(sim.dir, recursive = TRUE) 
}

### Quit if results file already exists
results.file <- file.path(sim.dir, paste0("run-", i, ".rds"))
if(file.exists(results.file)){
    quit('no')
}

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
setting.df <- expand.grid(run = 1:500,
                          country = c("India", "Tanzania"),
                          calib.size = c(50, 100, 200))



setting <- setting.df[i,]
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
### Going to keep Pneumonia, Diarrhea/Dysentery, Sepsis, Malaria, 
# other_infectious <- c("AIDS",
#                       "Malaria",
#                       "Other Infectious Diseases",
#                       "Encephalitis",
#                       "Hemorrhagic fever",
#                       "Measles",
#                       "Meningitis",
#                       "Sepsis")
other_infectious <- c("AIDS",
                      "Other Infectious Diseases",
                      "Encephalitis",
                      "Hemorrhagic fever",
                      "Measles",
                      "Meningitis")
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

### Read in the CSMFs that we will be using for ptest
### this should correspond to NCOD

### Get pcalib
top.cause.df$pcalib <- top.cause.df$ptest
calib_ind <- sample(seq_len(nrow(country.data)), size = setting$calib.size,
                    replace = F)


test.final <- country.data
calib.final <- country.data[calib_ind,]


### Get causes (these will be the ordering we will use)
causes <- as.character(top.cause.df$cause)
################################

################################
### Now the 6 methods w/o calibration
### Tariff with just training set (will also predict on calibration data)

tariff_file <- here("CSMFPredictionsWeightedResampling",
                    "child_data", paste0("tariff_", setting$country, ".rds"))
tariff_cod <- readRDS(tariff_file)
tariff.train.cod <- cause.df$broad[match(tariff_cod, cause.df$va34)]

### Insilico with just training set
insilico_file <- here("CSMFPredictionsWeightedResampling", "child_data", paste0("insilico_", setting$country, ".rds"))
insilico_cod <- readRDS(insilico_file)
insilico.train.cod <- cause.df$broad[match(insilico_cod, cause.df$va34)]


### Get the COD estimates for test & calibration set from tariff & insilicova
tariff.train.cod.test <- tariff.train.cod
tariff.train.cod.calib <- tariff.train.cod[calib_ind]

insilico.train.cod.test <- insilico.train.cod
insilico.train.cod.calib <- insilico.train.cod[calib_ind]

#Calibration truth
calib.truth <- changeTopCOD(calib.final$Cause)

### CalibratedVA

alpha <- 5
beta <- .5
tau <- .5
epsilon <- .001
tau.vec <- rep(tau, length(causes))
delta <- 1
gamma.init <- 1
ndraws <- 50E3
nchains <- 3

### calibva with Tariff
set.seed(123)
calibva.seeds <- sample(1e6, nchains, replace = F)
tariff.calibva <-  calibva.sampler(test.cod = tariff.train.cod.test,
                                   calib.cod = tariff.train.cod.calib,
                                   calib.truth = calib.truth, causes = causes,
                                   epsilon = epsilon, alpha=alpha, beta=beta,
                                   tau.vec=tau.vec, delta=delta,
                                   gamma.init=gamma.init, ndraws = ndraws,
                                   nchains = nchains,
                                   init.seed = calibva.seeds)

### CalibratedVA with InSilicoVA
insilico.calibva <- calibva.sampler(test.cod = insilico.train.cod.test,
                                    calib.cod = insilico.train.cod.calib,
                                    calib.truth = calib.truth, causes = causes,
                                    epsilon = epsilon, alpha=alpha, beta=beta,
                                    tau.vec=tau.vec, delta=delta,
                                    gamma.init=gamma.init, ndraws = ndraws,
                                    nchains = nchains,
                                    init.seed = calibva.seeds)

### calibva with Ensemble lite sampler
test.cod.mat <- matrix(c(tariff.train.cod.test,  insilico.train.cod.test), ncol = 2)
calib.cod.mat <- matrix(c(tariff.train.cod.calib,  insilico.train.cod.calib), ncol = 2)

ensemble.calibva <- calibva.ensemble.lite.sampler(test.cod.mat = test.cod.mat,
                                                  calib.cod.mat = calib.cod.mat,
                                                  calib.truth = calib.truth,
                                                  causes = causes,
                                                  epsilon = epsilon, alpha=alpha,
                                                  beta=beta,
                                                  tau.vec=tau.vec, delta=delta,
                                                  gamma.init=gamma.init,
                                                  ndraws = ndraws,
                                                  nchains = nchains,
                                                  init.seed = calibva.seeds)

### Get posterior draws for CSMF parameters
burnin <- 10E3
thin <- 10
tariff.calibva <- window(tariff.calibva, start = burnin, thin = thin)
insilico.calibva <- window(insilico.calibva, start = burnin, thin = thin)
ensemble.calibva <- window(ensemble.calibva, start = burnin, thin = thin)


tariff.calibva.csmf.df <- calibvaCSMFPosteriorSamples(tariff.calibva, causes = causes)
insilico.calibva.csmf.df <- calibvaCSMFPosteriorSamples(insilico.calibva, causes = causes)
ensemble.calibva.csmf.df <- calibvaCSMFPosteriorSamples(ensemble.calibva, causes = causes)

### Function to get CSMF estimate from CalibratedVA
calibratedVACSMF <- function(calib_df, causes) {
    csmf_df <- calibvaCSMFPosteriorSummary(calib_df, percentile.L = .5)
    m <- match(causes, csmf_df$cause)
    #return(csmf_df$ci.L[m])
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
ensemble.calibva.csmf <- calibratedVACSMF(ensemble.calibva.csmf.df, causes)


tariff.train.csmf <- openVACSMF(tariff.train.cod.test, causes)
insilico.train.csmf <- openVACSMF(insilico.train.cod.test, causes)
calib.csmf <- openVACSMF(calib.truth, causes)

### CSMF  data frame
methods <- c("tariff_calibva",
             "insilico_calibva",
             "ensemble_calibva",
             "tariff_train",
             "insilico_train",
             "calib")
csmf.df.expanded <- data.frame(csmf.est = c(tariff.calibva.csmf,
                                            insilico.calibva.csmf,
                                            ensemble.calibva.csmf,
                                            tariff.train.csmf,
                                            insilico.train.csmf,
                                            calib.csmf),
                               cause = rep(causes, length(methods)),
                               cause.text = rep(top.cause.df$cause, length(methods)),
                               true.csmf = rep(top.cause.df$ptest, length(methods)),
                               method = rep(methods, each = nrow(top.cause.df)))
csmf.acc.expanded <-
    csmf.df.expanded %>%
    group_by(method) %>%
    summarize(csmfa = getCSMF_accuracy(csmf.est, true.csmf))

### Now re-classify sepsis & malaria to Other Infectious
csmf.df <- csmf.df.expanded
csmf.df$cause <- as.character(csmf.df$cause)
csmf.df$cause <- ifelse(csmf.df$cause %in% c("Sepsis", "Malaria"),
                                 "Other Infectious", csmf.df$cause)
csmf.df <-
    csmf.df %>%
    group_by(cause, method) %>%
    summarize(csmf.est = sum(csmf.est),
              true.csmf = sum(true.csmf)) %>%
    ungroup()


### Now CSMF Accuracy
csmf.acc.df <-
    csmf.df %>%
    group_by(method) %>%
    summarize(csmfa = getCSMF_accuracy(csmf.est, true.csmf))

results <- list(topcause.df = top.cause.df,
                csmf = csmf.df,
                csmf.expanded = csmf.df.expanded,
                csmf.acc = csmf.acc.df,
                csmf.acc.expanded = csmf.acc.expanded,
                setting = setting) 
results.file <- file.path(sim.dir, paste0("run-", i, ".rds"))
saveRDS(results, results.file)
quit('no')

