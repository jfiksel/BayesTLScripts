library(here)
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
setting.df <- expand.grid(run = 1,
                          country = c("Tanzania"),
                          calib.size = c(50, 100, 200),
                          cod = c("top3"),
                          calib.csfma = c("representative"))

posteriors <- lapply(1:nrow(setting.df), function(i) {
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
    if(setting$cod == "top3") {
        top_causes <- names(sort(table(country.data$Cause), decreasing = TRUE))[1:3]
        cause.df$broad <- ifelse(cause.df$gs_text34 %in% top_causes,
                                 as.character(cause.df$gs_text34), "Other")
    } else {
        cause.df$broad <- ifelse(cause.df$gs_text34 %in% external, "External",
                                 ifelse(cause.df$gs_text34 %in% other, "Other",
                                        ifelse(cause.df$gs_text34 %in% other_infectious, "Other Infectious",
                                               cause.df$gs_text34)))
    }
    
    
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
   # ndraws <- 50E3
    ndraws <- 30E3
    nchains <- 3
    
    ### calibva with Tariff
    set.seed(123)
    calibva.seeds <- sample(1e6, nchains, replace = F)
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
    ensemble.calibva <- window(ensemble.calibva, start = burnin, thin = thin)
    ensemble.calibva.csmf.df <- calibvaCSMFPosteriorSamples(ensemble.calibva, causes = causes) %>%
        mutate(ncalib = setting$calib.size)
    return(ensemble.calibva.csmf.df)
})

all_posteriors <- do.call(bind_rows, posteriors)
levels(all_posteriors$cause) <- c("Pneumonia",
                                  "Diarrhea/Dysentery",
                                  "Malaria",
                                  "Other")

posteriors_plot <-
    ggplot(all_posteriors, aes(x = value, color= factor(ncalib))) +
    geom_line(stat = 'density') +
    xlim(0, 1) +
    facet_wrap(~cause, nrow = 1, scales = 'free_y') +
    xlab("CSMF") +
    guides(color=guide_legend(title="|H|")) +
    theme(legend.position="bottom")
  
fig.dir <- here("CSMFPredictionsWeightedResampling", "figs")
ggsave(file.path(fig.dir, "posteriors.pdf"),
       posteriors_plot, width = 12, height = 6)


#########################################################################
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
setting.df <- expand.grid(run = 2,
                          country = c("India"),
                          calib.size = c(50, 100, 200),
                          cod = c("broad"),
                          calib.csfma = c("low"))

posteriors <- lapply(1:nrow(setting.df), function(i) {
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
    if(setting$cod == "top3") {
        top_causes <- names(sort(table(country.data$Cause), decreasing = TRUE))[1:3]
        cause.df$broad <- ifelse(cause.df$gs_text34 %in% top_causes,
                                 as.character(cause.df$gs_text34), "Other")
    } else {
        cause.df$broad <- ifelse(cause.df$gs_text34 %in% external, "External",
                                 ifelse(cause.df$gs_text34 %in% other, "Other",
                                        ifelse(cause.df$gs_text34 %in% other_infectious, "Other Infectious",
                                               cause.df$gs_text34)))
    }
    
    
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
    beta <- .5
    tau <- .5
    epsilon <- .001
    #epsilon <- .01
    tau.vec <- rep(tau, length(causes))
    delta <- 1
    gamma.init <- 1
    # ndraws <- 50E3
    ndraws <- 30E3
    nchains <- 3
    #nchains <- 1
    
    set.seed(123)
    calibva.seeds <- sample(1e6, nchains, replace = F)
    insilico_list <- lapply(c(5, 10, 15), function(a) {
        insilico.calibva <- calibva.sampler(test.cod = insilico.train.cod.test,
                                            calib.cod = insilico.train.cod.calib,
                                            calib.truth = calib.truth, causes = causes,
                                            epsilon = epsilon, alpha=a, beta=beta,
                                            tau.vec=tau.vec, delta=delta,
                                            gamma.init=gamma.init, ndraws = ndraws,
                                            nchains = nchains,
                                            init.seed = calibva.seeds)
        
        ### Get posterior draws for CSMF parameters
        burnin <- 10E3
        thin <- 10
        insilico.calibva <- window(insilico.calibva, start = burnin, thin = thin)
        insilico.calibva.csmf.df <- calibvaCSMFPosteriorSamples(insilico.calibva, causes = causes) %>%
            mutate(ncalib = setting$calib.size,
                   alpha = a)
        return(insilico.calibva.csmf.df)
    })
    insilico.calibva.csmf.df <- do.call(rbind, insilico_list)
    pinsilico <- sapply(top.cause.df$cause, function(c) mean(insilico.train.cod.test == c))
    top.cause.df$pinsilico <- pinsilico
    insilico.calibva.csmf.df <- inner_join(insilico.calibva.csmf.df,
                                           top.cause.df,
                                           by = "cause")
    return(cbind(insilico.calibva.csmf.df, setting))
})

all_posteriors <- do.call(bind_rows, posteriors)
all_posteriors$cause <- factor(all_posteriors$cause,
                               levels = c("Pneumonia",
                                  "Diarrhea/Dysentery",
                                  "External",
                                  "Other Infectious",
                                  "Other"))


plots <- lapply(c(50, 100, 200), function(ncal) {
    p <-
        all_posteriors %>%
        filter(ncalib == ncal) %>%
        ggplot(aes(x = value, color = factor(alpha))) +
        geom_line(stat = 'density') +
        #facet_wrap(~cause, nrow = 1, scales = 'free_y') +
        xlab("CSMF") +
        facet_wrap( ~ cause, scales = "free_y", nrow = 1) +
        #geom_vline(aes(xintercept = ptest, col = 'True CSMF')) +
        geom_vline(aes(xintercept = pinsilico)) +
        xlim(0, 1) +
        ylab(paste0("NCalib ", ncal)) +
        guides(color = guide_legend(title = "Alpha"))
    return(p)
})
library(ggpubr)
insilicio_prior_plot <- ggarrange(plotlist = plots, nrow = 3, common.legend = TRUE)
ggsave(file.path(fig.dir, "InsilicoPriorComparison.pdf"),
       plot = insilicio_prior_plot, width = 10, height = 10)

fig.dir <- here("CSMFPredictionsWeightedResampling", "figs")
ggsave(file.path(fig.dir, "posteriors.pdf"),
       posteriors_plot, width = 12, height = 6)







