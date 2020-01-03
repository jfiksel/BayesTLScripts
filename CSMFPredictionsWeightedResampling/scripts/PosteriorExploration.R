library(here)
library(openVA)
library(CalibratedVA)
library(tidyverse)
library(coda)
library(ggmcmc)
library(gtools)
library(tools)
library(ggpubr)
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




### Get distribution of these causes in the train set
### first, function to change all guesses not equal to top 3 causes to other
changeTopCOD <- function(topcod, cause.df) {
    topcod <- as.character(topcod)
    m <- match(topcod, cause.df$gs_text34)
    return(cause.df$broad[m])
}

### Function to give posterior distribution of M matrix and CSMF
get_posteriors <- function(setting, calib_ind) {
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
                                    ifelse(cause.df$gs_text34 %in% other_infectious,
                                           "Other Infectious",
                                           cause.df$gs_text34)))
    if(setting$cod == 7) {
        cause.df$broad[cause.df$gs_text34 == "Malaria"] <- "Malaria"
        cause.df$broad[cause.df$gs_text34 == "Sepsis"] <- "Sepsis"
    }
    
    top.cause.df <- data.frame(cause = unique(cause.df$broad))
    causes <- as.character(top.cause.df$cause)
    true_causes <- changeTopCOD(child.clean[countries == setting$country,]$Cause,
                                cause.df = cause.df)
    insilico_file <- here("CSMFPredictionsWeightedResampling", "child_data",
                          paste0("insilico_", setting$country, ".rds"))
    insilico_cod <- readRDS(insilico_file)
    insilico.train.cod <- cause.df$broad[match(insilico_cod, cause.df$va34)]
    
    # Get "true"  misclassification matrix 
    causes <- top.cause.df$cause
    true_t <- rawMisclassificationMatrix(insilico.train.cod, true_causes, causes)
    true_m <- normalizedMisclassificationMatrix(true_t)
    
    
    ### Calibration
    set.seed(all.seeds[setting$run])
    ### split data into country we are using and other
    country.data <- child.clean[countries == setting$country,]
    train.final <- child.clean[countries != setting$country,]
    
    
    top.cause.df$ptrain <- sapply(top.cause.df$cause, function(c) mean(changeTopCOD(train.final$Cause, cause.df) == c))
    top.cause.df$ptest <- sapply(top.cause.df$cause, function(c) mean(changeTopCOD(country.data$Cause, cause.df) == c))
    top.cause.df$ntest <- sapply(top.cause.df$cause, function(c) sum(changeTopCOD(country.data$Cause, cause.df) == c))
    


    test.final <- country.data
    calib.final <- country.data[calib_ind,]
    
    
    ### Get causes (these will be the ordering we will use)
    causes <- as.character(top.cause.df$cause)
    calib.truth <- changeTopCOD(calib.final$Cause, cause.df)
    top.cause.df$pcalib <- sapply(top.cause.df$cause, function(c) mean(calib.truth == c))
    
    
    ### Insilico with just training set
    insilico_file <- here("CSMFPredictionsWeightedResampling", "child_data", paste0("insilico_", setting$country, ".rds"))
    insilico_cod <- readRDS(insilico_file)
    insilico.train.cod <- cause.df$broad[match(insilico_cod, cause.df$va34)]
    
    insilico.train.cod.test <- insilico.train.cod
    insilico.train.cod.calib <- insilico.train.cod[calib_ind]
    
    top.cause.df$q <- sapply(top.cause.df$cause, function(c) mean(insilico.train.cod.test == c))
    
    
    alpha <- setting$alpha
    beta <- .5
    tau <- .5
    epsilon <- .001
    tau.vec <- rep(tau, length(causes))
    delta <- 1
    gamma.init <- 1
    ndraws <- 30E3
    nchains <- 3
    
    ### CalibratedVA with InSilicoVA
    set.seed(123)
    calibva.seeds <- sample(1e6, nchains, replace = F)
    insilico.calibva <- calibva.sampler(test.cod = insilico.train.cod.test,
                                        calib.cod = insilico.train.cod.calib,
                                        calib.truth = calib.truth, causes = causes,
                                        epsilon = epsilon, alpha=alpha, beta=beta,
                                        tau.vec=tau.vec, delta=delta,
                                        gamma.init=gamma.init, ndraws = ndraws,
                                        nchains = nchains,
                                        init.seed = calibva.seeds)
    burnin <- 10E3
    thin <- 10
    insilico.calibva <- window(insilico.calibva, start = burnin, thin = thin)
    
    ### Posterior samples of M-matrix
    insilico_m <- ggs(insilico.calibva, family = "M")
    
    insilico_m_mean <-
        insilico_m  %>%
        group_by(Parameter) %>%
        summarize(mean = mean(value))
    
    insilico_m_mean$true_value <- rep(NA, nrow(insilico_m_mean))
    insilico_m_mean$cause <- rep(NA, nrow(insilico_m_mean))
    for(i in 1:nrow(insilico_m_mean)) {
        row <- as.numeric(substr(insilico_m_mean$Parameter[i], 3, 3))
        col <- as.numeric(substr(insilico_m_mean$Parameter[i], 5, 5))
        insilico_m_mean$true_value[i] <- true_m[causes[row], causes[col]]
        insilico_m_mean$cause[i] <- paste0("M[", causes[row], ",", causes[col], "]")
    }
    
    
    ### Now look at estimation of CSMFs
    insilico.calibva.csmf.df <- calibvaCSMFPosteriorSamples(insilico.calibva, causes = causes)
    
    insilico_csmf_mean <-
        insilico.calibva.csmf.df %>%
        group_by(cause) %>%
        summarize(posterior_mean = mean(value)) 
    top.cause.df <- inner_join(top.cause.df, insilico_csmf_mean, by = "cause")
    setting <- select(setting, country, calib.size, cod, alpha, calib.csfma)
    return(list(setting = setting,
                m_posterior = insilico_m,
                m_df = insilico_m_mean,
                csmf_posterior = insilico.calibva.csmf.df,
                csmf_df =  top.cause.df ))
}


### Read in list with settings & calibration indices
calib_list <- readRDS(here("CSMFPredictionsWeightedResampling", "child_data",
                           "calib_indices_for_eda.rds"))



india_settings <- calib_list[1:2]
### Add setting where we use everyone for calibration
india_settings[[3]] <- india_settings[[2]]
india_settings[[3]]$setting$calib.size <- sum(countries == "India")
india_settings[[3]]$setting$calib.csmfa <- "All" 
india_settings[[3]]$calib_ind <- 1:sum(countries == "India")

####### Just get one posterior for India
india_insilico_posterior <-  get_posteriors(india_settings[[1]]$setting, india_settings[[1]]$calib_ind)

csmf_posterior <- india_insilico_posterior$csmf_posterior
fig.dir <- here("CSMFPredictionsWeightedResampling", "figs")
posterior_plot <-
    csmf_posterior %>%
    ggplot(aes(x = value)) +
    geom_line(stat = 'density') +
    facet_wrap(~cause, scales = "free_y", nrow = 1) +
    xlim(0,1)
ggsave(file.path(fig.dir, "insilico_posteriors_200_example.pdf"),
       posterior_plot, width = 8, height = 5)



### First explore India, low vs representative calibration set  

india_posteriors <- lapply(seq_along(india_settings), function(i) {
    get_posteriors(india_settings[[i]]$setting, india_settings[[i]]$calib_ind)
})

### CSMF plots
csmf_plots <- lapply(seq_along(india_posteriors), function(i) {
    posterior <- india_posteriors[[i]]
    title <- paste("Calibration Set", toTitleCase(as.character(posterior$setting$calib.csfma)), "CSMFA")
    if(i == 3) {
        title <- "Calibration set all data"
    }
    plot <-
        ggplot(posterior$csmf_posterior, aes(x = value)) +
        geom_line(stat = 'density') +
        geom_vline(data = posterior$csmf_df, aes(xintercept = ptest, colour = 'Truth')) +
        geom_vline(data = posterior$csmf_df, aes(xintercept = q, colour = 'Raw')) +
        geom_vline(data = posterior$csmf_df, aes(xintercept = posterior_mean, colour = 'Posterior Mean')) +
        facet_wrap(~cause, nrow = 5,
                   scales = "free_y") +
        scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, by = .2)) +
        scale_colour_manual(name = '', 
                            values =c('Truth'='red', 'Raw' = 'green', 'Posterior Mean'='blue')) +
        theme(legend.position = 'bottom') +
        ggtitle(title)
    return(plot)
})

csmf_comparison_plot <- ggarrange(plotlist = csmf_plots, ncol = 3, common.legend = TRUE)
fig.dir <- here("CSMFPredictionsWeightedResampling", "figs")
ggsave(file.path(fig.dir, "insilico_rep_vs_low_csmf.pdf"),
       csmf_comparison_plot, width = 10, height = 10)


### Now do M matrix plots
insilico_m_df <- do.call(rbind, lapply(seq_along(india_posteriors), function(i) {
    posterior <- india_posteriors[[i]]
    m_posterior <- posterior$m_posterior
    m_df <- posterior$m_df
    m_posterior$cause <- m_df$cause[match(m_posterior$Parameter, m_df$Parameter)]
    m_posterior$calib.csmfa <- posterior$setting$calib.csfma
    if(i == 3) {
        m_posterior$calib.csmfa <- "All"  
    }
    return(m_posterior)
}))



m_comparison_plot <-
    ggplot(insilico_m_df, aes(x = value, color = calib.csmfa)) +
    geom_line(stat = 'density') +
    geom_vline(data = india_posteriors[[1]]$m_df, aes(xintercept = true_value)) +
    facet_wrap(~cause, scales = "free_y") +
    theme(legend.position = "bottom")
ggsave(file.path(fig.dir, "insilico_rep_vs_low_m_matrices_combined.pdf"),
       m_comparison_plot, width = 14, height = 10)

############################

### Now compare Tanzania posteriors for # cod
### Need new function because we're comparing ensemble
get_posteriors_ensemble <- function(setting, calib_ind) {
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
                                    ifelse(cause.df$gs_text34 %in% other_infectious,
                                           "Other Infectious",
                                           cause.df$gs_text34)))
    if(setting$cod == 7) {
        cause.df$broad[cause.df$gs_text34 == "Malaria"] <- "Malaria"
        cause.df$broad[cause.df$gs_text34 == "Sepsis"] <- "Sepsis"
    }
    
    top.cause.df <- data.frame(cause = unique(cause.df$broad))
    causes <- as.character(top.cause.df$cause)
    true_causes <- changeTopCOD(child.clean[countries == setting$country,]$Cause,
                                cause.df = cause.df)
    insilico_file <- here("CSMFPredictionsWeightedResampling", "child_data",
                          paste0("insilico_", setting$country, ".rds"))
    insilico_cod <- readRDS(insilico_file)
    insilico.train.cod <- cause.df$broad[match(insilico_cod, cause.df$va34)]
    
    # Get "true"  misclassification matrix 
    causes <- top.cause.df$cause
    true_t <- rawMisclassificationMatrix(insilico.train.cod, true_causes, causes)
    true_m <- normalizedMisclassificationMatrix(true_t)
    
    
    ### Calibration
    set.seed(all.seeds[setting$run])
    ### split data into country we are using and other
    country.data <- child.clean[countries == setting$country,]
    train.final <- child.clean[countries != setting$country,]
    
    
    top.cause.df$ptrain <- sapply(top.cause.df$cause, function(c) mean(changeTopCOD(train.final$Cause, cause.df) == c))
    top.cause.df$ptest <- sapply(top.cause.df$cause, function(c) mean(changeTopCOD(country.data$Cause, cause.df) == c))
    top.cause.df$ntest <- sapply(top.cause.df$cause, function(c) sum(changeTopCOD(country.data$Cause, cause.df) == c))
    
    
    
    test.final <- country.data
    calib.final <- country.data[calib_ind,]
    
    
    ### Get causes (these will be the ordering we will use)
    causes <- as.character(top.cause.df$cause)
    calib.truth <- changeTopCOD(calib.final$Cause, cause.df)
    top.cause.df$pcalib <- sapply(top.cause.df$cause, function(c) mean(calib.truth == c))
    
    
    
    ### Insilico 
    insilico_file <- here("CSMFPredictionsWeightedResampling", "child_data", paste0("insilico_", setting$country, ".rds"))
    insilico_cod <- readRDS(insilico_file)
    insilico.train.cod <- cause.df$broad[match(insilico_cod, cause.df$va34)]
    
    insilico.train.cod.test <- insilico.train.cod
    insilico.train.cod.calib <- insilico.train.cod[calib_ind]
    
    ### Tariff
    tariff_file <- here("CSMFPredictionsWeightedResampling",
                        "child_data", paste0("tariff_", setting$country, ".rds"))
    tariff_cod <- readRDS(tariff_file)
    tariff.train.cod <- cause.df$broad[match(tariff_cod, cause.df$va34)]
    tariff.train.cod.test <- tariff.train.cod
    tariff.train.cod.calib <- tariff.train.cod[calib_ind]
    
    top.cause.df$qinsilico <- sapply(top.cause.df$cause, function(c) mean(insilico.train.cod.test == c))
    top.cause.df$qtariff <- sapply(top.cause.df$cause, function(c) mean(tariff.train.cod.test == c))
    
    alpha <- setting$alpha
    beta <- .5
    tau <- .5
    epsilon <- .001
    tau.vec <- rep(tau, length(causes))
    delta <- 1
    gamma.init <- 1
    ndraws <- 30E3
    nchains <- 3
    
    ### CalibratedVA with Ensemble
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
    burnin <- 10E3
    thin <- 10
    ensemble.calibva <- window(ensemble.calibva, start = burnin, thin = thin)
    
    ### Get M matrix posterior
    m_matrices <- ggs(ensemble.calibva, family = "M")
    
    ### Change column names
    m_matrix_names <- data.frame(Parameter = unique(m_matrices$Parameter))
    m_matrix_names$Algorithm <- rep(NA, nrow(m_matrix_names))
    m_matrix_names$cause <- rep(NA, nrow(m_matrix_names))
    
    for(i in 1:nrow(m_matrix_names)) {
        row <- as.numeric(substr(m_matrix_names$Parameter[i], 3, 3))
        col <- as.numeric(substr(m_matrix_names$Parameter[i], 5, 5))
        alg <- as.numeric(substr(m_matrix_names$Parameter[i], 7, 7))
        m_matrix_names$Algorithm[i] <- ifelse(alg == 1, "Tariff", "InSilico")
        m_matrix_names$cause[i] <-  paste0("M[", causes[row], ",", causes[col], "]")
    }
    m_matrices$cause <- m_matrix_names$cause[match(m_matrices$Parameter, m_matrix_names$Parameter)]
    m_matrices$algorithm <- m_matrix_names$Algorithm[match(m_matrices$Parameter, m_matrix_names$Parameter)]
    ### Now look at estimation of CSMFs
    ensemble.calibva.csmf.df <- calibvaCSMFPosteriorSamples(ensemble.calibva, causes = causes)
    
    ensemble_csmf_mean <-
        ensemble.calibva.csmf.df %>%
        group_by(cause) %>%
        summarize(posterior_mean = mean(value)) 
    top.cause.df <- inner_join(top.cause.df, ensemble_csmf_mean, by = "cause")
    setting <- select(setting, country, calib.size, cod, alpha, calib.csfma)
    return(list(setting = setting,
                csmf_posterior = ensemble.calibva.csmf.df,
                csmf_df =  top.cause.df,
                m_matrix_posterior = m_matrices))
}


tanzania_settings <- calib_list[3:6]
tanzania_posteriors <- lapply(seq_along(tanzania_settings), function(i) {
    posteriors <- get_posteriors_ensemble(tanzania_settings[[i]]$setting, tanzania_settings[[i]]$calib_ind)
    ### Adjust estimates for COD = 7
    csmf_posterior <-  posteriors$csmf_posterior
    csmf_posterior$cause <- as.character(csmf_posterior$cause)
    csmf_posterior$cause[csmf_posterior$cause %in% c('Malaria', 'Sepsis')] <- "Other Infectious"
    csmf_posterior <-
        csmf_posterior %>%
        group_by(Iteration, Chain, cause) %>%
        summarize(value = sum(value))
    posteriors$csmf_posterior <- csmf_posterior
    csmf_df <- posteriors$csmf_df
    csmf_df$cause <- as.character(csmf_df$cause)
    csmf_df$cause[csmf_df$cause %in% c('Malaria', 'Sepsis')] <- "Other Infectious"
    csmf_df <-
        csmf_df %>%
        group_by(cause) %>%
        summarize(ptest = sum(ptest),
                  qinsilico = sum(qinsilico),
                  qtariff = sum(qtariff),
                  posterior_mean = sum(posterior_mean))
    posteriors$csmf_df <- csmf_df
    return(posteriors)
})

tanzania_df <- do.call(rbind, lapply(seq_along(tanzania_posteriors), function(i) {
    csmf_posterior <-  tanzania_posteriors[[i]]$csmf_posterior
    csmf_posterior$cod <-  tanzania_posteriors[[i]]$setting$cod
    csmf_posterior$calib.size <- tanzania_posteriors[[i]]$setting$calib.size
    return(csmf_posterior)
}))

tanzania_csmf_df <- do.call(rbind, lapply(seq_along(tanzania_posteriors), function(i) {
    csmf_df <-  tanzania_posteriors[[i]]$csmf_df
    csmf_df$cod <-  tanzania_posteriors[[i]]$setting$cod
    csmf_df$calib.size <- tanzania_posteriors[[i]]$setting$calib.size
    return(csmf_df)
}))

tanzania_cod_plot <-
    ggplot(tanzania_df, aes(x = value, linetype = factor(cod))) +
    geom_line(stat = 'density') +
    facet_grid(calib.size ~ cause) +
    geom_vline(data = tanzania_csmf_df, aes(xintercept = qtariff, color = 'Tariff Raw')) +
    geom_vline(data = tanzania_csmf_df, aes(xintercept = qinsilico, color = 'InSilico Raw')) +
    geom_vline(data = tanzania_csmf_df, aes(xintercept = ptest, color = 'Truth')) +
    scale_colour_manual(name = '', 
                        values =c('Tariff Raw'='blue', 'InSilico Raw' = 'green', 'Truth' = 'red'))  +
    labs(linetype = "Number COD") +
    theme(legend.position = 'bottom')
ggsave(file.path(fig.dir, "tanzania_cod_5_vs_7_posterior.pdf"),
       tanzania_cod_plot, width = 8, height = 8)

### Get m matrix  posteriors
tanzania_m_df <- do.call(rbind, lapply(seq_along(tanzania_posteriors), function(i) {
    m_matrix <- tanzania_posteriors[[i]]$m_matrix_posterior
    m_matrix$cod <-  tanzania_posteriors[[i]]$setting$cod
    m_matrix$calib.size <- tanzania_posteriors[[i]]$setting$calib.size
    return(m_matrix)
}))

setting.df <- expand.grid(cod = c(5,7), algorithm = c("Tariff", "InSilico"))
pdf(file.path(fig.dir, "tanzania_cod_5_vs_7_m_matrices.pdf"), width = 14, height = 10)
for(i in 1:nrow(setting.df)) {
    setting <- setting.df[i,]
    title <- paste(setting$cod, "COD", setting$algorithm)
    plot <-
        tanzania_m_df %>%
        filter(cod == setting$cod, algorithm == setting$algorithm) %>%
        ggplot(aes(x = value, color = factor(calib.size))) +
        geom_line(stat = 'density') +
        facet_wrap(~cause, scales = 'free_y') +
        ggtitle(title)+
        labs(color = "Calibration size") +
        theme(legend.position = 'bottom')
    print(plot)
}
dev.off()








