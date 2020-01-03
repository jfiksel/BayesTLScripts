### Plot true CSMF for different cause categories
library(tidyverse)
library(openVA)

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


### Get distribution of these causes in the train set
### first, function to change all guesses not equal to top 3 causes to other
changeTopCOD <- function(topcod) {
    topcod <- as.character(topcod)
    m <- match(topcod, cause.df$gs_text34)
    return(cause.df$broad[m])
}
fig.dir <- here("CSMFPredictionsWeightedResampling", "figs")

pdf(file.path(fig.dir, "CSMF_by_country.pdf"), width = 8, height = 8)
for(cod in c('broad5', 'broad7')) {
    for(country in c('India', 'Tanzania')) {
        country.data <- child.clean[countries == country,]
        train.final <- child.clean[countries != country,]
        if(cod == "broad5"){
            cause.df$broad <- ifelse(cause.df$gs_text34 %in% external, "External",
                                     ifelse(cause.df$gs_text34 %in% other, "Other",
                                            ifelse(cause.df$gs_text34 %in% other_infectious, "Other Infectious",
                                                   cause.df$gs_text34)))
        } else {
            cause.df$broad <- ifelse(cause.df$gs_text34 %in% external, "External",
                                     ifelse(cause.df$gs_text34 %in% other, "Other",
                                            ifelse(cause.df$gs_text34 %in% other_infectious, "Other Infectious",
                                                   cause.df$gs_text34)))
            cause.df$broad[cause.df$gs_text34 %in% c("Malaria", "Sepsis")] <- cause.df$gs_text34[cause.df$gs_text34 %in% c("Malaria", "Sepsis")]
        }
        causes <- changeTopCOD(as.character(country.data$Cause))
        train.causes <- changeTopCOD(as.character(train.final$Cause))
        
        insilico_file <- here("CSMFPredictionsWeightedResampling", "child_data", paste0("insilico_", country, ".rds"))
        insilico_cod <- readRDS(insilico_file)
        insilico.train.cod <- cause.df$broad[match(insilico_cod, cause.df$va34)]
        
        tariff_file <- here("CSMFPredictionsWeightedResampling", "child_data", paste0("tariff_", country, ".rds"))
        tariff_cod <- readRDS(tariff_file)
        tariff.train.cod <- cause.df$broad[match(tariff_cod, cause.df$va34)]
        
        if(cod == "broad5") {
            cause_order <- c("Pneumonia", "Diarrhea/Dysentery", "External", "Other Infectious",
                             "Other") 
        } else {
            cause_order <- c("Pneumonia", "Diarrhea/Dysentery", "External",
                             "Malaria", "Sepsis", "Other Infectious",
                             "Other") 
        }
        
        csmf_true <- sapply(cause_order, function(c) mean(causes == c))
        csmf_train <- sapply(cause_order, function(c) mean(train.causes == c))
        csmf_insilico <- sapply(cause_order, function(c) mean(insilico.train.cod == c))
        csmf_tariff <- sapply(cause_order, function(c) mean(tariff.train.cod == c))
        csmf_all <- c(csmf_true, csmf_train, csmf_insilico, csmf_tariff)
        plot_df <- data.frame(cause = names(csmf_all), csmf = unname(as.vector(csmf_all)),
                              method = rep(c("CSMF_P", "CSMF_G",
                                             "CSMF_Insilico", "CSMF_Tariff"),
                                           each = length(cause_order)))
        plot_df$cause = factor(plot_df$cause, levels = cause_order)
        levels(plot_df$method) <- c("CSMF_P", "CSMF_G",
                                    "CSMF_Insilico", "CSMF_Tariff")
        
        title <- paste0("CSMF ", country)
        csmf_plot <-
            ggplot(plot_df, aes(x = cause, y = csmf, fill = method)) +
            geom_bar(stat = 'identity', position = "dodge") +
            ylab("CSMF") +
            ylim(0,.6) +
            xlab("Cause") +
            theme(legend.title = element_blank(),
                  legend.position = "bottom") +
            ggtitle(title) 
        print(csmf_plot)
    }
}
dev.off()


#### True M Matrices for India & Tanzania with InSilico 
m_mats <- lapply(c("India", "Tanzania"), function(country) {
    insilico_file <- here("CSMFPredictionsWeightedResampling", "child_data", paste0("insilico_", country, ".rds"))
    insilico_cod <- readRDS(insilico_file)
    cause.df$broad <- ifelse(cause.df$gs_text34 %in% external, "External",
                             ifelse(cause.df$gs_text34 %in% other, "Other",
                                    ifelse(cause.df$gs_text34 %in% other_infectious, "Other Infectious",
                                           cause.df$gs_text34)))
    insilico.train.cod <- cause.df$broad[match(insilico_cod, cause.df$va34)]
    country.data <- child.clean[countries == country,]
    causes_truth <- changeTopCOD(as.character(country.data$Cause))
    causes <- unique(cause.df$broad)
    T_mat <- rawMisclassificationMatrix(insilico.train.cod, causes_truth, causes)
    M_mat <- normalizedMisclassificationMatrix(T_mat)
    return(M_mat)
})

### Different shrinkage

