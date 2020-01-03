library(tidyverse)
library(here)
library(openVA)

### First CSMFA between G and P

set.seed(123)
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

changeTopCOD <- function(topcod) {
    topcod <- as.character(topcod)
    m <- match(topcod, cause.df$gs_text34)
    return(cause.df$broad[m])
}

### get causes
top.cause.df <- data.frame(cause = unique(cause.df$broad))

for(country in c('India', 'Tanzania')) {
    country.data <- child.clean[countries == country,]
    train.final <- child.clean[countries != country,]
    top.cause.df$ptrain <- sapply(top.cause.df$cause, function(c) mean(changeTopCOD(train.final$Cause) == c))
    top.cause.df$ptest <- sapply(top.cause.df$cause, function(c) mean(changeTopCOD(country.data$Cause) == c)) 
    print(getCSMF_accuracy(top.cause.df$ptrain, top.cause.df$ptest))
}

#############
### CSMFA using just calibration set
child_results <- readRDS(here("CSMFPredictionsWeightedResampling", "child_data", "csmfa_results.rds"))
child_results %>%
    filter(method == "calib", cod == 5, calib.csfma == "low") %>%
    group_by(country, calib.size) %>%
    summarize(mean_csmfa = mean(csmfa))

### CSMFA when using calibration set to train
child_results %>%
    filter(grepl("_and_calib", method), cod == 5, calib.csfma == "low") %>%
    group_by(method, country, calib.size) %>%
    summarize(mean_csmfa = mean(csmfa))
