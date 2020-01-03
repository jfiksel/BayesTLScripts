library(openVA)
library(here)
models_dir <- here("CSMFPredictionsWeightedResampling", "adult_data")
if(!dir.exists(models_dir)){
    dir.create(models_dir, recursive = TRUE)
}
if(!file.exists(here("CSMFPredictionsWeightedResampling", "adult_data", "phmrc_adult.rds"))){
    phmrc <- read.csv(getPHMRC_url("adult"))
    saveRDS(phmrc, here("CSMFPredictionsWeightedResampling", "adult_data", "phmrc_adult.rds"))
} else {
    phmrc <- readRDS(here("CSMFPredictionsWeightedResampling", "adult_data", "phmrc_adult.rds"))
}

country.df <- data.frame(site = c("AP", "Bohol", "Dar", "Mexico", "Pemba", "UP"),
                         country = c("India", 
                                     "Philippines",
                                     "Tanzania",
                                     "Mexico",
                                     "Tanzania",
                                     "India"))
country <- country.df$country[match(phmrc$site, country.df$site)]

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

cause_map <- data.frame(causes = c(external, circulatory, non_communicable,
                                   infectious, maternal),
                        broad_cause = rep(c("external", "circulatory",
                                            "non_communicable", "infectious", "maternal"),
                                          c(length(external), length(circulatory),
                                            length(non_communicable), length(infectious),
                                            length(maternal))))

### Train for each country
id <- as.numeric(commandArgs(trailingOnly = TRUE))
c <- unique(as.character(country))[id]
test <- phmrc[country == c,]
train <- phmrc[country != c,]
if(!file.exists(file.path(models_dir, paste0("insilico_model_", c, ".rds")))){
    set.seed(123)
    phmrc.insilicova <- codeVA(data = test, data.type = "PHMRC",
                               model = "InSilicoVA",
                               data.train = train, causes.train = "gs_text34",
                               phmrc.type = "adult",
                               jump.scale = 0.05, convert.type = "fixed",
                               Nsim=10000, auto.length = FALSE)
    saveRDS(phmrc.insilicova, file.path(models_dir, paste0("insilico_model_", c, ".rds")))
} else {
    phmrc.insilicova <- readRDS(file.path(models_dir, paste0("insilico_model_", c, ".rds")))
}
if(!file.exists(file.path(models_dir, paste0("tariff_model_", c, ".rds")))){
    set.seed(123)
    tariff <- codeVA(data = test, data.type = "PHMRC", model = "Tariff",
                     data.train = train, causes.train = "gs_text34",
                     phmrc.type = "adult")
    saveRDS(tariff, file.path(models_dir, paste0("tariff_model_", c, ".rds")))
} else {
    tariff <- readRDS(file.path(models_dir, paste0("tariff_model_", c, ".rds")))
} 
if(!file.exists(file.path(models_dir, paste0("interva_model_", c, ".rds")))){
    set.seed(123)
    interva <- codeVA(data = test, data.type = "PHMRC", model = "InterVA",
                      data.train = train, causes.train = "gs_text34",
                      phmrc.type = "adult")
    saveRDS(interva, file.path(models_dir, paste0("interva_model_", c, ".rds")))
} else {
    interva <- readRDS(file.path(models_dir, paste0("interva_model_", c, ".rds")))
} 
if(!file.exists(file.path(models_dir, paste0("nbc_model_", c, ".rds")))){
    set.seed(123)
    nbc <- codeVA(data = test, data.type = "PHMRC", model = "NBC",
                  data.train = train, causes.train = "gs_text34",
                  phmrc.type = "adult")
    saveRDS(nbc, file.path(models_dir, paste0("nbc_model_", c, ".rds")))
} else {
    nbc <- readRDS(file.path(models_dir, paste0("nbc_model_", c, ".rds")))
} 

### Top COD
### InSilicoVA
insilico_cod <- getTopCOD(phmrc.insilicova)[,2]
insilico_broad_cod <- cause_map$broad_cause[match(insilico_cod, cause_map$causes)]
saveRDS(insilico_broad_cod, file.path(models_dir, paste0("insilico_model_", c, "_cod.rds")))

### Tariff 
tariff_cod <- getTopCOD(tariff)[,2]
tariff_broad_cod <- cause_map$broad_cause[match(tariff_cod, cause_map$causes)]
saveRDS(tariff_broad_cod, file.path(models_dir, paste0("tariff_model_", c, "_cod.rds")))

### InterVA 
interva_cod <- getTopCOD(interva)[,2]
interva_broad_cod <- cause_map$broad_cause[match(interva_cod, cause_map$causes)]
saveRDS(interva_broad_cod, file.path(models_dir, paste0("interva_model_", c, "_cod.rds")))


### NBC 
nbc_cod <- getTopCOD(nbc)[,2]
nbc_broad_cod <- cause_map$broad_cause[match(nbc_cod, cause_map$causes)]
saveRDS(nbc_broad_cod, file.path(models_dir, paste0("nbc_model_", c, "_cod.rds")))
quit('no')

