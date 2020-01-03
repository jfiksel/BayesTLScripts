library(openVA)
library(here)
set.seed(123)
### Read in child data
child.raw <- read.csv(getPHMRC_url("child"))
### Get data frame matching va34 code to actual COD
cause.df <- unique(child.raw[,c("gs_text34", "va34")])
cause.df$va34 <- as.character(as.numeric(cause.df$va34))
### Clean data into output usable by Tariff & InsilicoVA
child.clean <- ConvertData.phmrc(child.raw, phmrc.type = "child")$output
### Assign countries
countries <- rep(NA, nrow(child.raw))
countries[child.raw$site %in% c("AP", "UP")] <- "India"
countries[child.raw$site %in% c("Mexico")] <- "Mexico"
countries[child.raw$site %in% c("Dar", "Pemba")] <- "Tanzania"
countries[child.raw$site %in% c("Bohol")] <- "Philippines"

### Create data directory
data_dir <- here("CSMFPredictionsWeightedResampling", "child_data")
if(!dir.exists(data_dir)){
    dir.create(data_dir)
}

for(country in c('India', "Tanzania")){
    ### split data into country we are using and other
    country.data <- child.clean[countries == country,]
    train.final <- child.clean[countries != country,]

    ################################
    ### Tariff with just training set (will also predict on calibration data)
    
    set.seed(123)
    tariff.train <- codeVA(data = country.data,
                           data.type = "customize", model = "Tariff",
                           data.train = train.final, causes.train = "Cause")
    tariff.train.cod <- getTopCOD(tariff.train)[,2]
    saveRDS(tariff.train.cod, file.path(data_dir, paste0("tariff_", country, ".rds")))
   
    #### Now InsilicoVA
    ### Set number of simulations for insilico
    insilico.nsim <- 5000
    set.seed(123)
    insilico.train <- codeVA(data = country.data,
                             data.type = "customize", model = "InSilicoVA",
                             data.train = train.final, causes.train = "Cause",
                             jump.scale = 0.05, Nsim=insilico.nsim, auto.length = FALSE)
    insilico.train.cod <- getTopCOD(insilico.train)[,2]
    saveRDS(insilico.train.cod, file.path(data_dir, paste0("insilico_", country, ".rds")))
}
