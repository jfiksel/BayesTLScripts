library(here)
library(openVA)
library(CalibratedVA)
library(tidyverse)
library(coda)
library(ggmcmc)
library(gtools)


fig.dir <- here("CSMFPredictionsWeightedResampling", "figs")
if(!dir.exists(fig.dir)) {
    dir.create(fig.dir, recursive = TRUE)
}
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
setting.df <- expand.grid(country = c("India", "Tanzania"),
                          calib.size = c(200),
                          cod = c(5),
                          calib.csfma = c("low"))

for(i in 1:nrow(setting.df)) {
    setting <- setting.df[i,]
    ptestsymp <- c()
    pcalibsymp <- c()
    for(run in seq(1, 500, by = 50)) {
        ### split data into country we are using and other
        country.data <- child.clean[countries == setting$country,]
        train.final <- child.clean[countries != setting$country,]
        calib.dir <- here("CSMFPredictionsWeightedResampling", "cluster_output_child", "calib_indices")
        calib.file <- file.path(calib.dir,
                                paste0(setting$country,
                                       "run", run,
                                       "calibsize", setting$calib.size,
                                       "csfma", setting$calib.csfma, ".rds"))
        calib_ind <- readRDS(calib.file)
        
        test.final <- country.data[,-(1:2)]
        calib.final <- country.data[calib_ind,-(1:2)]
        ptestsymp <- c(ptestsymp, colMeans(test.final == "Y"))
        pcalibsymp <- c(pcalibsymp, colMeans(calib.final == "Y"))
    }
    fig.file <- file.path(fig.dir, paste0(setting$country, "_L_U_sympdiff.pdf"))
    plot_df <- data.frame(ptest = ptestsymp, pcalib = pcalibsymp)
    sympdiff_plot <-
        ggplot(plot_df, aes(x = ptest, y = pcalib)) +
        geom_point(alpha = .3) +
        geom_abline(col = 'red') +
        xlab("P(Symptom_U = Yes)") +
        ylab("P(Symptom_L = Yes)") +
        ggtitle(setting$country)
    ggsave(fig.file, sympdiff_plot, width = 6, height = 6)
}


### Now plot difference between G and U
for(country in c('India', 'Tanzania')) {
    country.data <- child.clean[countries == country,-(1:2)]
    train.final <- child.clean[countries != country,-(1:2)]
    ptrainsymp <- colMeans(train.final == "Y")
    ptestsymp <- colMeans(country.data == "Y")
    fig.file <- file.path(fig.dir, paste0(country, "_G_U_sympdiff.pdf"))
    plot_df <- data.frame(ptest = ptestsymp, ptrainsymp = ptrainsymp)
    sympdiff_plot <-
        ggplot(plot_df, aes(x = ptest, y = ptrainsymp)) +
        geom_point(alpha = .3) +
        geom_abline(col = 'red') +
        xlab("P(Symptom_U = Yes)") +
        ylab("P(Symptom_G = Yes)") +
        ggtitle(country)
    ggsave(fig.file, sympdiff_plot, width = 6, height = 6)
}

