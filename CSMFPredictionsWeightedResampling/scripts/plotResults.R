library(tidyverse)
library(here)

fig.dir <- here("CSMFPredictionsWeightedResampling", "figs")
if(!dir.exists(fig.dir)) {
    dir.create(fig.dir, recursive = TRUE)
}

### CSMF ccuracy boxplots
### First child data
child_results <- readRDS(here("CSMFPredictionsWeightedResampling", "child_data", "csmfa_results.rds"))
child_results$method <- factor(child_results$method,
                         levels = c('tariff_train', 'insilico_train',
                                    'tariff_calibva', 'insilico_calibva',
                                    'ensemble_calibva', 'tariff_train_and_calib',
                                    'insilico_train_and_calib', 'calib'),
                         labels = c('Tariff_G', 'InSilicoVA_G', 'Tariff_BTL', 'InSilicoVA_BTL',
                                    'Ensemble_I', 'Tariff_G_and_L', 'InSilicoVA_G_and_L', 'L_avg'))
################# Main figure
plot_df <-
    child_results %>%
    filter(cod == 5, alpha == 5) %>%
    filter(grepl("(*_G$|*_BTL$|*_I$)", method)) %>%
    filter(calib.csfma == "low") %>%
    group_by(country, method, calib.size) %>%
    summarize(mean_csmfa = mean(csmfa)) %>%
    ungroup()

main_fig <- 
    plot_df %>%
    ggplot(aes(x = method, y = mean_csmfa)) +
    geom_point(aes(shape = method), size = 2.0) +
    facet_grid(country ~ calib.size)+
    theme_bw() +
    theme(text = element_text(size=14),
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") +
    xlab("Method") +
    ylab("Mean CSMFA") +
    scale_y_continuous(limits=c(.6, 1), breaks = seq(.6, 1, by = .1)) +
    scale_shape_manual(values = c('Tariff_G' = 0,
                                  'InSilicoVA_G' = 1,
                                  'Tariff_BTL' = 15,
                                  'InSilicoVA_BTL' = 16,
                                  'Ensemble_I' = 17)) 
ggsave(file.path(fig.dir, "main_csmfa_fig.jpg"),
       plot = main_fig, width = 6, height = 4)

#### Figure comparing 7 COD to 5
plot_df <-
    child_results %>%
    filter(alpha == 5) %>%
    filter(grepl("(*_G$|*_BTL$|*_I$)", method)) %>%
    filter(calib.csfma == "low") %>%
    group_by(country, method, calib.size, cod) %>%
    summarize(mean_csmfa = mean(csmfa)) %>%
    ungroup()
cod_7_fig <- 
    plot_df %>%
    ggplot(aes(x = method, y = mean_csmfa, color = factor(cod))) +
    #geom_jitter(aes(shape = method), width = .5) +
    geom_point(aes(shape = method), position = position_dodge(width=0.5)) +
    facet_grid(country ~ calib.size)+
    theme(text = element_text(size=14),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("Method") +
    ylab("Mean CSMFA") +
    labs(color = "Number COD") +
    scale_y_continuous(limits=c(.6, 1), breaks = seq(.6, 1, by = .1))+
    scale_shape_manual(values = c('Tariff_G' = 0,
                                  'InSilicoVA_G' = 1,
                                  'Tariff_BTL' = 15,
                                  'InSilicoVA_BTL' = 16,
                                  'Ensemble_I' = 17)) +
    guides(shape = FALSE)
ggsave(file.path(fig.dir, "cod_5_vs_7_csmfa_fig.jpg"),
       plot = cod_7_fig, width = 6, height = 4)

######################## 
### Add calibration average to figure
plot_df <-
    child_results %>%
    filter(cod == 5, alpha == 5) %>%
    filter(grepl("(*_G$|*_BTL$|*_I$|L_avg)", method)) %>%
    filter(calib.csfma == "low") %>%
    group_by(country, method, calib.size) %>%
    summarize(mean_csmfa = mean(csmfa)) %>%
    ungroup()

main_fig_with_calib <- 
    plot_df %>%
    ggplot(aes(x = method, y = mean_csmfa)) +
    geom_point(aes(shape = method)) +
    facet_grid(country ~ calib.size)+
    theme(text = element_text(size=14),
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") +
    xlab("Method") +
    ylab("Mean CSMFA") +
    scale_y_continuous(limits=c(.3, 1), breaks = seq(.3, 1, by = .1)) +
    scale_shape_manual(values = c('Tariff_G' = 0,
                                  'InSilicoVA_G' = 1,
                                  'Tariff_BTL' = 15,
                                  'InSilicoVA_BTL' = 16,
                                  'Ensemble_I' = 17,
                                  'L_avg' = 5)) 
ggsave(file.path(fig.dir, "main_with_calib_csmfa_fig.jpg"),
       plot = main_fig_with_calib, width = 6, height = 4)

##########################
#### Plot showing CSMFA compared to using calibration set in training
plot_df <-
    child_results %>%
    filter(cod == 5, alpha == 5) %>%
    filter(method != "L_avg") %>%
    filter(calib.csfma == "low") %>%
    group_by(country, method, calib.size) %>%
    summarize(mean_csmfa = mean(csmfa)) %>%
    ungroup()

main_fig_with_calib_train <- 
    plot_df %>%
    ggplot(aes(x = method, y = mean_csmfa)) +
    geom_point(aes(shape = method)) +
    facet_grid(country ~ calib.size)+
    theme(text = element_text(size=14),
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") +
    xlab("Method") +
    ylab("Mean CSMFA") +
    scale_y_continuous(limits=c(.6, 1), breaks = seq(.6, 1, by = .1)) +
    scale_shape_manual(values = c('Tariff_G' = 0,
                                  'InSilicoVA_G' = 1,
                                  'Tariff_BTL' = 15,
                                  'InSilicoVA_BTL' = 16,
                                  'Ensemble_I' = 17,
                                  'Tariff_G_and_L' = 12,
                                  'InSilicoVA_G_and_L' = 10)) 
ggsave(file.path(fig.dir, "main_with_calib_train_csmfa_fig.jpg"),
       plot = main_fig_with_calib_train, width = 6, height = 4)

###########################
setting.df <- as.data.frame(expand.grid(calib.csfma = c('low', 'representative'),
                          stringsAsFactors = FALSE))
pdf(file.path(fig.dir, "child_csmfa_5_cod_alpha_5_mean.pdf"), width = 6, height = 6)
for(i in 1:nrow(setting.df)) {
    setting <- setting.df[i,]
    title <- paste(setting, "CSMFA")
    plot_df <-
        child_results %>%
        filter(cod == 5, alpha == 5) %>%
        filter(grepl("(*_G$|*_BTL$|*_I$)", method)) %>%
        filter(calib.csfma == setting) %>%
        group_by(country, method, calib.size) %>%
        summarize(mean_csmfa = mean(csmfa)) %>%
        ungroup()
    plot <- 
        plot_df %>%
        ggplot(aes(x = method, y = mean_csmfa)) +
        geom_point() +
        facet_grid(country ~ calib.size)+
        theme(text = element_text(size=14),
              axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(title) +
        ylab("Mean CSMFA") +
        scale_y_continuous(limits=c(.6, 1), breaks = seq(.6, 1, by = .1))
    print(plot)
}
dev.off()

pdf(file.path(fig.dir, "child_csmfa_5_cod_alpha_5_boxplot.pdf"), width = 6, height = 6)
for(i in 1:nrow(setting.df)) {
    setting <- setting.df[i,]
    title <- paste(setting, "CSMFA")
    
    plot <-
        child_results %>%
        filter(cod == 5, alpha == 5) %>%
        filter(grepl("(*_G$|*_BTL$|*_I$)", method)) %>%
        filter(calib.csfma == setting) %>%
        ggplot(aes(x = method, y = csmfa)) +
        geom_boxplot(outlier.size = .25) +
        facet_grid(country ~ calib.size)+
        theme(text = element_text(size=14),
              axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(title) +
        ylim(0.4, 1)
    print(plot)
}
dev.off()





############## Look at India, low calib CSMFA for actual CSMF estimates
csmf.results <- readRDS(here("CSMFPredictionsWeightedResampling", "child_data", "csmf_results.rds"))
csmf.results$method <- factor(csmf.results$method,
                                         levels = c('tariff_train', 'insilico_train',
                                                    'tariff_calibva', 'insilico_calibva',
                                                    'ensemble_calibva', 'tariff_train_and_calib',
                                                    'insilico_train_and_calib', 'calib'),
                                         labels = c('Tariff_G', 'Insilico_G', 'Tariff_BTL', 'Insilico_BTL',
                                                    'Ensemble_I', 'Tariff_G_and_L', 'Insilico_G_and_L', 'H_avg'))

india_csmf_plot_low_broad <-
    csmf.results %>%
    filter(!(method %in% c('Tariff_G_and_L', 'Insilico_G_and_L'))) %>%
    filter(cod == "broad", calib.csfma == "low", country == "India") %>%
    ggplot(aes(x = method, y = csmf.est)) +
    geom_boxplot() +
    facet_grid(calib.size ~ cause) +
    geom_hline(aes(yintercept = true.csmf), col = 'red') +
    theme(text = element_text(size=14),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(0, .6)
ggsave(file.path(fig.dir, "India_CSMF_Est_Low_Broad.pdf"),
       india_csmf_plot_low_broad,
       width = 10, height = 10)

india_csmf_plot_rep_broad <-
    csmf.results %>%
    filter(!(method %in% c('Tariff_G_and_L', 'Insilico_G_and_L'))) %>%
    filter(cod == "broad", calib.csfma == "representative", country == "India") %>%
    ggplot(aes(x = method, y = csmf.est)) +
    geom_boxplot() +
    facet_grid(calib.size ~ cause) +
    geom_hline(aes(yintercept = true.csmf), col = 'red') +
    theme(text = element_text(size=14),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(0, .6)
ggsave(file.path(fig.dir, "India_CSMF_Est_Rep_Broad.pdf"),
       india_csmf_plot_rep_broad,
       width = 10, height = 10)

Tanzania_csmf_plot_low_broad <-
    csmf.results %>%
    filter(!(method %in% c('Tariff_G_and_L', 'Insilico_G_and_L'))) %>%
    filter(cod == "broad", calib.csfma == "low", country == "Tanzania") %>%
    ggplot(aes(x = method, y = csmf.est)) +
    geom_boxplot() +
    facet_grid(calib.size ~ cause) +
    geom_hline(aes(yintercept = true.csmf), col = 'red') +
    theme(text = element_text(size=14),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(0, .6)
ggsave(file.path(fig.dir, "Tanzania_CSMF_Est_Low_Broad.jpg"),
       Tanzania_csmf_plot_low_broad,
       width = 10, height = 10)

Tanzania_csmf_plot_rep_broad <-
    csmf.results %>%
    filter(!(method %in% c('Tariff_G_and_L', 'Insilico_G_and_L'))) %>%
    filter(cod == "broad", calib.csfma == "representative", country == "Tanzania") %>%
    ggplot(aes(x = method, y = csmf.est)) +
    geom_boxplot() +
    facet_grid(calib.size ~ cause) +
    geom_hline(aes(yintercept = true.csmf), col = 'red') +
    theme(text = element_text(size=14),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(0, .6)
ggsave(file.path(fig.dir, "Tanzania_CSMF_Est_Rep_Broad.pdf"),
       Tanzania_csmf_plot_rep_broad,
       width = 10, height = 10)


############# Compare expanded causes vs true causes
expanded_csmf <- readRDS(here("CSMFPredictionsWeightedResampling", "child_data", "csmfa_expanded_causes_results.rds"))

expanded_csmf <-
    expanded_csmf %>%
    filter(grepl("calibva", method)) %>%
    mutate(calibrated_causes = '7')

original_csmf <- readRDS(here("CSMFPredictionsWeightedResampling", "child_data", "csmfa_results.rds"))
original_csmf <- 
    original_csmf %>%
    filter(cod == "broad", calib.csfma == "representative", grepl("calibva", method)) %>%
    select(method, csmfa, run, country, calib.size)%>%
    mutate(calibrated_causes = '5')

plot_df <-
    bind_rows(expanded_csmf, original_csmf)
plot_df$method <- factor(plot_df$method,
                               levels = c('tariff_calibva', 'insilico_calibva',
                                          'ensemble_calibva'),
                               labels = c('Tariff_BTL', 'Insilico_BTL',
                                          'Ensemble_I'))

expanded_csmfa_plot <-
    ggplot(plot_df, aes(x = method, y = csmfa, color = calibrated_causes)) +
    geom_boxplot() +
    facet_grid(country ~ calib.size) +
    guides(color=guide_legend(title="Number of causes")) +
    theme(text = element_text(size=14),
          axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(fig.dir, "expanded_causes_csmfa.jpg"), expanded_csmfa_plot,
       width = 6, height = 4)

############## Now adult


adult_results <- readRDS(here("CSMFPredictionsWeightedResampling", "adult_data", "csmfa_results.rds"))
adult_results$method <- factor(adult_results$method,
                               levels = c('tariff_train', 'insilico_train', 'nbc_train',
                                          'tariff_calibva', 'insilico_calibva', 'nbc_calibva',
                                          'ensemble_calibva', 'calib'),
                               labels = c('Tariff_G', 'Insilico_G', "NBC_G",
                                          'Tariff_BTL', 'Insilico_BTL', "NBC_BTL",
                                          'Ensemble_I', 'H_avg'))



plot_df <- expand.grid(calib.csfma = c("low", "representative"))
for(i in 1:nrow(plot_df)){
    title <- paste0("adult_csmfa_", plot_df$calib.csfma[i], "calibcsfma.jpg")
    csmfa.plot <-
        adult_results %>%
        filter(calib.csfma == plot_df$calib.csfma[i]) %>%
        ggplot(aes(x = method, y = csmfa)) +
        facet_grid(country ~ calib.size) +
        geom_boxplot() +
        ylab("CSMFA") +
        theme(text = element_text(size=14),
              axis.text.x = element_text(angle = 90, hjust = 1))+
        ggtitle(title)  +
        ylim(0,1)
    file <- file.path(fig.dir,
                      paste0("adult_csmfa_", plot_df$calib.csfma[i], "calibcsfma.jpg"))
    ggsave(filename = file, plot = csmfa.plot, width = 6, height = 4)
}



