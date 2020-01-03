library(tidyverse)
library(here)
results <- readRDS(here("SimulationStudy", "data", "simulation_results.rds"))
viz.dir <- here("SimulationStudy", "figs")

if(!dir.exists(viz.dir)){
    dir.create(viz.dir, recursive = TRUE)
}

x=as.character(results$method)
x[which(x=='tariff_calib')]='Tariff_L'
x[which(x=='tariff_train_and_calib')]='Tariff_G_and_L'
x[which(x=='tariff_train')]='Tariff_G'
x[which(x=='tariff_mle')]='Tariff_NTL'
x[which(x=='tariff_revamp')]='Tariff_BTL'
x[which(x=='tariff_revamp_informed')]='Tariff_BTL_Informed'
x[which(x=='insilico_calib')]='InSilico_L'
x[which(x=='insilico_train_and_calib')]='InSilico_G_and_L'
x[which(x=='insilico_train')]='InSilico_G'
x[which(x=='insilico_mle')]='InSilico_NTL'
x[which(x=='insilico_revamp')]='InSilico_BTL'
x[which(x=='insilico_revamp_informed')]='InSilico_BTL_Informed'
x[which(x=="ensemble_revamp")]='Ensemble_J'
x[which(x=="ensemble_lite_revamp")]='Ensemble_I'
x[which(x=="ensemble_lite_revamp_informed")]='Ensemble_I_Informed'
#results$csmf.acc.split <- round(results$csmf.acc.split, 2)

results$method=as.factor(x)

x=as.character(results$fake)
x[which(x=="tariff")]="Tariff"
x[which(x=="insilico")]="InSilico"
results$fake=x

causes=c("Pneumonia", "Diarrhea/Dysentery", "Sepsis", "Other")

revamp.train.plot <-
    results %>%
    filter(measure == "csmf", method %in% paste0(fake,c("_G","_BTL"))) %>% 
    spread(method, accuracy) %>% 
    mutate(ratio=ifelse(fake=="Tariff",Tariff_BTL/Tariff_G,InSilico_BTL/InSilico_G)) %>% 
    ggplot(aes(x = csmfa.test.train, y = ratio, color = as.factor(Mtype))) +
    facet_grid(calib.size ~ fake) +
    geom_smooth() +
    # geom_point(alpha=0.2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom") +
    xlab("CSMFA between source & target domains") +
    ylab("CSMFA(CCVA calibrated by REVAMP)/CSMFA(CCVA)") +  
    guides(color=guide_legend("Misclassification matrix M"))


causes=c("Pneumonia", "Diarrhea/Dysentery", "Sepsis", "Other")
### Compare informed vs uninformed calibration

results$informed <- ifelse(grepl("Informed", results$method),
                           "Informed Shrinkage",
                           "Uninformed Shrinkage")
results$method <- gsub("_Informed", "", results$method)

level_order <- c("Tariff_BTL", "InSilico_BTL", "Ensemble_I")
fig_informed <-
    results %>%
    filter(method %in% c('Ensemble_I', 'Ensemble_I_Informed',
                         'InSilico_BTL', 'InSilico_BTL_Informed',
                         'Tariff_BTL', 'Tariff_BTL_Informed'),
           measure == "csmf", Mtype == 2, csmfa.calib == "low",
           calib.size %in% c(50, 100)) %>%
    ggplot(aes(x = method, y = accuracy, color = informed)) +
    geom_boxplot() +
    facet_grid(calib.size ~ fake) +
    theme(legend.position="bottom",
          legend.title = element_blank()) +
    xlab("Method") +
    ylab("CSMFA")
ggsave(file.path(viz.dir, "fig_informed_calib.pdf"),
       plot = fig_informed, width = 6, height = 6)

### Now get rid of informed
results <- filter(results, informed != "Informed Shrinkage")

revamp.train.plot <-
    results %>%
    filter(measure == "csmf", method %in% paste0(fake,c("_G","_BTL"))) %>% 
    spread(method, accuracy) %>% 
    mutate(ratio=ifelse(fake=="Tariff",Tariff_BTL/Tariff_G,InSilico_BTL/InSilico_G)) %>% 
    ggplot(aes(x = csmfa.test.train, y = ratio, color = as.factor(Mtype))) +
    facet_grid(calib.size ~ fake) +
    geom_smooth() +
    geom_point(alpha=0.2, size = .5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom") +
    xlab("CSMFA between source & target domains") +
    ylab("CSMFA(CCVA calibrated by REVAMP)/CSMFA(CCVA)") +  
    guides(color=guide_legend("Misclassification matrix M"))
ggsave(filename = file.path(viz.dir, "revamp_performance.pdf"),
       plot = revamp.train.plot, width = 6, height = 9)


revamp.train.400.plot <-
    results %>%
    filter(measure == "csmf", method %in% paste0(fake,c("_G","_BTL")), calib.size==400) %>% 
    spread(method, accuracy) %>% 
    mutate(ratio=ifelse(fake=="Tariff",Tariff_BTL/Tariff_G,InSilico_BTL/InSilico_G)) %>% 
    ggplot(aes(x = csmfa.test.train, y = ratio, color = as.factor(Mtype))) +
    facet_grid(~ fake) +
    geom_smooth() +
    geom_point(alpha=0.2, size = .5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom") +
    xlab("CSMFA between source & target domains") +
    #ylab("CSMFA(CCVA calibrated by REVAMP)/CSMFA(CCVA)") +  
    guides(color=guide_legend("Misclassification matrix M"))

ggsave(filename = file.path(viz.dir, "revamp400_performance.pdf"),
       plot = revamp.train.400.plot, width = 6, height = 3)

for(f in c("Tariff","InSilico")) {
    csmf.plot <-
        results %>% 
        filter(measure == "csmf", fake==f, grepl(f, method) & !grepl("NTL", method)) %>%
        #group_by(method,Mtype,calib.size,csmfa.calib) %>%
        #summarize(med.acc=median(accuracy),low.acc=quantile(accuracy,0.025),high.acc=quantile(accuracy,0.975)) %>%
        ggplot(aes(x = csmfa.calib, y = accuracy, color = method)) +
        facet_grid(calib.size ~ Mtype) +
        geom_boxplot() +
        ylab("CSMFA") +
        xlab("") +   
        #geom_errorbar(mapping=aes(ymin=low.acc, ymax=high.acc),width=0.2,position = position_dodge(0.9)) +
        theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position = "bottom")
    
    ccc.plot <-
        results %>% 
        filter(measure == "ccc", fake==f, method %in% c("Tariff_G","InSilico_G","Tariff_BTL","InSilico_BTL","Ensemble_I")) %>%
        #group_by(method,Mtype,calib.size,csmfa.calib) %>%
        #summarize(med.acc=median(accuracy),low.acc=quantile(accuracy,0.025),high.acc=quantile(accuracy,0.975)) %>%
        ggplot(aes(x = csmfa.calib, y = accuracy, color = method)) +
        facet_grid(calib.size ~ Mtype) +
        geom_boxplot() +
        ylab("CCC") +
        xlab("") +   
        #geom_errorbar(mapping=aes(ymin=low.acc, ymax=high.acc),width=0.2,position = position_dodge(0.9)) +
        theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position = "bottom")
    
    ens.csmf.plot <-
        results %>% 
        filter(measure == "csmf", fake==f, method %in% c("Tariff_G","InSilico_G","Tariff_BTL","InSilico_BTL","Ensemble_I")) %>%
        #group_by(method,Mtype,calib.size,csmfa.calib) %>%
        #summarize(med.acc=median(accuracy),low.acc=quantile(accuracy,0.025),high.acc=quantile(accuracy,0.975)) %>%
        ggplot(aes(x = csmfa.calib, y = accuracy, color = method)) +
        facet_grid(calib.size ~ Mtype) +
        geom_boxplot() +
        ylab("CSMFA") +
        xlab("") +   
        #geom_errorbar(mapping=aes(ymin=low.acc, ymax=high.acc),width=0.2,position = position_dodge(0.9)) +
        theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position = "bottom")
    
    
    ggsave(filename = file.path(viz.dir, paste0(f,"_csmf_results.pdf")),
           plot = csmf.plot, width = 7, height = 9)
    ggsave(filename = file.path(viz.dir, paste0(f,"_ccc_results.pdf")),
           plot = ccc.plot, width = 7, height = 9)
    ggsave(filename = file.path(viz.dir, paste0(f,"_ensemble_csmf_results.pdf")),
           plot = ens.csmf.plot, width = 7, height = 9)
}

#### Remake Fig S2 with mean of CSMF
ens.csmf.plot <-
    results %>% 
    filter(measure == "csmf", fake=='Tariff', method %in% c("Tariff_G","InSilico_G","Tariff_BTL","InSilico_BTL","Ensemble_I")) %>%
    mutate(Mtype = paste0("M", Mtype)) %>%
    group_by(method,Mtype,calib.size,csmfa.calib) %>%
    summarize(mean_csmfa = mean(accuracy)) %>%
    ggplot(aes(x = csmfa.calib, y = mean_csmfa, shape = method)) +
    facet_grid(calib.size ~ Mtype) +
    geom_point(position = position_dodge(width = 0.7), size = 1.5) +
    ylab("Mean CSMFA") +
    ylim(0.6, 1) +
    xlab("") +   
    #geom_errorbar(mapping=aes(ymin=low.acc, ymax=high.acc),width=0.2,position = position_dodge(0.9)) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position = "bottom")+
    scale_shape_manual(values = c('Tariff_G' = 0,
                                  'InSilico_G' = 1,
                                  'Tariff_BTL' = 15,
                                  'InSilico_BTL' = 16,
                                  'Ensemble_I' = 17))
ggsave(filename = file.path(viz.dir, "Tariff_ensemble_csmf_mean_results.pdf"),
       plot = ens.csmf.plot, width = 7, height = 9)

### Do this for just sample size 400
#### Remake Fig S2 with mean of CSMF
for(f in c("Tariff", "InSilico")) {
    ens.csmf.plot <-
        results %>% 
        filter(measure == "csmf", fake==f, method %in% c("Tariff_G","InSilico_G","Tariff_BTL","InSilico_BTL","Ensemble_I"), calib.size == 400) %>%
        mutate(Mtype = paste0("M", Mtype)) %>%
        group_by(method,Mtype,calib.size,csmfa.calib) %>%
        summarize(mean_csmfa = mean(accuracy)) %>%
        ggplot(aes(x = csmfa.calib, y = mean_csmfa, shape = method)) +
        facet_wrap( ~ Mtype) +
        geom_point(position = position_dodge(width = 0.7), size = 2.0) +
        ylab("Mean CSMFA") +
        ylim(0.6, 1) +
        xlab("") +   
        #geom_errorbar(mapping=aes(ymin=low.acc, ymax=high.acc),width=0.2,position = position_dodge(0.9)) +
        theme_bw()
        theme(axis.text.x = element_text(angle = 0, hjust = 1),
              legend.position = "bottom")+
        scale_shape_manual(values = c('Tariff_G' = 0,
                                      'InSilico_G' = 1,
                                      'Tariff_BTL' = 15,
                                      'InSilico_BTL' = 16,
                                      'Ensemble_I' = 17))
    ggsave(filename = file.path(viz.dir, paste0(f, "_ensemble_csmf_mean_results_400.pdf")),
           plot = ens.csmf.plot, width = 9, height = 4)
}

##### cause specific barplots #####
cause.plot <-
    results %>%
    filter(measure %in% causes, str_detect(as.character(method),as.character(fake)), 
           !grepl("_L",method), !grepl("_NTL",method), calib.size==400) %>%
    group_by(method,measure,Mtype,fake) %>%
    summarize(med.acc=median(accuracy),low.acc=quantile(accuracy,0.025),high.acc=quantile(accuracy,0.975)) %>%
    ggplot(aes(x = measure, y = med.acc, fill = method)) +
    facet_grid(fake ~ Mtype) +
    geom_bar(stat="identity",position=position_dodge()) +
    ylab("bias") +
    #geom_errorbar(mapping=aes(ymin=low.acc, ymax=high.acc),width=0.2,position = position_dodge(0.9)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "bottom")

ggsave(filename = file.path(viz.dir, "cause_biases.pdf"),
       plot = cause.plot, width = 6, height = 4)


csmf.mle.plot <-
    results %>%
    filter(measure=="csmf", str_detect(as.character(method),as.character(fake)), 
           grepl("_BTL",method) | grepl("_NTL",method)) %>%
    ggplot(aes(x = as.factor(calib.size), y = accuracy, color = method)) +
    facet_grid(fake ~ Mtype) +
    geom_boxplot() +
    ylab("CSMFA") +
    xlab("n") + 
    #geom_errorbar(mapping=aes(ymin=low.acc, ymax=high.acc),width=0.2,position = position_dodge(0.9)) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position = "bottom")

ggsave(filename = file.path(viz.dir, "mle_vs_revamp.pdf"),
       plot = csmf.mle.plot, width = 7, height = 4)

##### Ensemble

revamp.ensemble.plot <-
    results %>% 
    filter(measure == "csmf", grepl("_BTL",method) | grepl("Ensemble",method)) %>% 
    spread(method, accuracy) %>% 
    # mutate(lambda=ifelse(fake=="tariff",(ensemble_lite_revamp-insilico_revamp)/(tariff_revamp-insilico_revamp),
    #   (ensemble_revamp-tariff_revamp)/(insilico_revamp-tariff_revamp)),
    #   true_revamp=ifelse(fake=="tariff",tariff_revamp,insilico_revamp),
    #   false_revamp=ifelse(fake=="insilico",tariff_revamp,insilico_revamp)) %>% 
    transmute(Mtype,fake,calib.size,max_revamp=pmax(Tariff_BTL,InSilico_BTL),
              min_revamp=pmin(Tariff_BTL,InSilico_BTL),ind=Ensemble_I,
              joint=Ensemble_J) %>%
    gather(value=accuracy,key=ensemble_model,-Mtype,-fake,-calib.size,-min_revamp,-max_revamp)  %>%
    ggplot(aes(x = max_revamp - min_revamp ,
               y = accuracy - min_revamp,
               color = as.factor(Mtype), linetype=as.factor(ensemble_model))) +
    facet_grid(fake ~ calib.size) +
    geom_smooth() +
    #geom_point(alpha=0.2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom") +
    #xlab("CSMFA(REVAMP with Best CCVA) - CSMFA(REVAMP with Worst CCVA)") +
    #ylab("CSMFA(Ensemble REVAMP) - CSMFA(REVAMP with Worst CCVA)") +  
    xlab(expression(delta)) + 
    ylab(expression(nu)) + 
    guides(color=guide_legend("Misclassification matrix M"),
           linetype=guide_legend("Ensemble model",override.aes = list(col = 'black')))

ggsave(filename = file.path(viz.dir, "ensemble_performance.pdf"),
       plot = revamp.ensemble.plot, width = 8, height = 4)
##########################################
results$csfma.calib <- factor(results$csmfa.calib, levels = c("low", "medium", "high"))

### Fig. 2
fig2 <-
    results %>%
    filter(method %in% c('Ensemble_I',
                         'InSilico_BTL', 'InSilico_G',
                         'Tariff_BTL', 'Tariff_G'),
           measure == "csmf", fake == "tariff") %>%
    ggplot(aes(x = csfma.calib, y = accuracy, color = method)) +
    geom_boxplot() +
    facet_grid(calib.size ~ Mtype) +
    theme(legend.position="bottom")
ggsave(file.path(viz.dir, "fig2.pdf"),
       plot = fig2, width = 14, height = 12)

### Fig S1
s1_df <-
    results %>%
    filter(measure == "csmf", calib.size == 400) %>%
    group_by(run, Mtype, fake, csmfa.calib) %>%
    summarize(csmfa.test.train = mean(csmfa.test.train),
              insilico_ratio = accuracy[method == "InSilico_BTL"] /accuracy[method == "InSilico_G"],
              tariff_ratio = accuracy[method == "Tariff_BTL"] /accuracy[method == "Tariff_G"]) %>%
    mutate(ratio = ifelse(fake == 'tariff', tariff_ratio , insilico_ratio))
figS1 <-
    ggplot(s1_df, aes(x = csmfa.test.train, y = ratio, color = factor(Mtype))) +
    facet_wrap(~fake) +
    geom_smooth()
ggsave(file.path(viz.dir, "figS1.pdf"),
       plot = figS1, width = 14, height = 12)

### Fig S2
figS2 <-
    results %>%
    filter(method %in% c('InSilico_BTL', 'InSilico_G', 'Tariff_BTL', 'Tariff_G'),
           !(measure %in% c('ccc', 'csmf'))) %>%
    group_by(method, measure, method, Mtype, fake) %>%
    summarize(bias = mean(accuracy)) %>%
    ggplot(aes(x = measure, y = bias, fill = method)) +
    geom_bar(stat="identity") +
    facet_grid(fake ~ Mtype)

### Fig S3
figS3 <-
    results %>%
    filter(method %in% c('InSilico_BTL', 'InSilico_G',
                         'InSilico_G_and_L', 'InSilico_L'),
           measure == "csmf", fake == "insilico") %>%
    ggplot(aes(x = csfma.calib, y = accuracy, color = method)) +
    geom_boxplot() +
    facet_grid(calib.size ~ Mtype) +
    theme(legend.position="bottom")

### Fig S6
figS6 <-
    results %>%
    filter(method %in% c('Ensemble_I',
                         'InSilico_BTL', 'InSilico_G',
                         'Tariff_BTL', 'Tariff_G'),
           measure == 'ccc', fake == 'insilico' )%>%
    ggplot(aes(x = csfma.calib, y = accuracy, color = method)) +
    geom_boxplot() +
    facet_grid(calib.size ~ Mtype) +
    theme(legend.position="bottom")



    




for(m in 1:3) {
    for(f in c("tariff","insilico")) {
        csmf.plot <-
            results %>% 
            filter(measure == "csmf", Mtype==m, fake==f) %>%
            ggplot(aes(x = method, y = accuracy, color = is.revamp)) +
            facet_grid(calib.size ~ csmfa.calib) +
            geom_boxplot() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        ccc.plot <-
            results %>%
            filter(measure == "ccc", Mtype==m, fake==f) %>%
            ggplot(aes(x = method, y = accuracy, color = is.revamp)) +
            facet_grid(calib.size ~ csmfa.calib) +
            geom_boxplot() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        cause.plot <-
            results %>%
            filter(measure %in% causes, Mtype==m, fake==f) %>%
            ggplot(aes(x = method, y = accuracy, color = is.revamp)) +
            facet_grid(calib.size ~ measure) +
            geom_boxplot() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        
        ggsave(filename = file.path(viz.dir, paste0("M",m,"_",f,"_csmf_results.pdf")),
               plot = csmf.plot, width = 14, height = 12)
        ggsave(filename = file.path(viz.dir, paste0("M",m,"_",f,"_ccc_results.pdf")),
               plot = ccc.plot, width = 14, height = 12)
        ggsave(filename = file.path(viz.dir, paste0("M",m,"_",f,"_cause_results.pdf")),
               plot = cause.plot, width = 14, height = 12)
    }
}


### if we just need to plot the cause boxplots 
for(m in 1:3) for(f in c("tariff","insilico")) {
     cause.plot <-
             results %>%
             filter(measure %in% causes, Mtype==m, fake==f) %>%
             ggplot(aes(x = method, y = accuracy, color = is.revamp)) +
             facet_grid(calib.size ~ measure) +
             geom_boxplot() + 
             theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(filename = file.path(viz.dir, paste0("M",m,"_",f,"_cause_results.pdf")),
                        plot = cause.plot, width = 14, height = 12)
    }


### median numbers ##
for(f in c("tariff","insilico")) {
    cause.barplot <- results %>%
    filter(measure %in% causes, fake==f, grepl(f, method) | grepl("revamp", method)) %>% 
    group_by(Mtype,measure,method,calib.size) %>% 
    summarise(low=quantile(accuracy,0.25),high=quantile(accuracy,0.75),accuracy=median(accuracy)) %>% data.frame() %>%
        mutate(is.revamp=grepl("revamp", method)) %>% 
        ggplot(aes(x = method, y = accuracy, fill=as.factor(Mtype))) +
        facet_grid(calib.size ~ measure) +
        geom_bar(stat="identity",position="dodge") + geom_hline(yintercept = 0) +
        #geom_errorbar(aes(ymin=low,ymax=high)) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(filename = file.path(viz.dir, paste0(f,"_cause_barplots.pdf")),
           plot = cause.barplot, width = 14, height = 12)
        
}

for(f in c("tariff","insilico")) {
    cause.barplot <- results %>%
        filter(measure=="csmf", fake==f, grepl(f, method) | grepl("revamp", method)) %>% 
        group_by(Mtype,measure,method,calib.size,csmfa.calib) %>% 
        summarise(low=quantile(accuracy,0.25),high=quantile(accuracy,0.75),accuracy=median(accuracy)) %>% data.frame() %>%
        mutate(is.revamp=grepl("revamp", method)) %>% 
        ggplot(aes(x = method, y = accuracy, fill=as.factor(Mtype))) +
        facet_grid(calib.size ~ csmfa.calib) + 
        geom_bar(stat="identity",position="dodge") + geom_hline(yintercept = 0) +
        #geom_errorbar(aes(ymin=low,ymax=high)) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(filename = file.path(viz.dir, paste0(f,"_csmf_barplots.pdf")),
           plot = cause.barplot, width = 14, height = 12)
    
}

f="tariff"
x=results %>% filter(csmfa.calib=="high",calib.size==50,Mtype==1,fake==f,
    method %in% paste0(f,"_",c("train","revamp"))) 


x = x %>% mutate(id = group_indices_(x,.dots=c("measure","run"))) %>% select(accuracy,method,measure,run,id)
y= x %>% spread(method,accuracy)
    
    group_by(measure) %>% summarize() %>% data.frame()

library(tidyr)
x %>% spread(method,accuracy)
x
group_by(measure) %>% summarize(median(accuracy)) %>% data.frame()
