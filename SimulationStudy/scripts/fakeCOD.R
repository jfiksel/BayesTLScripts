library(here)
setting.df <- expand.grid(run = 1:100, Mtype=1:3, calib.size = c(50, 100, 200, 400), 
                          fake=c("tariff","insilico"), csmfa.calib=c("low","medium","high"))
### i is row of setting.df
#i <- as.numeric(commandArgs(trailingOnly = TRUE))
i <- as.numeric(commandArgs(trailingOnly = TRUE))
#i <- as.numeric(Sys.getenv('SGE_TASK_ID'))
#i=1601
#i = 1207 ### for cluster runs set i=as.numeric(Sys.getenv('SGE_TASK_ID')) or i=as.numeric(commandArgs(trailingOnly = TRUE))
setting <- setting.df[i,]
### Make data directory for this run
sim.dir <- here("SimulationStudy", "cluster_output", paste0("M",setting$Mtype),
                setting$fake,
                paste0("size_", setting$calib.size),
                paste0("calib_",setting$csmfa.calib),
                paste0("run_", setting$run))
if(!dir.exists(sim.dir)){
    dir.create(sim.dir, recursive = TRUE) 
}

## Quit if results file already exists
if(file.exists(file.path(sim.dir, "results.rds"))){
    quit('no')
}

source(here("SimulationStudy", "scripts", "calibration.R"))
#library(rJava)
library(openVA)
library(dplyr)
library(Rsolnp)
library(MCMCpack)


set.seed(123)
### Read in child data
child.raw <- read.csv(getPHMRC_url("child"))
### Get data frame matching va34 code to actual COD
cause.df <- unique(child.raw[,c("gs_text34", "va34")])
### Reduce these to the top 3 COD + other (but not including Other defined causes in top 3)
### VA34 code for other defined COD is 14
top3cod <- names(sort(table(child.raw$gs_text34[child.raw$va34 != 14]),
                      decreasing = TRUE)[1:3])
top3cause.df <- cause.df[cause.df$gs_text34 %in% top3cod,]
### Add row for other
top3cause.df$gs_text34 <- as.character(top3cause.df$gs_text34)
top3cause.df <- rbind(top3cause.df, c("Other", 99))
### Clean data into output usable by Tariff & InsilicoVA
child.clean <- ConvertData.phmrc(child.raw, phmrc.type = "child")$output

### Function to change all guesses not equal to top 3 causes to other
changeTopCOD <- function(topcod, topcause.df = top3cause.df) {
  topcod <- as.character(topcod)
  ncauses <- sum(topcause.df$va34 != "99")
  topcod <- ifelse(topcod %in% topcause.df$va34[1:ncauses], topcod, "99")
  return(topcod)
}

### shifting to only top 3 causes (due to use of fakeCause for these runs)
child.clean$Cause=changeTopCOD(child.clean$Cause)

### Create list of seeds
all.seeds <- sample(1e6, size = nrow(setting.df), replace = F)
set.seed(all.seeds[i])

### defining csmfa ranges (0 - 0.4 =low, 0.4 - 0.6 = medium, > 0.6 = high)
csmfa.range=matrix(c(0,0.4,0.6,0.4,0.6,1),3,2)
colnames(csmfa.range)=c("down","up")
row.names(csmfa.range)=c("low","medium","high")

## geenrates prob vectors with minimum csmfs > 0.05 and < 0.5
rdirichlet_mod=function(v) {
  p=rdirichlet(1,v)
  while(min(p) < 0.05 | max(p) > 0.5) p=rdirichlet(1,v)
  return(p)
}


flag=0
while(flag==0){
  ptest=as.vector(rdirichlet_mod(rep(1,4)))
  print("Generating CSMFs with desired CSMFA")
  pcalib.mat=sapply(1:100000,function(i) rdirichlet_mod(rep(1,4)))
  csmfa.train=apply(pcalib.mat,2,getCSMF_accuracy,ptest)
  ind.calib=which(csmfa.train > csmfa.range[setting$csmfa.calib,"down"] & 
                    csmfa.train < csmfa.range[setting$csmfa.calib,"up"])
  if(length(ind.calib) > 0) flag=1
 }

pcalib=pcalib.mat[,sample(ind.calib,1)]

top3cause.df$ptest <- ptest
top3cause.df$pcalib <- pcalib

### error Matrix
M1=diag(4) ## perfect alignment
M2=rbind(c(1,0,0,0),c(0.65,0.35,0,0),c(0,0,0.5,0.5),c(0,0,0,1)) ## two large misclassification
M3=0.6*diag(4)+0.1 ## many small misclassifications
M=get(paste0("M",setting$Mtype))
colnames(M)=row.names(M)=top3cause.df$va34

## generating q vectors according to our model
qcalib=as.vector(t(M)%*%pcalib)
qtest=as.vector(t(M)%*%ptest)

### initial resampling of the PHMRC data into train, calib and test
### making sure calib has atleast 400 samples
ind=as.vector(t(rmultinom(nrow(child.clean),1,c(800,nrow(child.clean)-1600,800)))%*%(1:3))
while(table(ind)[2] < 400) {
  print("resampling ot ensure calibration set has atleast 400 samples")
  ind=as.vector(t(rmultinom(nrow(child.clean),1,c(800,nrow(child.clean)-1600,800)))%*%(1:3))
}

train.init <- child.clean[which(ind==1),]
calib.init <- child.clean[which(ind==2),]
test.init <- child.clean[which(ind==3),]

### training a VA algorithm on train to predict VA cause for test + calib 
## the VA cause will be used to predict fake cause for test + calib using the matrix M
if(setting$fake=="tariff"){
    set.seed(123)
    tariff.init.train <- codeVA(data = rbind(calib.init, test.init),
                           data.type = "customize", model = "Tariff",
                           data.train = train.init, causes.train = "Cause")
    init.train.cod <- changeTopCOD(getTopCOD(tariff.init.train)[,2])
}else{
  insilico.nsim <- 5000
  set.seed(123)
  insilico.init.train <- codeVA(data = rbind(calib.init, test.init),
                           data.type = "customize", model = "InSilicoVA",
                           data.train = train.init, causes.train = "Cause",
                           jump.scale = 0.05, Nsim=insilico.nsim, auto.length = FALSE)
  init.train.cod <- changeTopCOD(getTopCOD(insilico.init.train)[,2])
}

#### VA predicted causes for test + calib
train.init$VAcause = rep(NA,nrow(train.init))
calib.init$VAcause = init.train.cod[1:nrow(calib.init)]
test.init$VAcause = init.train.cod[-(1:nrow(calib.init))]

### Split test data frames by VAcause
calib.split <- split(calib.init, calib.init$VAcause)

# calib VAcause distribution (calib VA ccause distribution should be qcalib)
top3cause.df$calib.n <- floor(qcalib * setting$calib.size)
top3cause.df$calib.n[4] <- setting$calib.size - sum(top3cause.df$calib.n[1:3])

### resampling for calibration set
calib.list <- lapply(seq_along(calib.split), function(i) {
  calib.cause <- calib.split[[i]]
  cause <- names(calib.split)[i]
  ncause <- top3cause.df[top3cause.df$va34 == cause, "calib.n"]
  return(calib.cause[sample(nrow(calib.cause), ncause, replace = TRUE),]) 
})
calib.final <- do.call(rbind, calib.list)

### Split test data frames by VA cause (test VA ccause distribution should be qtest)
test.split <- split(test.init, test.init$VAcause)

# Test cause distribution
top3cause.df$test.n <- floor(qtest * (nrow(test.init)))
top3cause.df$test.n[4] <- nrow(test.init) - sum(top3cause.df$test.n[1:3])

### resampling for test set
test.list <- lapply(seq_along(test.split), function(i) {
  test.cause <- test.split[[i]]
  cause <- names(test.split)[i]
  ncause <- top3cause.df[top3cause.df$va34 == cause, "test.n"]
  return(test.cause[sample(nrow(test.cause), ncause, replace = TRUE),]) 
})
test.final <- do.call(rbind, test.list)

### generating fake causes
train.init$fakeCause = train.init$Cause
train.final=train.init

## since marginal distribution of VA cause in calib is qcalib = t(M) %*% pcalib, using this formula ensures 
## marginal distribution of fake cause in calib will be pcalib
## similarly marginal distribution of fake cause in test will be ptest
calib.final$fakeCause = sapply(calib.final$VAcause,
  function(i) top3cause.df$va34[sum(rmultinom(1,1,M[,i]*top3cause.df$pcalib)*(1:nrow(top3cause.df)))])
test.final$fakeCause = sapply(test.final$VAcause,
  function(i) top3cause.df$va34[sum(rmultinom(1,1,M[,i]*top3cause.df$ptest)*(1:nrow(top3cause.df)))])

### Get causes (these will be the ordering we will use)
causes <- top3cause.df$va34

## removing extra cause columns 
train.final = train.final %>% dplyr::select(- c(Cause,VAcause))
test.final = test.final %>% dplyr::select(- c(Cause,VAcause))
calib.final = calib.final %>% dplyr::select(- c(Cause,VAcause))

top3cause.df$ptrain.final <- table(changeTopCOD(train.final$fakeCause))[top3cause.df$va34]/nrow(train.final)
top3cause.df$ptest.final <- table(changeTopCOD(test.final$fakeCause))[top3cause.df$va34]/nrow(test.final)
top3cause.df$pcalib.final <- table(changeTopCOD(calib.final$fakeCause))[top3cause.df$va34]/nrow(calib.final)



################################
### Now the 6 methods w/o calibration
### Tariff with just training set (will also predict on calibration data)

set.seed(123)
tariff.train <- codeVA(data = rbind(calib.final, test.final),
                       data.type = "customize", model = "Tariff",
                       data.train = train.final, causes.train = "fakeCause")
tariff.train.cod <- changeTopCOD(getTopCOD(tariff.train)[,2])

### Get the COD estimates for test & calibration set from tariff & insilicova
tariff.train.cod.test <- tariff.train.cod[-(1:setting$calib.size)]
tariff.train.cod.calib <- tariff.train.cod[1:setting$calib.size]

### Tariff with training + calibration
set.seed(123)
tariff.train.calib <- codeVA(data = test.final,
                             data.type = "customize", model = "Tariff",
                             data.train = rbind(train.final, calib.final), causes.train = "fakeCause")
tariff.train.calib.cod <- changeTopCOD(getTopCOD(tariff.train.calib)[,2])

### Tariff with just calibration
set.seed(123)
tariff.calib <- codeVA(data = test.final,
                       data.type = "customize", model = "Tariff",
                       data.train = calib.final, causes.train = "fakeCause")
tariff.calib.cod <- changeTopCOD(getTopCOD(tariff.calib)[,2])


### Set number of simulations for insilico
insilico.nsim <- 5000
### Commented out line for interactive testing
#insilico.nsim <- 500

### Insilico with just training set
set.seed(123)
insilico.train <- codeVA(data = rbind(calib.final, test.final),
                         data.type = "customize", model = "InSilicoVA",
                         data.train = train.final, causes.train = "fakeCause",
                         jump.scale = 0.05, Nsim=insilico.nsim, auto.length = FALSE)
insilico.train.cod <- changeTopCOD(getTopCOD(insilico.train)[,2])

### Insilico with training + calibration set
set.seed(123)
insilico.train.calib <- codeVA(data = test.final,
                               data.type = "customize", model = "InSilicoVA",
                               data.train = rbind(train.final, calib.final), causes.train = "fakeCause",
                               jump.scale = 0.05, Nsim=insilico.nsim, auto.length = FALSE)
insilico.train.calib.cod <- changeTopCOD(getTopCOD(insilico.train.calib)[,2])


### Insilico with just calibration set
set.seed(123)
insilico.calib <- codeVA(data = test.final,
                         data.type = "customize", model = "InSilicoVA",
                         data.train = calib.final, causes.train = "fakeCause",
                         jump.scale = 0.05, Nsim=insilico.nsim, auto.length = FALSE) 
insilico.calib.cod <- changeTopCOD(getTopCOD(insilico.calib)[,2])



### VA predicted causes
insilico.train.cod.test <- insilico.train.cod[-(1:setting$calib.size)]
insilico.train.cod.calib <- insilico.train.cod[1:setting$calib.size]

#Calibration truth
calib.truth <- changeTopCOD(calib.final$fakeCause)

### MLE Methods
set.seed(123)
tariff.mle <- mle.calibration(tariff.train.cod.test, tariff.train.cod.calib, calib.truth, causes)

set.seed(123)
insilico.mle <- mle.calibration(insilico.train.cod.test, insilico.train.cod.calib, calib.truth, causes)

### ReVAMP
epsilon <- .001
alpha <- .001
beta <- .001
tau <- .1
tau.vec <- rep(tau, length(causes))
delta <- 1
gamma.init <- 1
ndraws <- 50E3
nchains <- 3

### ReVAMP with Tariff
set.seed(123)
revamp.seeds <- sample(1e6, nchains, replace = F)
tariff.revamp <- lapply(1:nchains, function(i) {
  set.seed(revamp.seeds[i])
  revamp.sampler(test.cod = tariff.train.cod.test, calib.cod = tariff.train.cod.calib,
                 calib.truth = calib.truth, causes = causes,
                 epsilon = epsilon, alpha=alpha, beta=beta,
                 tau.vec=tau.vec, delta=delta,
                 gamma.init=gamma.init, ndraws = ndraws)
})

### ReVAMP with InSilicoVA
set.seed(123)
revamp.seeds <- sample(1e6, nchains, replace = F)
insilico.revamp <- lapply(1:nchains, function(i) {
  set.seed(revamp.seeds[i])
  revamp.sampler(test.cod = insilico.train.cod.test, calib.cod = insilico.train.cod.calib,
                 calib.truth = calib.truth, causes = causes,
                 epsilon = epsilon, alpha=alpha, beta=beta,
                 tau.vec=tau.vec, delta=delta,
                 gamma.init=gamma.init, ndraws = ndraws)
})

### ensemble revamp
calib.cod.tariff <- tariff.train.cod[(1:setting$calib.size)]
test.cod.tariff <- tariff.train.cod[-(1:setting$calib.size)]

calib.cod.insilico <- insilico.train.cod[(1:setting$calib.size)]
test.cod.insilico <- insilico.train.cod[-(1:setting$calib.size)]

test.cod.mat <- matrix(c(test.cod.tariff, test.cod.insilico), ncol = 2)
calib.cod.mat <- matrix(c(calib.cod.tariff, calib.cod.insilico), ncol = 2)


set.seed(123)
revamp.seeds <- sample(1e6, nchains, replace = F)
ensemble.revamp <- lapply(1:nchains, function(i) {
  set.seed(revamp.seeds[i])
  revamp.ensemble.sampler(test.cod.mat = test.cod.mat, calib.cod.mat = calib.cod.mat,
                          calib.truth = calib.truth, causes = causes,
                          epsilon = epsilon, alpha=alpha, beta=beta,
                          tau.vec=tau.vec, delta=delta,
                          gamma.init=gamma.init, ndraws = ndraws)
})


### ensemble_lite
set.seed(123)
revamp.seeds <- sample(1e6, nchains, replace = F)
ensemble.lite.revamp <- lapply(1:nchains, function(i) {
  set.seed(revamp.seeds[i])
  revamp.ensemble.lite.sampler(test.cod.mat = test.cod.mat, calib.cod.mat = calib.cod.mat,
                          calib.truth = calib.truth, causes = causes,
                          epsilon = epsilon, alpha=alpha, beta=beta,
                          tau.vec=tau.vec, delta=delta,
                          gamma.init=gamma.init, ndraws = ndraws)
})


###### Informed calibration
if(setting$Mtype == 2) {
    X <- matrix(c(1, 0, 0, 0,
                  1, 1, 0, 0,
                  0, 0, 1, 1,
                  0, 0, 1, 1),
                4, 4, byrow = TRUE)
    ### ReVAMP with Tariff
    set.seed(123)
    revamp.seeds <- sample(1e6, nchains, replace = F)
    tariff.revamp.informed <- lapply(1:nchains, function(i) {
        set.seed(revamp.seeds[i])
        revamp.sampler(test.cod = tariff.train.cod.test, calib.cod = tariff.train.cod.calib,
                       calib.truth = calib.truth, causes = causes,
                       epsilon = epsilon, alpha=alpha, beta=beta,
                       tau.vec=tau.vec, delta=delta,
                       gamma.init=gamma.init, ndraws = ndraws,
                       informed = TRUE, X = X)
    })
    ### ReVAMP with InSilicoVA
    set.seed(123)
    revamp.seeds <- sample(1e6, nchains, replace = F)
    insilico.revamp.informed <- lapply(1:nchains, function(i) {
        set.seed(revamp.seeds[i])
        revamp.sampler(test.cod = insilico.train.cod.test, calib.cod = insilico.train.cod.calib,
                       calib.truth = calib.truth, causes = causes,
                       epsilon = epsilon, alpha=alpha, beta=beta,
                       tau.vec=tau.vec, delta=delta,
                       gamma.init=gamma.init, ndraws = ndraws,
                       informed = TRUE, X = X)
    })
    
    ### ensemble lite revamp informed
    calib.cod.tariff <- tariff.train.cod[(1:setting$calib.size)]
    test.cod.tariff <- tariff.train.cod[-(1:setting$calib.size)]
    
    calib.cod.insilico <- insilico.train.cod[(1:setting$calib.size)]
    test.cod.insilico <- insilico.train.cod[-(1:setting$calib.size)]
    
    test.cod.mat <- matrix(c(test.cod.tariff, test.cod.insilico), ncol = 2)
    calib.cod.mat <- matrix(c(calib.cod.tariff, calib.cod.insilico), ncol = 2)
    
    
    ### ensemble_lite
    set.seed(123)
    revamp.seeds <- sample(1e6, nchains, replace = F)
    ensemble.lite.revamp.informed <- lapply(1:nchains, function(i) {
        set.seed(revamp.seeds[i])
        revamp.ensemble.lite.sampler(test.cod.mat = test.cod.mat, calib.cod.mat = calib.cod.mat,
                                     calib.truth = calib.truth, causes = causes,
                                     epsilon = epsilon, alpha=alpha, beta=beta,
                                     tau.vec=tau.vec, delta=delta,
                                     gamma.init=gamma.init, ndraws = ndraws,
                                     informed = TRUE, X = X)
    })
    ### Get posterior draws for CSMF parameters
    burnin <- 10E3
    thin <- 10
    
    tariff.revamp.csmf.list <- lapply(1:nchains, function(i) {
        revampCSMF(tariff.revamp[[i]], burnin = burnin, thin = thin)
    })
    tariff.revamp.csmf.df <- do.call(rbind, tariff.revamp.csmf.list)
    
    insilico.revamp.csmf.list <- lapply(1:nchains, function(i) {
        revampCSMF(insilico.revamp[[i]], burnin = burnin, thin = thin)
    })
    insilico.revamp.csmf.df <- do.call(rbind, insilico.revamp.csmf.list)
    
    ensemble.revamp.csmf.list <- lapply(1:nchains, function(i) {
        revampCSMF(ensemble.revamp[[i]], burnin = burnin, thin = thin)
    })
    ensemble.revamp.csmf.df <- do.call(rbind, ensemble.revamp.csmf.list)
    
    ensemble.lite.revamp.csmf.list <- lapply(1:nchains, function(i) {
        revampCSMF(ensemble.lite.revamp[[i]], burnin = burnin, thin = thin)
    })
    ensemble.lite.revamp.csmf.df <- do.call(rbind, ensemble.lite.revamp.csmf.list)
    
    tariff.revamp.informed.csmf.list <- lapply(1:nchains, function(i) {
        revampCSMF(tariff.revamp.informed[[i]], burnin = burnin, thin = thin)
    })
    tariff.revamp.informed.csmf.df <- do.call(rbind, tariff.revamp.informed.csmf.list)
    
    insilico.revamp.informed.csmf.list <- lapply(1:nchains, function(i) {
        revampCSMF(insilico.revamp.informed[[i]], burnin = burnin, thin = thin)
    })
    insilico.revamp.informed.csmf.df <- do.call(rbind, insilico.revamp.informed.csmf.list)
    
    ensemble.lite.revamp.informed.csmf.list <- lapply(1:nchains, function(i) {
        revampCSMF(ensemble.lite.revamp.informed[[i]], burnin = burnin, thin = thin)
    })
    ensemble.lite.revamp.informed.csmf.df <- do.call(rbind, ensemble.lite.revamp.informed.csmf.list)
    
    
    ### Function to get CSMF estimate for ReVAMP
    revampCSMFMeanEstimate <- function(revamp.csmf.df, causes) {
        csmf.init <- sapply(causes, function(c) {
            revamp.cause <- revamp.csmf.df[revamp.csmf.df$cause == c,]
            return(mean(revamp.cause$p))
        })
        return(csmf.init)
    }
    
    ### Function to get CSMF estimate from openVA
    openVACSMF <- function(topcod, causes) {
        csmf <- sapply(causes, function(c) mean(topcod == c))
        return(csmf)
    }
    
    ### Get CSMF estimates 
    tariff.revamp.csmf <- revampCSMFMeanEstimate(tariff.revamp.csmf.df, causes)
    insilico.revamp.csmf <- revampCSMFMeanEstimate(insilico.revamp.csmf.df, causes)
    ensemble.revamp.csmf <- revampCSMFMeanEstimate(ensemble.revamp.csmf.df, causes)
    ensemble.lite.revamp.csmf <- revampCSMFMeanEstimate(ensemble.lite.revamp.csmf.df, causes)
    tariff.revamp.informed.csmf <- revampCSMFMeanEstimate(tariff.revamp.informed.csmf.df, causes)
    insilico.revamp.informed.csmf <- revampCSMFMeanEstimate(insilico.revamp.informed.csmf.df, causes)
    ensemble.lite.revamp.informed.csmf <- revampCSMFMeanEstimate(ensemble.lite.revamp.informed.csmf.df, causes)
    tariff.mle.csmf <- tariff.mle$p
    insilico.mle.csmf <- insilico.mle$p
    tariff.train.csmf <- openVACSMF(tariff.train.cod.test, causes)
    tariff.train.calib.csmf <- openVACSMF(tariff.train.calib.cod, causes)
    tariff.calib.csmf <- openVACSMF(tariff.calib.cod, causes)
    insilico.train.csmf <- openVACSMF(insilico.train.cod.test, causes)
    insilico.train.calib.csmf <- openVACSMF(insilico.train.calib.cod, causes)
    insilico.calib.csmf <- openVACSMF(insilico.calib.cod, causes)
    
    ### CSMF  data frame
    methods <- c("tariff_revamp",
                 "tariff_revamp_informed",
                 "insilico_revamp",
                 "insilico_revamp_informed",
                 "tariff_mle",
                 "insilico_mle",
                 "tariff_train",
                 "tariff_train_and_calib",
                 "tariff_calib",
                 "insilico_train",
                 "insilico_train_and_calib",
                 "insilico_calib",
                 "ensemble_revamp",
                 "ensemble_lite_revamp",
                 "ensemble_lite_revamp_informed")
    csmf.df <- data.frame(csmf.est = c(tariff.revamp.csmf,
                                       tariff.revamp.informed.csmf,
                                       insilico.revamp.csmf,
                                       insilico.revamp.informed.csmf,
                                       tariff.mle.csmf,
                                       insilico.mle.csmf,
                                       tariff.train.csmf,
                                       tariff.train.calib.csmf,
                                       tariff.calib.csmf,
                                       insilico.train.csmf,
                                       insilico.train.calib.csmf,
                                       insilico.calib.csmf,
                                       ensemble.revamp.csmf,
                                       ensemble.lite.revamp.csmf,
                                       ensemble.lite.revamp.informed.csmf),
                          cause = rep(causes, length(methods)),
                          cause.text = rep(top3cause.df$gs_text34, length(methods)),
                          true.csmf = rep(top3cause.df$ptest, length(methods)),
                          method = rep(methods, each = nrow(top3cause.df)))
    ### Now CCC 
    ### Need custom function
    ccc <- function(cod.est, cod.truth, causes) {
        C <- length(causes)
        ccc <- sapply(seq_along(causes), function(j) {
            cause.j <- causes[j]
            correct.assign <- sum(cod.est == cause.j & cod.truth == cause.j)
            total <- sum(cod.truth == cause.j)
            if(total == 0) {
                total <- 1
            }
            numerator <- (correct.assign / total) - (1 / C)
            denominator <- 1 - (1 / C)
            return(numerator / denominator)
        })
    }
    
    cod.truth <- changeTopCOD(test.final$fakeCause)
    
    ### Get CCC estimates for each cause for each method
    ### Only use first chain for this
    tariff.revamp.cod <- revampIndPredictions(tariff.revamp[[1]],
                                              test.cod = tariff.train.cod.test,
                                              causes = causes, burnin = burnin,
                                              thin = thin)$topCOD
    tariff.revamp.ccc <- ccc(tariff.revamp.cod, cod.truth, causes)
    
    tariff.revamp.informed.cod <- revampIndPredictions(tariff.revamp.informed[[1]],
                                              test.cod = tariff.train.cod.test,
                                              causes = causes, burnin = burnin,
                                              thin = thin)$topCOD
    tariff.revamp.informed.ccc <- ccc(tariff.revamp.informed.cod, cod.truth, causes)
    
    insilico.revamp.cod <- revampIndPredictions(insilico.revamp[[1]],
                                                test.cod = insilico.train.cod.test,
                                                causes = causes, burnin = burnin,
                                                thin = thin)$topCOD
    insilico.revamp.ccc <- ccc(insilico.revamp.cod, cod.truth, causes)
    
    insilico.revamp.informed.cod <- revampIndPredictions(insilico.revamp.informed[[1]],
                                                test.cod = insilico.train.cod.test,
                                                causes = causes, burnin = burnin,
                                                thin = thin)$topCOD
    insilico.revamp.informed.ccc <- ccc(insilico.revamp.informed.cod, cod.truth, causes)
    
    ensemble.revamp.cod <- revampEnsembleIndPredictions(ensemble.revamp[[1]],
                                                        test.cod.mat = test.cod.mat,
                                                        causes = causes, burnin = burnin,
                                                        thin = thin)$topCOD
    ensemble.revamp.ccc <- ccc(ensemble.revamp.cod, cod.truth, causes)
    
    
    ensemble.lite.revamp.cod <- revampEnsembleLiteIndPredictions(ensemble.lite.revamp[[1]],
                                                                 test.cod.mat = test.cod.mat,
                                                                 causes = causes, burnin = burnin,
                                                                 thin = thin)$topCOD
    ensemble.lite.revamp.ccc <- ccc(ensemble.lite.revamp.cod, cod.truth, causes)
    
    ensemble.lite.revamp.informed.cod <- revampEnsembleLiteIndPredictions(ensemble.lite.revamp.informed[[1]],
                                                                 test.cod.mat = test.cod.mat,
                                                                 causes = causes, burnin = burnin,
                                                                 thin = thin)$topCOD
    ensemble.lite.revamp.informed.ccc <- ccc(ensemble.lite.revamp.informed.cod, cod.truth, causes)
    
    
    tariff.train.ccc <- ccc(tariff.train.cod.test, cod.truth, causes)
    tariff.train.calib.ccc <- ccc(tariff.train.calib.cod, cod.truth, causes)
    tariff.calib.ccc <- ccc(tariff.calib.cod, cod.truth, causes)
    insilico.train.ccc <- ccc(insilico.train.cod.test, cod.truth, causes)
    insilico.train.calib.ccc <- ccc(insilico.train.calib.cod, cod.truth, causes)
    insilico.calib.ccc <- ccc(insilico.calib.cod, cod.truth, causes)
    
    methods.ccc <- methods[!grepl("mle", methods)]
    
    ccc.df <- data.frame(ccc.cause = c(tariff.revamp.ccc,
                                       tariff.revamp.informed.ccc,
                                       insilico.revamp.ccc,
                                       insilico.revamp.informed.ccc,
                                       tariff.train.ccc,
                                       tariff.train.calib.ccc,
                                       tariff.calib.ccc,
                                       insilico.train.ccc,
                                       insilico.train.calib.ccc,
                                       insilico.calib.ccc,
                                       ensemble.revamp.ccc,
                                       ensemble.lite.revamp.ccc,
                                       ensemble.lite.revamp.informed.ccc),
                         cause = rep(causes, length(methods.ccc)),
                         cause.text = rep(top3cause.df$gs_text34, length(methods.ccc)),
                         true.ccc = rep(top3cause.df$ptest, length(methods.ccc)),
                         method = rep(methods.ccc, each = nrow(top3cause.df)))
    
    avg.ccc.df <- data.frame(ccc.mean = c(mean(tariff.revamp.ccc),
                                          mean(tariff.revamp.informed.ccc),
                                          mean(insilico.revamp.ccc),
                                          mean(insilico.revamp.informed.ccc),
                                          mean(tariff.train.ccc),
                                          mean(tariff.train.calib.ccc),
                                          mean(tariff.calib.ccc),
                                          mean(insilico.train.ccc),
                                          mean(insilico.train.calib.ccc),
                                          mean(insilico.calib.ccc),
                                          mean(ensemble.revamp.ccc),
                                          mean(ensemble.lite.revamp.ccc),
                                          mean(ensemble.lite.revamp.informed.ccc)),
                             method = methods.ccc)
} else {
    ### Get posterior draws for CSMF parameters
    burnin <- 10E3
    thin <- 10
    
    tariff.revamp.csmf.list <- lapply(1:nchains, function(i) {
        revampCSMF(tariff.revamp[[i]], burnin = burnin, thin = thin)
    })
    tariff.revamp.csmf.df <- do.call(rbind, tariff.revamp.csmf.list)
    
    insilico.revamp.csmf.list <- lapply(1:nchains, function(i) {
        revampCSMF(insilico.revamp[[i]], burnin = burnin, thin = thin)
    })
    insilico.revamp.csmf.df <- do.call(rbind, insilico.revamp.csmf.list)
    
    ensemble.revamp.csmf.list <- lapply(1:nchains, function(i) {
        revampCSMF(ensemble.revamp[[i]], burnin = burnin, thin = thin)
    })
    ensemble.revamp.csmf.df <- do.call(rbind, ensemble.revamp.csmf.list)
    
    ensemble.lite.revamp.csmf.list <- lapply(1:nchains, function(i) {
        revampCSMF(ensemble.lite.revamp[[i]], burnin = burnin, thin = thin)
    })
    ensemble.lite.revamp.csmf.df <- do.call(rbind, ensemble.lite.revamp.csmf.list)
    
    
    ### Function to get CSMF estimate for ReVAMP
    revampCSMFMeanEstimate <- function(revamp.csmf.df, causes) {
        csmf.init <- sapply(causes, function(c) {
            revamp.cause <- revamp.csmf.df[revamp.csmf.df$cause == c,]
            return(mean(revamp.cause$p))
        })
        return(csmf.init)
    }
    
    ### Function to get CSMF estimate from openVA
    openVACSMF <- function(topcod, causes) {
        csmf <- sapply(causes, function(c) mean(topcod == c))
        return(csmf)
    }
    
    ### Get CSMF estimates 
    tariff.revamp.csmf <- revampCSMFMeanEstimate(tariff.revamp.csmf.df, causes)
    insilico.revamp.csmf <- revampCSMFMeanEstimate(insilico.revamp.csmf.df, causes)
    ensemble.revamp.csmf <- revampCSMFMeanEstimate(ensemble.revamp.csmf.df, causes)
    ensemble.lite.revamp.csmf <- revampCSMFMeanEstimate(ensemble.lite.revamp.csmf.df, causes)
    tariff.mle.csmf <- tariff.mle$p
    insilico.mle.csmf <- insilico.mle$p
    tariff.train.csmf <- openVACSMF(tariff.train.cod.test, causes)
    tariff.train.calib.csmf <- openVACSMF(tariff.train.calib.cod, causes)
    tariff.calib.csmf <- openVACSMF(tariff.calib.cod, causes)
    insilico.train.csmf <- openVACSMF(insilico.train.cod.test, causes)
    insilico.train.calib.csmf <- openVACSMF(insilico.train.calib.cod, causes)
    insilico.calib.csmf <- openVACSMF(insilico.calib.cod, causes)
    
    ### CSMF  data frame
    methods <- c("tariff_revamp",
                 "insilico_revamp",
                 "tariff_mle",
                 "insilico_mle",
                 "tariff_train",
                 "tariff_train_and_calib",
                 "tariff_calib",
                 "insilico_train",
                 "insilico_train_and_calib",
                 "insilico_calib",
                 "ensemble_revamp",
                 "ensemble_lite_revamp")
    csmf.df <- data.frame(csmf.est = c(tariff.revamp.csmf,
                                       insilico.revamp.csmf,
                                       tariff.mle.csmf,
                                       insilico.mle.csmf,
                                       tariff.train.csmf,
                                       tariff.train.calib.csmf,
                                       tariff.calib.csmf,
                                       insilico.train.csmf,
                                       insilico.train.calib.csmf,
                                       insilico.calib.csmf,
                                       ensemble.revamp.csmf,
                                       ensemble.lite.revamp.csmf),
                          cause = rep(causes, length(methods)),
                          cause.text = rep(top3cause.df$gs_text34, length(methods)),
                          true.csmf = rep(top3cause.df$ptest, length(methods)),
                          method = rep(methods, each = nrow(top3cause.df)))
    ### Now CCC 
    ### Need custom function
    ccc <- function(cod.est, cod.truth, causes) {
        C <- length(causes)
        ccc <- sapply(seq_along(causes), function(j) {
            cause.j <- causes[j]
            correct.assign <- sum(cod.est == cause.j & cod.truth == cause.j)
            total <- sum(cod.truth == cause.j)
            if(total == 0) {
                total <- 1
            }
            numerator <- (correct.assign / total) - (1 / C)
            denominator <- 1 - (1 / C)
            return(numerator / denominator)
        })
    }
    
    cod.truth <- changeTopCOD(test.final$fakeCause)
    
    ### Get CCC estimates for each cause for each method
    ### Only use first chain for this
    tariff.revamp.cod <- revampIndPredictions(tariff.revamp[[1]],
                                              test.cod = tariff.train.cod.test,
                                              causes = causes, burnin = burnin,
                                              thin = thin)$topCOD
    tariff.revamp.ccc <- ccc(tariff.revamp.cod, cod.truth, causes)
    
    
    insilico.revamp.cod <- revampIndPredictions(insilico.revamp[[1]],
                                                test.cod = insilico.train.cod.test,
                                                causes = causes, burnin = burnin,
                                                thin = thin)$topCOD
    insilico.revamp.ccc <- ccc(insilico.revamp.cod, cod.truth, causes)
    
    
    ensemble.revamp.cod <- revampEnsembleIndPredictions(ensemble.revamp[[1]],
                                                        test.cod.mat = test.cod.mat,
                                                        causes = causes, burnin = burnin,
                                                        thin = thin)$topCOD
    ensemble.revamp.ccc <- ccc(ensemble.revamp.cod, cod.truth, causes)
    
    
    ensemble.lite.revamp.cod <- revampEnsembleLiteIndPredictions(ensemble.lite.revamp[[1]],
                                                                 test.cod.mat = test.cod.mat,
                                                                 causes = causes, burnin = burnin,
                                                                 thin = thin)$topCOD
    ensemble.lite.revamp.ccc <- ccc(ensemble.lite.revamp.cod, cod.truth, causes)
    
    
    tariff.train.ccc <- ccc(tariff.train.cod.test, cod.truth, causes)
    tariff.train.calib.ccc <- ccc(tariff.train.calib.cod, cod.truth, causes)
    tariff.calib.ccc <- ccc(tariff.calib.cod, cod.truth, causes)
    insilico.train.ccc <- ccc(insilico.train.cod.test, cod.truth, causes)
    insilico.train.calib.ccc <- ccc(insilico.train.calib.cod, cod.truth, causes)
    insilico.calib.ccc <- ccc(insilico.calib.cod, cod.truth, causes)
    
    methods.ccc <- methods[!grepl("mle", methods)]
    
    ccc.df <- data.frame(ccc.cause = c(tariff.revamp.ccc,
                                       insilico.revamp.ccc,
                                       tariff.train.ccc,
                                       tariff.train.calib.ccc,
                                       tariff.calib.ccc,
                                       insilico.train.ccc,
                                       insilico.train.calib.ccc,
                                       insilico.calib.ccc,
                                       ensemble.revamp.ccc,
                                       ensemble.lite.revamp.ccc),
                         cause = rep(causes, length(methods.ccc)),
                         cause.text = rep(top3cause.df$gs_text34, length(methods.ccc)),
                         true.ccc = rep(top3cause.df$ptest, length(methods.ccc)),
                         method = rep(methods.ccc, each = nrow(top3cause.df)))
    
    avg.ccc.df <- data.frame(ccc.mean = c(mean(tariff.revamp.ccc),
                                          mean(insilico.revamp.ccc),
                                          mean(tariff.train.ccc),
                                          mean(tariff.train.calib.ccc),
                                          mean(tariff.calib.ccc),
                                          mean(insilico.train.ccc),
                                          mean(insilico.train.calib.ccc),
                                          mean(insilico.calib.ccc),
                                          mean(ensemble.revamp.ccc),
                                          mean(ensemble.lite.revamp.ccc)),
                             method = methods.ccc)
}


csmf.acc.df <-
    csmf.df %>%
    group_by(method) %>%
    summarize(csmf.acc = getCSMF_accuracy(csmf.est, true.csmf))


results <- list(topcause.df = top3cause.df,
                csmf = csmf.df,
                csmf.acc = csmf.acc.df,
                csmf.test.train = getCSMF_accuracy(top3cause.df$ptrain,top3cause.df$ptest),
                csmf.test.calib = getCSMF_accuracy(top3cause.df$pcalib,top3cause.df$ptest),
                csmf.test.train.final = getCSMF_accuracy(top3cause.df$ptrain.final,top3cause.df$ptest.final),
                csmf.test.calib.final = getCSMF_accuracy(top3cause.df$pcalib.final,top3cause.df$ptest.final),
                ccc.cause = ccc.df,
                avg.ccc = avg.ccc.df)

results.file <- file.path(sim.dir, "results.rds")
saveRDS(results, results.file)

# csmf.acc.df[order(csmf.acc.df$method),]
# avg.ccc.df[order(avg.ccc.df$method),]
quit('no')