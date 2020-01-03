library(openVA)
library(gridExtra)

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

### First plot InSilico & Tariff trained on all countries

country.data <- child.clean[countries == "Tanzania",]
train.data <- child.clean[countries != "Tanzania",]
top.cod <- names(sort(table(country.data$Cause[country.data$Cause != 14]), decreasing = TRUE)[1:7])
top.cause.df <- cause.df[cause.df$va34 %in% top.cod,]
### Add row for other
top.cause.df$gs_text34 <- as.character(top.cause.df$gs_text34)
top.cause.df <- rbind(top.cause.df, c("Other", 99))

set.seed(123)
tariff <- codeVA(data = country.data,
                       data.type = "customize", model = "Tariff",
                       data.train = train.data, causes.train = "Cause")

set.seed(123)
insilico <- codeVA(data = country.data,
                         data.type = "customize", model = "InSilicoVA",
                         data.train = train.data, causes.train = "Cause",
                         jump.scale = 0.05, Nsim=5000, auto.length = FALSE)


changeTopCOD <- function(topcod, topcause.df = top.cause.df) {
    topcod <- as.character(topcod)
    ncauses <- sum(topcause.df$va34 != "99")
    topcod <- ifelse(topcod %in% topcause.df$va34[1:ncauses], topcod, "99")
    topcod <- sapply(topcod, function(c) top.cause.df$gs_text34[top.cause.df$va34==c])
    return(topcod)
}

cod.tariff <- tariff$causes.test[,2]
cod.insilico <- getTopCOD(insilico)[,2]

cod.tariff <- changeTopCOD(cod.tariff)
cod.insilico <- changeTopCOD(cod.insilico)

library(lattice)
library(gplots)
cod_names <- top.cause.df$gs_text34
cod.true <- factor(changeTopCOD(country.data$Cause), levels = cod_names)
cod.all <- list(Tariff = cod.tariff, InSilicoVA = cod.insilico) 
#pdf("../visualizations/tanzania-misclassification-all.pdf")
plot.tanzania <- lapply(1:length(cod.all), function(i){
    cod.fit <- factor(cod.all[[i]], levels = cod_names)
    tab <- table(cod.true, cod.fit)
    misclass <- t(scale(t(tab), center = FALSE, 
                      scale = colSums(t(tab))))
    rownames(misclass)[rownames(misclass) == "Other Cardiovascular Diseases"] <- "Other CVD"
    colnames(misclass)[colnames(misclass) == "Other Cardiovascular Diseases"] <- "Other CVD"
    acc <-  round(sum(diag(tab)) / sum(tab), 4) * 100
    more <- ifelse(i == 1, TRUE, FALSE)
    mat.plot <- t(misclass)
    mat.plot <- apply(mat.plot,2,rev)
    if(i == 1){
        print(
            levelplot(mat.plot,
                      scales=list(tck=0, x=list(rot=45, cex = .75),
                                  y = list(cex = .75)),
                      col.regions=colorpanel(11, "white", "grey10"),
                      at=seq(0, 1, len = 11),
                      #colorkey = FALSE,
                      main=paste0(names(cod.all)[i], " - Accuracy: ", acc, "%"),
                      xlab="Predicted Causes", ylab="True Causes", border = "black")
        ) 
    } else {
        print(
            levelplot(mat.plot,
                      scales=list(tck=0, x=list(rot=45, cex = .75),
                                  y = list(cex = .75)),
                      col.regions=colorpanel(11, "white", "grey10"),
                      at=seq(0, 1, len = 11),
                      main=paste0(names(cod.all)[i], " - Accuracy: ", acc, "%"),
                      xlab="Predicted Causes", ylab=NULL, border = "black")
        )
    }
    
})
pdf("tanzania-misclassification-all.pdf",
    width = 10, height = 10)
grid.arrange(plot.tanzania[[1]], plot.tanzania[[2]], nrow = 1)
dev.off()