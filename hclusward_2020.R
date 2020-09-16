rm(list = ls()) #vider l'environnement
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
library(readxl)
library(tidyverse)
library(ggthemes)
library(agricolae) #pour t.test par paires --> HSD.test()
library(Rmisc)     #pour la fonction summarySE()
# library(reshape)   #pour la fonction melt()
setwd("C:/Users/a0148148/Documents/Dossiers Olivier HARLE/R") #retrouve les fichiers - definition de l'espace de travail

panG <- read_xlsx("~/Dossiers Olivier HARLE/INRA/EtudeConsortia/ARN/search_LB_Microcyc.xlsx")
B45B50 <-read.table("~/Dossiers Olivier HARLE/INRA/EtudeConsortia/ARN/DESeq2_B/reports_01/e1/tables/B45vsB50.complete.txt",header=T)
B45B50 <- (B45B50[,1:54])
summary(B45B50)

for(n in 2:ncol(B45B50)){
  newcol <- B45B50[,n][match(panG$ID,B45B50$Id)]
  panG <- cbind(panG,newcol)
}
names(panG)[11:ncol(panG)] <- names(B45B50)[2:ncol(B45B50)]
# B45B50_712<-B45B50[grep("met712",row.names(B45B50)),]
# B45B50_715<-B45B50[grep("met715",row.names(B45B50)),]
panG <- panG[,c(8,11:ncol(panG))] %>%
  group_by(Pathway) %>%
  dplyr::summarise_all(.funs= c(mean="mean"))
row.names(panG) <- panG$Pathway

plot(hclust(as.dist(1-cor(t(panG[,-1]))),method="ward.D"),hang=-1)
dev.off()
