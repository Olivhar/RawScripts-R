#chargement des packages
library(readxl)
library(writexl)
library(SensoMineR)
library(data.table)
library(tidyverse)
library(reshape) 
library(agricolae)
library(ggthemes)
library(Rmisc)     #pour la fonction summarySE()
rm(list = ls()) #vider l'environnement

#####TestCorrélation#####
cor(fr$Sucrose, fr$Acetate, method = c("pearson", "kendall", "spearman"))
cor.test(fr$Sucrose, fr$Acetate, method=c("pearson"))
cor.test(fr$Sucrose, fr$Acetate, method=c("kendall"))
cor.test(fr$Sucrose, fr$Acetate, method=c("spearman"))


##############################################
##################OU##########################
##############################################
pad <-read_xlsx("~/Dossiers Olivier HARLE/INRA/Tabtravail/AromeSucresAcides_1218.xlsx", sheet=3, col_names = T) 
pad<- filter(pad, date == 6)
pad <- select(pad,-date)
pad <- cbind(pad[,5:203],pad[,c(1:4,204:262)])
names(pad)
Suci <- pad %>%
  select(RefCIRM,Saccharose,Acetate,Lactate,pH)
# TRAITEMENT DES DONNEES POUR LES SUCRES, LES ACIDES ET LES DONNEES DE GC
Suc <- NULL
for(ref in c(18,20,23,26,34,36,67,251,254,257,258,261,313,
             653,772,777,845,855,1035,1046,1051,1053,1056,1108,
             1111,1128,1135,1358,1363,1364,1420,1490,1568,1860,
             1864,2102,2103,2104,2107,2115,2169,2183,2184,2185,2186,2210,2239,"t30","t43")){
  Suc1 <- filter(Suci, RefCIRM==c(ref))
  Suc <- rbind(Suc,Suc1)
}
Suc$RefCIRM[102:107] <- "t"
Suc <- filter(Suc, !is.na(Saccharose))

ggplot(Suc,aes(x= Saccharose, y= pH)) +
  geom_point(stat = "identity")

cor.test(Suc$Acetate,Suc$Saccharose, method=c("pearson"))
cor.test(Suc$Saccharose, Suc$pH, method=c("kendall"))
cor.test(Suc$Saccharose, Suc$pH, method=c("spearman"))

