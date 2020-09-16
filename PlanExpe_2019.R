#### Réalisation d'un plan d'expériences pour 100 produits et 40 juges ne pouvant pas en déguster plus de 10
rm(list = ls()) #vider l'environnement
#chargement des packages
library(writexl)
library(SensoMineR)
planoptimal <- optimaldesign(nbProd=8,nbPanelist =8, nbProdByPanelist = 8)
write_xlsx(as.data.frame(planoptimal$design), path ="~/Dossiers Olivier HARLE/INRA/Souches/TabR/planoptimal20-20-15.xlsx")

library(DoE.base)
plan1<-fac.design(nlevels=3,nfactors=2,factor.names=c("consor","ratios"), replications=1,repeat.only=FALSE, randomize=FALSE)

?design
