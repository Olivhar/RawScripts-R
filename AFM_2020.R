rm(list = ls()) #vider l'environnement
#chargement des packages
library(readxl)
library(writexl)
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(tidyverse) 
library(ggthemes) ##permet d'utiliser %>%
# library(reshape)
library(ggplot2)
library(DESeq2)         ## RNA-seq differential analysis
library(edgeR)          ## RNA-seq differential analysis
library(limma)          ## RNA-seq differential analysis via limma-voom
library(RColorBrewer)   ## Colors for plotting
library(venn)           ## Venn Diagrams
library(HTSFilter)       ## To remove not expressed genes
####	  Define working directory
setwd("C:/Users/a0148148/Documents/Dossiers Olivier HARLE/R") #retrouve les fichiers - definition de l'espace de travail
# tab<-read.delim("C:/Users/a0148148/Documents/Dossiers Olivier HARLE/INRA/EtudeConsortia/ADN/GFFcontigs/Galaxy1997-[Lactobacillus_plantarum_CIRM-BIA_777p.gff].gff3", header=F, stringsAsFactors=FALSE) 

#Definition des fonctions
erreur.std <- function(x)
{ x <- na.omit(x)
ET <- sd(x)/sqrt(length(x))
return(ET)
}
#____________, lecture du fichier resultat et verifications
count1 <-read_xls("~/Dossiers Olivier HARLE/INRA/EtudeConsortia/ARN/counts/counts130220/counts1702.xls", sheet=1, col_names = T) 
count1 <- as.data.frame(count1)
rownames(count1)<-count1$ID
a <- select(count1,ID,db_xref,ec_number,`function`,gene,locus_tag,note,product,gn_start,gn_end,strand,gnsz)[-c(1:4),]
count <- select(count1,-c(ID,db_xref,ec_number,`function`,gene,locus_tag,note,product,gn_start,gn_end,strand,gnsz))
count<- count[-c(0:4),]
count <- select(count,B2a,B2b,B2c,B3a,B3b,B3c,B4a,B4b,B4c,B6a,B6b,B6c,
                S2a,S2b,S2c,S2e,S3a,S3b2110,S3c,S3d,S3e,S4b,S4e,S4d,S4e1510,S4f,S5a,S5b,S5c, #sans S3b
                M2a,M2b,M2b1810,M2c,M2c1810,M3a,M3a1810,M3b,M3c,M3c1810,M4a,M4a1810,M4b,M4b1810,M4c,M4f,M5a,M5a1810,M5b,M5b1810,M5c)
for(i in 1:ncol(count)){
  count[,i] <- as.numeric(count[,i])
}
# tcount1 <- t(count1)
# tcount1 <- as.data.frame(tcount1)
# tcount1 <- select(tcount1,-LBPLA777__no_feature,-LBPLA777__ambiguous,-LBPLA777__too_low_aQual,-LBPLA777__not_aligned,-LBPLA777__alignment_not_unique,-LBPLA777_sumcount) %>%
# select(-LBDEL865__no_feature,-LBDEL865__ambiguous,-LBDEL865__too_low_aQual,-LBDEL865__not_aligned,-LBDEL865__alignment_not_unique,-LBDEL865_sumcount)
# tcount1 <- tcount1[1:10,1:10]  ####pour test####
# summary(tcount1[1:10,1:10])
# for(i in 5:ncol(tcount1)){
  # tcount1[,i] <- as.numeric(tcount1[,i])
# }
# summary(tcount1)
# PCA(tcount1,quali.sup =c(1:4))
# tcountf <- tcount1[-which(rownames(tcount1)=="S3b"),]  ###ON enl�ve S3b qui apparait comme un T
# PCA(tcountf,quali.sup =c(1:4))
# blop <- filter(tcount1,ech=="S"&(pH=="pH5"|pH=="pH6"|pH=="pH6.5")|(ech=="M"&pH=="pH5"))
# blop <- tcount1[which(tcount1$ech=="S"&(tcount1$pH=="pH5"|tcount1$pH=="pH6"|tcount1$pH=="pH6.5")|(tcount1$ech=="M"&tcount1$pH=="pH5")),]
# blopf <- blop[-which(rownames(blop)=="S3b"),] 
# PCA(blop,quali.sup =c(1:4))
# blop$essai <- rownames(blop)
# count <- t(tcount1)%>%
  # as.data.frame()
row.names((count))
head(count)
dim(count)
str(count)
### rowmeans, colmeans, rowsum... Pour faire les sommes ou moyennes par lignes ou colonnes
head(rowMeans(count)) #Permet de connaitre le nombre de read moyen par g�nes
mean(rowMeans(count)) # Permet de connaitre le nombre de read moyen pour l'ensemble des g�nes
median(rowMeans(count)) 
#----------------------------------------------------------------------------------#
#  A2 Définition des facteurs/modalités et couleurs associées
#----------------------------------------------------------------------------------#
#### vector "group"
data<- count
group <- substring (colnames(data),1,2) ; group ; table(group)  #substring() pour prendre les 2 premiers caract�res de chaque colonnes
#### vector "couleur" => go the website: http://colorbrewer2.org
couleur<-as.factor(group)  ; couleur ; levels(couleur) <- c("gold3","gold2","gold1","gold","darkorange3","darkorange2","darkorange1","darkorange","red3","red2","red1","red")   #,"#f7fcb9","#fed98e")
couleur<-as.character(couleur) ; couleur
#----------------------------------------------------------------------------------#
# A3 exploration des données brutes 
#----------------------------------------------------------------------------------#
FSJ<-"FSJ" ; 
dim(data)
summary(colSums(data, na.rm=T))   #3304624 de read fait en moyenne par �chantillons (pour les analyse stats; c'est bien entre 25 et 50 millions de read)

#### plot des tailles des libraries - reads alignés sur le genome - (barplot)
par(mfcol=c(1,1)) ; par(mar=c(7,7,3,2))    ### mfcol -> pour mettre plusieurs graph par fen�tres, par(mar=c(,,,)) pour pr�ciser les marges 
barplot(colSums(data, na.rm=T),  col=couleur,  main = paste(FSJ," ", nrow(data)," genes - (readCount) Library sizes ",sep=""), las=2, cex.axis=1.1, cex.names=0.7, space=c(0,1)) #names.arg=""  las= -> pour �crire noms des barres en vertical
abline(h=mean(colSums(data, na.rm=T))  , col="black")
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 5288  1671743  2990896  3304624  4663834 13643536


#### Distribution des comptages par echantillon (boxplot) Un boxplot fonctionne toujours sur les colonnes
par(mfcol=c(2,1)) ; par(mar=c(6,5,3,2))
boxplot((data), col=couleur, las=1 , main = paste(FSJ," ", nrow(data)," genes -  (readCount) distribution",sep=""), cex=0.1, cex.axis=1,  xaxt="n" )
boxplot((data),col=couleur, las=2 , main = paste(FSJ," ", nrow(data)," genes -  (readCount) distribution - ZOOM",sep=""), cex=0.1, cex.axis=1, ylim=c(0, 1000) ,  xaxt="n" )   # zoom sur les ordonnées entre 0 et 500 xaxt="n" Supprime les labels si on les veut, faire las=2
#possibilit� de faire log10 mais il faut alors le faire data+1 � chaque donn�es pour ne pas avoir de 0


#### Distribution des comptages par echantillon (histogramme) 
par(mfcol=c(2,1))
hist(log10(unlist(data) + 1), col = "darkgrey", main = paste(FSJ,nrow(data),  " genes - log10(read count) "), xlab = "log10 (counts + 1)", ylab="",breaks=20, las=1)
#### Distribution pour le 1er echantillon (histogramme) 
hist(log10(data[,1] + 1), col = "darkgrey", main = paste(FSJ,nrow(data),  " genes - log10(read count) "), xlab = "log10 (counts + 1)", ylab="",breaks=20, las=1)

#### Plot des gènes les plus exprimés (somme  reads sur les 32 ech) - sélection du top 10 
## s�lection du top10
index <- order(rowMeans(data), decreasing=TRUE)[1:10]
top10<- cbind(a[index,c(1,3,5,8)] , moyRead=rowMeans(data)[index] )
top10
## poids des reads du top10 sur les reads totaux
sum(top10$moyRead) / sum (rowMeans(data)) %>%   ###Necessite magrittr si on veut utiliser %>% 
  round(2)# 8888584 / 25388801 = 35%

#### Combien de donn�es =0 sur tous les échantillons (par exemple comme le gène en position 3?  
head(data)
data[3,]
length(which(rowSums(data)==0))  ;  nrow(data)
round(length(which(rowSums(data)==0))  /  nrow(data),4  )
#### quid des gènes comme le gène n°2
data[2,]

#=> conclusion : existence de nombreux gènes non exprimés : nécessité de les enlever


#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
#              B - Normalisation des Comptages en RPKM (Normalisation la plus ancienne)
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
# 	B1 Calul des RPKM (Normalisation la plus ancienne)
#----------------------------------------------------------------------------------#
# RPKM : Read par Kilobase (du gène) et millions (de reads de l'échantillon)
#=> normalisation par la taille des gènes (RPK) et  sommes des reads par sample /10M (taille de la librairie) 
#=> read  / par tailleGene =RPK (ligne)  / tailleLibrairie /10M (somme de la colonne)
#-----------------------------------------------------------------------------------#
dim(count) 
####  RPK = normalisation des reads pour la taille du gène
## Preparation des données
size<-a$gnsz/1000; length(size) ;head(size)  # en kb   gnsz = gene size (diff�rentes mani�res de le calculer --> peut �tre fait selon base de donn�e gdf, cf. script YUNA)
summary(size)   #Certains g�nes peuvent faire 26 kb, ils sont alors plus fragile
dataETsize<-cbind.data.frame(size,count) 
dataETsize [1:3,1:4] ;dim(dataETsize)  # visualisation
## normalisation by the size
normByTailleGene<-function(x){x[2:length(x)]/x[1]} # creation de la fontion
rpk <- t(apply(dataETsize,1, normByTailleGene) )    ###??? t() pour transposer les lignes et les colonnes :D
rpk[1:2,1:3] ; dim(rpk) 
dataETsize[1,2]/dataETsize[1,1] # pour verification sur la première cellule
####  RPKM = normalisation des reads pour la taille du gène (RPK) et  sommes des reads par sample (taille de la librairie) *10M 
## data Preparation
rpkSumLibrary <- rbind.data.frame(libSize=colSums(count)/1000000, rpk) 
rpkSumLibrary[1:2,1:3] ; dim(rpkSumLibrary) 
## normalisation by the library
normByLibrary<-function(x){x[2:length(x)]/x[1]} # creation de la fontion
rpkm <- apply(rpkSumLibrary,2, normByLibrary) 
rpkm[1:2,1:3] ; dim(rpkm) 
rpkSumLibrary[2,1]/rpkSumLibrary[1,1] # pour verification sur la première cellule
#----------------------------##
#  B2 Exploration des données RPKM (plot)
#----------------------------##
#### summary 
data<-rpkm ; summary(colSums(data,na.rm=T))
data<-rpkm ; summary(data[,1],na.rm=T)
#### plot des RPKM (ajout de 0.01 car log10(0) = inf)
par(mfcol=c(2,1)) ; par(mar=c(5,4,3,2),las=1)
data<- count ; barplot(colSums(data, na.rm=T),col=couleur, main = paste(FSJ," ",nrow(data)," genes - COUNT Library sizes ",sep=""), las=2, cex.axis=0.8, cex.names=0.7) #, names.arg="")
abline(h=mean(colSums(data, na.rm=T))  , col="red")
data<-rpkm ; barplot(colSums(data, na.rm=T),col=couleur, main = paste(FSJ," ",nrow(data)," genes - RPKM Library sizes ",sep=""), las=2, cex.axis=0.8, cex.names=0.7) #, names.arg="")
abline(h=mean(colSums(data, na.rm=T))  , col="red")

#----------------------------------------------------------------------------------#
# B3 Sélection des gènes dits exprimés 
#----------------------------------------------------------------------------------#
# Selection des genes ayant un TPM >=0.1 dans au moins X% des N echantillon d'une condition
data<-rpkm
group
####  Définissons les seuils
seuilExpr <- 0.1  ### @@@    seuilExpr <- 0.1
seuilNbr <- 0.75  #  0.75 x 8 samples =6 samples
####  Simplifions dans un premier temps en prenant un tableau � 4 lignes (genes)
tt=data[1:4,]  ; tt[,1:6]
geneOK1 <-apply(tt ,1, function(x) {tapply (x, group, function(x)  sum(x>= seuilExpr) / length(x) ) }) ; head(t(geneOK1))   ###tapply(appelle les groupes des facteurs)
geneOK2<- apply(t(geneOK1), 1,function(x) sum(x >= seuilNbr)) ; head(geneOK2)
geneOK3<-geneOK2 >=1 ; geneOK3
tt[which(geneOK3),] 
####  Généralisons aux N lignes (i.e. gènes) du tableau data
geneOK1=apply(data ,1, function(x) {tapply (x, group, function(x)  sum(x>= seuilExpr) / length(x) ) }) 
head(t(geneOK1))
geneOK2<- apply(t(geneOK1), 1,function(x) sum(x >= seuilNbr))
head(geneOK2)
geneOK3<-geneOK2 >=1 # genes with TPM >= 0.1 in at least seuil%  in at least one condition
c(sum(geneOK3,na.rm=T) , nrow(data)) 
# 13321 24879= 53,5% avec seuilNb = 0.75  et seuilExpr =0.1  soit 13669 / 24879 = 55% si seuil � 60%
# possibilité de changer les seuils : 13669 / 24879 = 55% si seuilNb � 60% et seuilExpr =0.1
# possibilité de changer les seuils : 9805 / 24879 = 39,4% si seuilNb � 75% et seuilExpr =1
rpkmExpr01 = data[which(geneOK3),] # genes with RPKM >= 0.1 in at least 75%  in at least one condition
rpkmExpr01  [1:3, 1:4] ;dim(rpkmExpr01)
#rpkmExpr1 = data[which(geneOK3),] # genes with RPKM >= 1 in at least 75%  in at least one condition
#rpkmExpr1  [1:3, 1:4] ;dim(rpkmExpr1)

##  Selection of the annotation 	+ count expression    			
# rpkmExprAnnot1 = a[match(rownames(rpkmExpr1),a$ID),] ; dim(rpkmExprAnnot1)  
rpkmExprAnnot01 = a[match(rownames(rpkmExpr01),a$ID),] ; dim(rpkmExprAnnot01) ###OU rpkmExprAnnot01 = tab[which(geneOK3),]

# table(rpkmExprAnnot1$gnSimpleBiotype)    #13669 / 24879 = 55% with RPKM >= 0.1 in at least 75%  in at least one condition
table(rpkmExprAnnot01$gnSimpleBiotype)  # 9805 / 24879 = 39,4% with RPKM >= 1 in at least 75%  in at least one condition
#   lnc   mis   mtr   pcg    pse   rbz   sca   sno 
#  330      3     2      9467     3  
#  806      3     2    12011     8     2       2     1 
# => choix du seuil 0.1

#----------------------------------------------------------------------------------#
# B4 Exploration des gènes dits exprimés 
#----------------------------------------------------------------------------------#
FSJ<-"FSJ" ; 
par(mfcol=c(1,2)) ; par(mar=c(5,4,3,2),las=1)
### plot avant et apres sélection sur l'expression

##  plot des gènes exprimés
data<- rpkmExpr01
barplot(colSums(data, na.rm=T),col=couleur, main = paste(FSJ,"-",nrow(data)," expressedGenes: RPKM Library sizes ",sep=""), las=2, cex.axis=0.8, cex.names=0.7, ylim=c(0,2000000)) # names.arg=""
## plot de tous les gènes (avant sélection)
data<-count
barplot(colSums(data, na.rm=T),col=couleur, main = paste(tissue,"-",nrow(data)," - readCount Library sizes ",sep=""), las=2, cex.axis=0.8, cex.names=0.7) #, names.arg="")

### distribution d'expression
summary(rpkmExpr01[,1])
summary(rpkmExpr01[,2])
moyenne<-rowMeans(rpkmExpr01)

#----------------------------------------------------------------------------------#
# B5 Recherche d éventuels échantillons Outliers par Analyse ACP
#---------------------------------------------------------------------------------#
library(FactoMineR) # les variables devront être en colonne !
par(mfrow=c(1,1))
levels(group);group;couleur
#### RPKM
data<-rpkmExpr01
a<-log10(data+0.01) ; head(a) ;dim(a)
dataPCA = cbind.data.frame(t(a), group) ; rownames(dataPCA) ; dim(dataPCA)
l<- ncol(dataPCA) ; head(dataPCA[, c(l-1, l)])    #visualization of the two last columns
res.pca <- PCA(dataPCA, quali.sup=ncol(dataPCA), graph = F,scale.unit=T) ;head(res.pca$eig)
plot.PCA(res.pca, choix="ind", title="PCA with log10(RPKM+0.01)",habillage=ncol(dataPCA), col.hab=c("gold3","gold2","gold1","gold","darkorange3","darkorange2","darkorange1","darkorange","red3","red2","red1","red") , cex=1, 
         legend = list(bty = "y", x = "topright"))
###-> Il n'y a pas un individus qui se s�pare particuli�rement des autres, mais on peut commencer � voir que les Rm sont un peu s�par�s des Rp (pas de outlier)
#### COUNT
a<-log10(count+0.01) ; head(a) ;dim(a)
dataPCA = cbind.data.frame(t(a), group) ; rownames(dataPCA) ; dim(dataPCA)
l<- ncol(dataPCA) ; head(dataPCA[, c(l-1, l)])    #visualization of the two last columns
res.pca <- PCA(dataPCA, quali.sup=ncol(dataPCA), graph = F,scale.unit=T) ;head(res.pca$eig)
plot.PCA(res.pca, choix="ind", title="PCA with log10(read COUNT+0.01)",habillage=ncol(dataPCA), col.hab=c("gold3","gold2","gold1","gold","darkorange3","darkorange2","darkorange1","darkorange","red3","red2","red1","red") , cex=1, legend = list(bty = "y", x = "bottomright"))

## le outlier 103 nest plus retrouvé lorsqu on a normalisé les données
##### optionnel RPKM
a<-rpkmExpr ; head(a) ;dim(a)
dataPCA = cbind.data.frame(t(a), group) ; rownames(dataPCA) ; dim(dataPCA)
l<- ncol(dataPCA) ; head(dataPCA[, c(l-1, l)])    #visualization of the two last columns
res.pca <- PCA(dataPCA, quali.sup=ncol(dataPCA), graph = F,scale.unit=T) ;head(res.pca$eig)
plot.PCA(res.pca, choix="ind", title="PCA with RPKM",habillage=ncol(dataPCA), col.hab=c("gold3","gold2","gold1","gold","darkorange3","darkorange2","darkorange1","darkorange","red3","red2","red1","red") , cex=1)

#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#          C - Introduisons une nouvelle métrique d'expression normalisées- le TPM 
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#### TPM :   transcripts/gènes par millions (de reads par genes déja normalisés pour la taille)
#=> normalisation par la taille des gènes (rpk) puis somme des rpk /1M(enemble des gènes normalisés par leur taille)  :
#=> read  / par tailleGene =RPK (ligne) / som(RPK) /1M (somme des RPK par colonne)
#=> sum et moy identique per ech contrary to rpkm

#----------------------------------------------------------------------------------#
# C1 Calcul des TPM
#----------------------------------------------------------------------------------#

####  RPK = normalisation des reads pour la taille du gène
## Preparation des données
data<- count

 size<-a$gnsz /1000; length(size) ;head(size)  # en kb
 dataETsize<-cbind.data.frame(size,data) 
 dataETsize [1:3,1:4] ;dim(dataETsize)  # visualisation
## normalisation by the size
normByTailleGene<-function(x){x[2:length(x)]/x[1]} # creation de la fontion
rpk <- t (apply(dataETsize,1, normByTailleGene) )
rpk[1:2,1:3] ; dim(rpk) 
dataETsize[1,2]/dataETsize[1,1] # pour verification sur la première cellule

#### TPM  = normalisation des reads pour la taille du gène (RPK) et puis som(RPK)/M(ensemble des gènes normalisés par leur taille) 
dim(data) ; length(size)
normBySomRPK<-function(x){x/(sum(x)/1000000)}
tpm <- apply(rpk,2, normBySomRPK) ; dim(tpm) ; tpm[1:2,1:3]
tpm[1:2,1:3] ; dim(tpm) 

## plot des TPM (ajout de 0.01 car log10(0) = inf)
data<-tpm
barplot(colSums(data, na.rm=T),col=couleur, main = paste(FSJ," - TPM Library sizes ",sep=""), las=2, cex.axis=0.8, cex.names=0.7) # names.arg=""
abline(h=log10(0.01),col="red")  # ligne pour symboliser les valeurs nulles (gènes non exprimés)


#----------------------------------------------------------------------------------#
# C2 Sélection des gènes dits exprimés 
#----------------------------------------------------------------------------------#
#### Selection des genes ayant un TPM >=0.1 dans au moins X% des N echantillon d'une condition
data<-tpm
seuilExpr <- 0.1
seuilNb <- 0.75  #  0.75 x 8 samples =6 samples
geneOK1=apply(data ,1, function(x) {tapply (x, group, function(x)  sum(x>= seuilExpr) / length(x) ) }) 
head(t(geneOK1))
geneOK2<- apply(t(geneOK1), 1,function(x) sum(x >= seuilNb))
head(geneOK2)
geneOK3<-geneOK2 >=1 # genes with TPM >= 0.1 in at least seuil%  in at least one condition
c(sum(geneOK3,na.rm=T) , nrow(data)) 
# 13321 24879= 53,5% avec seuil = 0.75   # possibilité de changer de seuil : 13669 24879 = 55% si seuil à 60%

tpmExpr01 = tpm[which(geneOK3),] # genes with TPM >= 0.1 in at least 80%  in at least one condition
tpmExpr01  [1:3, 1:4] ;dim(tpmExpr01)

##  Selection of the annotation 	  +   count	
tpmExprAnnot01= a[match(rownames(tpmExpr01),a$ID),] ; dim(tpmExprAnnot01)

##plot
FSJ<-"FSJ" ; data<- tpmExpr01
par(mfrow=c(2,1)) ; par(mar=c(5,4,3,2),las=1)
barplot(colSums(data, na.rm=T),col=couleur, main = paste(FSJ,"-",nrow(data)," expressedGenes: TPM Library sizes ",sep=""), las=2, cex.axis=0.8, cex.names=0.7, ylim=c(0,2000000)) # names.arg=""
boxplot(log10(data +0.01),col=couleur,las=3 , main = paste(FSJ,"-",nrow(data)," expressedGenes: log10(TPM+0.01) distribution",sep=""), cex=0.1,  ylim=c(-2,3),cex.axis=0.6)
summary(colSums(data,na.rm=T))
abline(h=log10(0.01),col="red")

#----------------------------------------------------------------------------------#
# C3 Recherche d éventuels échantillons Outliers par Analyse ACP
#----------------------------------------------------------------------------------#
library(FactoMineR) # les variables devront être en colonne !

par(mfrow=c(1,1))
levels(group);group;couleur
data<- tpmExpr01

#PASSAGE EN LOG car ACP basé sur les correlations (données normales)
#### TPM
a<-log10(data+0.01) ; head(a) ;dim(a)
dataPCA = cbind.data.frame(t(a), group) ; rownames(dataPCA) ; dim(dataPCA)
l<- ncol(dataPCA) ; head(dataPCA[, c(l-1, l)])    #visualization of the two last columns
res.pca <- PCA(dataPCA, quali.sup=ncol(dataPCA), graph = F,scale.unit=T) ;head(res.pca$eig)
plot.PCA(res.pca,  title="PCA with log10(TPM+0.01)", choix="ind", habillage=ncol(dataPCA), col.hab=c("gold3","gold2","gold1","gold","darkorange3","darkorange2","darkorange1","darkorange","red3","red2","red1","red","black") , cex=1)

=> pas d echantillon outlier

#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#     D Comparaison des gènes exprimés selon RPKM et TPM
#----------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------#
#### Comparaison des valeurs
## summary 
data<-rpkmExpr01 ; summary(colSums(data,na.rm=T))
data<-tpmExpr01 ; summary(colSums(data,na.rm=T))
data<-rpkmExpr01 ; summary(data[,1],na.rm=T) ; dim(data)
data<-tpmExpr01 ; summary(data[,1],na.rm=T) ; dim(data)
#		    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# rpkm	 0.00     1.19     4.75    	46.79    14.40 		54697.89 
# tpm    0.00     1.62     7.30    	75.06    22.82 		91066.29





