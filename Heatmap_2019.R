rm(list = ls()) #vider l'environnement
#chargement des packages
library(gplots)  ###Pour heatmaps.2
library(pheatmap)
library(readxl)
library(writexl)
library(FactoMineR)
library(factoextra)
library(tidyverse) 
library(ggthemes)
library(extrafont)
# library(agricolae) #pour t.test par paires --> HSD.test()
# library(Rmisc)     #pour la fonction summarySE()
library(reshape)   #pour la fonction melt()
setwd("C:/Users/a0148148/Documents/Dossiers Olivier HARLE/R") #retrouve les fichiers - definition de l'espace de travail


pad <-read_xlsx("~/Dossiers Olivier HARLE/INRA/Tabtravail/AromeSucresAcides_1218.xlsx", sheet=3, col_names = T) 
# pad <-read_xlsx("~/Dossiers Olivier HARLE/INRA/Tabtravail/Ibis-juillet2018_Abarome_publi.xlsx", sheet=3, col_names = T) 
pad<- filter(pad, date == 6)
pad <- select(pad,-date)
GCi <- cbind(pad[,5:42],pad[,c(1:4,205:208)],pad[,c(209)])
# filter(GC,RefCIRM == 23)
GC <- NULL #Souches 35, 1119, 1058, 651, 254 et 2168 enlevées pour article. Souche 865 ajoutée
  for(ref in c(18,20,23,26,34,36,67,251,257,258,261,313,653,772,777,845,855,865,1035,1046,1051,1053,1056,1108,1111,1128,1135,1358,1363,1364,1420,1490,1568,1860,1864,2102,2103,2104,2107,2115,2169,2183,2184,2185,2186,2210,2239,"t30","t43")){
GC1 <- filter(GCi, RefCIRM==c(ref))
GC <- rbind(GC,GC1)
}

####Heat map####
names(GC)
GC$Name
GC <- GC[,1:44]
# GClog <- GC
# GClog[-c(39:43)] <- log10(GC[-c(39:43)]+1)
GC[-c(39:44)]<- data.frame(lapply(GC[-c(39:44)], function(x) scale(x, center = T, scale = max(x, na.rm = TRUE)/100)))
# GC[-c(39:43)]<- data.frame(lapply(GC[-c(39:43)], function(x) log10(x)))
# devstd<- data.frame(apply(GC[-c(39:43)],2,sd))
# GC[-c(39:43)]<- data.frame(lapply(GC[-c(39:43)],function(x) (x/sqrt(devstd[,1]))))


GC1 <- GC %>%
  group_by(RefCIRM,Espece,Name) %>%   #ajouter fichier si on veut toute les données 
  summarise(n=n(),
            butanal = mean(Butanal),
            `methylacetate` = mean(`methyl acetate`),
            ethylacetate = mean(`Ethylacetate-61`),
            `3-methylbutanal` = mean(Butanal3methyl),
            ethanol = mean(Ethanol),
            `butane-2,3-dione` = mean(`23Butanedione`),
            `propan-1-ol` = mean(`1Propanol`),
            `pentane-2,3-dione` = mean(`23Pentanedione`),
            hexanal = mean(Hexanal),
            `2-methylpropan-1-ol` = mean(`1Propanol2methyl`),
             `2,4-dimethylpentan-3-one` = mean(`3pentanone24dimethylou3Hexanone2methyl`),
            `heptane-2,3-dione` = mean(Acetylvaleryl),
            `1-butanol` = mean(`1Butanol`),
             `heptan-2-one` = mean(`2heptanone`),
            `2-pentylfuran` = mean(Furan2pentyl),
            `1-pentanol` = mean(`1Pentanol`),
            `3-hydroxybutan-2-one` = mean(Acetoin),
            `1-hydroxypropan-2-one` = mean(`2Propanone1hydroxy`),
            `(2E)-hept-2-enal` = mean(`2Heptenal`),
            `4-methylpentan-1-ol` = mean(`1pentanol4methyl`),
            `2-hydroxypentan-3-one` = mean(`2hydroxy3pentanone`),
             aceticacid = mean(Aceticacid),
            `2-methylthiolan-3-one` = mean(`3(2H)Thiophenonedihydro2methyl`),
            benzaldehyde = mean(Benzaldehyde),
            butanoicacid = mean(Butanoicacid),
            `3-methylbutanoicacid` = mean(Butanoic3methyl),
            `furan-2-ylmethanol` = mean(`2Furanmethanol`),
            pentanoicacid = mean(PentanoicAcid),
            hexanoicacid = mean(Hexanoicacid),
            `3-hydroxy-2-methylpyran-4-one` = mean(Maltol),
            heptanoicacid = mean(HeptanoicAcid),
            octanoicacid = mean(OctanoicAcid),
            nonanoicacid = mean(NonanoicAcid),
            `1-phenylethan-1-ol` = mean(phenylethylAlcohol),
            `2,4-dimethylbenzaldehyde` = mean(benzaldehyde24dimethyl))

# GC1 <- GC %>%     Pour faire la même chose mais en affiachant le numero CAS.
#   group_by(RefCIRM,Espece,Name) %>%
#   summarise(n=n(), 
#             `123-72-8` = mean(Butanal),
#             `79-20-9` = mean(Aceticacidmethylester),
#             `141-78-6` = mean(EthylAcetate),
#             `78-93-3` = mean(`2Butanone`),
#             `590-86-3` = mean(Butanal3methyl),
#             `64-17-5` = mean(Ethanol),
#             `431-03-8` = mean(`23Butanedione`),
#             `71-23-8` = mean(`1Propanol`),
#             `600-14-6` = mean(`23Pentanedione`),
#             `66-25-1` = mean(Hexanal),
#             `78-83-1` = mean(`1Propanol2methyl`),
#             `108-38-3` = mean(`Benzene13dimethyl`),
#             `565-80-0` = mean(`3pentanone24dimethylou3Hexanone2methyl`),
#             `96-04-8` = mean(Acetylvaleryl),
#             `71-36-3` = mean(`1Butanol`),
#             `95-47-6` = mean(ooupxylene),
#             `110-43-0` = mean(`2heptanone`),
#             `1115-11-3` = mean(`2Butenal2methyl`),
#             `925-89-3` = mean(`1propene1thiol`),
#             `3777-69-3` = mean(Furan2pentyl),
#             `71-41-0` = mean(`1Pentanol`),
#             `513-86-0` = mean(Acetoin),
#             `116-09-6` = mean(`2Propanone1hydroxy`),
#             `2463-63-0` = mean(`2Heptenal`),
#             `288-16-4` = mean(Isothiazole),
#             `626-89-1` = mean(`1pentanol4methyl`),
#             `5704-20-1` = mean(`2hydroxy3pentanone`),
#             `17094-21-2` = mean(ButanoicAcid2methyl3oxomethylester),
#             `111-27-3` = mean(Hexanol),
#             `64-19-7` = mean(Aceticacid),
#             `6838-51-3` = mean(`4octanone5hydroxy27dimethyl`),
#             `13679-85-1` = mean(`3(2H)Thiophenonedihydro2methyl`),
#             `100-52-7` = mean(Benzaldehyde),
#             `107-92-6` = mean(Butanoicacid),
#             `79-31-2` = mean(propanoicacid2methyl),
#             `503-74-2` = mean(Butanoic3methyl),
#             `98-00-0` = mean(`2Furanmethanol`),
#             `109-52-4` = mean(PentanoicAcid),
#             `15764-16-6` = mean(benzaldehydedimethyl24),
#             `142-62-1` = mean(Hexanoicacid),
#             `118-71-8` = mean(Maltol),
#             `111-14-8` = mean(HeptanoicAcid),
#             `124-07-2` = mean(OctanoicAcid),
#             `112-05-0` = mean(NonanoicAcid),
#             `60-12-8` = mean(phenylethylAlcohol))
           
GC1 <- select(GC1, -n)
ord <- hclust( dist(GC1, method = "euclidean"), method = "ward.D" )$order           
hclust( dist(GC1, method = "euclidean"), method = "ward.D" )        
GC1[which(GC1$RefCIRM=="2183"),]$Name <- "L.pl2183"
GC1[which(GC1$RefCIRM=="2184"),]$Name <- "L.pl2184"
GC1[which(GC1$RefCIRM=="2185"),]$Name <- "L.pl2185"
GC1[which(GC1$RefCIRM=="2186"),]$Name <- "L.pl2186"
         

# GC3 <- filter(GC1, Espece == "thermophilus")  Si ou veut obtenir uniquement les thermophilus
GC2 <- as.data.frame(GC1)
GC2 <- GC2[,-c(1,2,3)]  ###???Avec 4 pour enlever fichier
rownames(GC2) <- GC1$Name  ###???Avec $fichier quand fichier
# my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
# col_breaks = c(seq(-100,-11,length=100),  # for red
#                seq(-10,10,length=100),           # for yellow
#                seq(11,100,length=100)) 



hm1 <-pheatmap(as.matrix(GC2), col=colorRampPalette(c("navy","white","red"))(n=100),scale="row",trace="none",
               border_color = NA, legend = T,Colv=as.dendrogram(hclust(dist(t(as.matrix(GC2))))),
               main = "Heatmap of correlation between the CIRM-BIA strains and the volatil compounds produced", width = 8, height = 12)
pheatmap(as.matrix(GC2), col=colorRampPalette(c("navy","white","red"))(n=100),trace="none",
         border_color = NA, legend = T,Colv=as.dendrogram(hclust(dist(t(as.matrix(GC2))))),
         main = "Heatmap of correlation between the CIRM-BIA strains and the volatil compounds produced", width = 8, height = 12)

jpeg("HM1119.jpeg", width = 8, height = 12, units = 'in', res = 1000)

pheatmap(as.matrix(GC2), col=colorRampPalette(c("navy","white","red"))(n=100),trace="none",
         border_color = NA, legend = T,Colv=as.dendrogram(hclust(dist(t(as.matrix(GC2))))),
         width =16, height =24)

dev.off()
###Calcul des ratios pour la table de la publi###
Res <- NA
Ta <- NULL
for(vol in 4:41){
  max <- max(GC1[,vol])
  min <- min(GC1[,vol])
  nmax <- which.max(unlist(GC1[,vol]))
  wmax<- GC1$RefCIRM[nmax]
  Mratio <- max/min
  Mratiot <- as.numeric(max/filter(GC1, RefCIRM == "t43")[,vol])
  Res<- cbind.data.frame(min,max,wmax,Mratio,Mratiot)
  rownames(Res) <- names(GC1)[vol]
  Ta <- rbind(Ta,Res)
}
names(GC1)

###ACP avec les données de GC###
A<-PCA(GC3, ncp=6, quali.sup = c(1:4) ) 
barplot(A$eig[,2],names=1:nrow(A$eig))
fviz_pca_biplot(A, axes = c(1, 2), element = "var")

a1 <- fviz_pca_ind(A, # Montre les points seulement (mais pas le "text")
                   col.ind = GC1$name, # colorer by groups
                   addEllipses = F)
# Ellipses de concentration
coordEsp <- A$quali.sup$coord[c(66:78),]
coordRef <- A$quali.sup$coord[c(1:65),]
fviz_add(a1,coordRef , color = coordEsp, repel=T)


plot.PCA(A, axes = c(1, 2), choix = "ind", invisible="quali",habillage=3, cex = 0.8)  #Individus coloré par type: Thermo ou meso
plot.PCA(A, axes = c(1, 2), choix = "var",invisible="var", cex = 0.8,select="cos2 0.5", autoLab = "auto")  #Représentation des variables illustratives (odeurs sniffings)
plot.PCA(A, axes = c(1, 2), choix = "var",invisible="quanti.sup", cex = 0.8, select= "cos2 0.5") #Représentation des molécules

plot.PCA(A, axes = c(3, 4), choix = "ind", cex = 0.8, invisible="quali",habillage=47,select="cos2 0.1")
plot.PCA(A, axes = c(3, 4), choix = "var",invisible="var", cex = 0.8,select="cos2 0.03")
plot.PCA(A, axes = c(3, 4), choix = "var",invisible="quanti.sup",cex = 0.8, select= "cos2 0.1")

plot.PCA(A, axes = c(5, 6), choix = "ind", cex = 0.7, invisible="quali",habillage=47,select="cos2 0.1")
plot.PCA(A, axes = c(5, 6), choix = "var",invisible="var", cex = 0.7,select="cos2 0.04")
plot.PCA(A, axes = c(5, 6), choix = "var",invisible="quanti.sup",cex = 0.7, select= "cos2 0.1")
