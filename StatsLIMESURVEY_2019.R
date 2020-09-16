#Script qui m'a surtout servi à assembler les descripteurs et à les transformer en fréquence
#via la fonction textual(). 
# 
#
rm(list = ls()) #vider l'environnement
#chargement des packages
library(readxl)
library(writexl)
library(FactoMineR)
library(factoextra)
library(tidyverse) 
library(reshape)   #pour la fonction melt()
library(data.table) ###Pour la fonction grep() et pour la fonction melt() --> package utilisé pour melt par patterns (utilisé setDT si table importé d'excel)
library(agricolae) #pour t.test par paires --> HSD.test()
library(Rmisc)     #pour la fonction summarySE()
library(ggthemes) ##permet d'utiliser %>%
# ?textual
# ?descfreq
setwd("C:/Users/a0148148/Documents/Dossiers Olivier HARLE/R") #retrouve les fichiers - definition de l'espace de travail


#____________ lecture du fichier resultat et verifications
fr <-read_xlsx("~/Dossiers Olivier HARLE/Triballat/Analyses/Senso/LimeSurvey/results-survey847622.xlsx", sheet=1, col_names = T) 
fr <- fr[,c(1,2,8:length(fr))]
names(fr)[3:ncol(fr)] <- c("P1","P2","P3","P4","Apparence P1","Apparence P2","Apparence P3","Apparence P4","Odeur P1","Odeur P2","Odeur P3","Odeur P4","Texturecuil P1","Texturecuil P2","Texturecuil P3","Texturecuil P4","Texturebouche P1","Texturebouche P2","Texturebouche P3","Texturebouche P4","Gout P1","Gout P2","Gout P3","Gout P4","ApparenceNote P1","ApparenceNote P2","ApparenceNote P3","ApparenceNote P4","OdeurNote P1","OdeurNote P2","OdeurNote P3","OdeurNote P4","TexturecuilNote P1","TexturecuilNote P2","TexturecuilNote P3","TexturecuilNote P4","TextureboucheNote P1","TextureboucheNote P2","TextureboucheNote P3","TextureboucheNote P4","GoutNote P1","GoutNote P2","GoutNote P3","GoutNote P4","Noteglobale P1","Noteglobale P2","Noteglobale P3","Noteglobale P4","commentaire")
fr <- filter(fr, !is.na(P1)) #Pour enlever les cases vides
fr <- fr[-1,]###Pour enelever première ligne (test)

####Treatment START#####
names(fr)
summary(fr) ###Vérif définition classes###

###Mise en forme data FRO pour analyser les descripteurs selon des critères différents
FRO <- fr
# FRO$pH21jT <- fr$pH21jT
# FRO$pH21j37 <- fr$pH21j37
# FRO$CtoutP1 <- paste(fr$`Apparence P1`,fr$`Odeur P1`,fr$`Texturecuil P1`,fr$`Texturebouche P1`,fr$`Gout P1`,sep=";")
# FRO$CtoutP2 <- paste(fr$`Apparence P2`,fr$`Odeur P2`,fr$`Texturecuil P2`,fr$`Texturebouche P2`,fr$`Gout P2`,sep=";")
# FRO$CtoutP3 <- paste(fr$`Apparence P3`,fr$`Odeur P3`,fr$`Texturecuil P3`,fr$`Texturebouche P3`,fr$`Gout P3`,sep=";")
# FRO$CtoutP4 <- paste(fr$`Apparence P4`,fr$`Odeur P4`,fr$`Texturecuil P4`,fr$`Texturebouche P4`,fr$`Gout P4`,sep=";")
meltFRO <- melt(setDT(as.data.frame(FRO)), measure = patterns( "^P","^Apparence P","^Odeur P","^Texturecuil P","^Texturebouche P","^Gout P","^ApparenceNote P","^OdeurNote P","^TexturecuilNote P","^TextureboucheNote P","^GoutNote P","^Noteglobale P"),
                value.name = c("Produit","Apparence","Odeur","TextCuil","Textbouche","Gout","NoteApparence","NoteOdeur","NoteTextCuil","NoteTextbouche","NoteGout","Noteglobale"))
meltFRO$Ctout <- paste(meltFRO$`Apparence`,meltFRO$`Odeur`,meltFRO$`Textcuil`,meltFRO$`Textbouche`,meltFRO$`Gout`,sep=";")
names(meltFRO)
meltFRO <- filter(meltFRO, !is.na(Produit))
levels(meltFRO$Produit) <- unique(meltFRO$Produit)
##############Analyses par produit##############
###GLOBAL###
res <-   filter(meltFRO) %>%
  textual(num.text=17,contingence.by=5,sep.word =";")
descriptT <- row.names(res$nb.words)
####Odeur###
 res <- filter(FRO) %>%
   textual(num.text=11,contingence.by=9,sep.word =";")
####Aspect####
res <- filter(FRO) %>%
  textual(num.text=10,contingence.by=9,sep.word =";")

library("wordcloud")
names(as.data.frame(res$cont.table))
comparison.cloud(t(as.data.frame((res$cont.table))))
library("wordcloud2")
res 
# res$nb.word$nb.list <- res$nb.words$words
res$nb.words$words <- row.names(res$nb.words)
# sum(res$nb.words) <- res$nb.words[-1,]
figPath <-  system.file("~/Dossiers Olivier HARLE/R/yogurt-and-spoon.png",package ="wordcloud2")
wordcloud2(res$nb.words,figPath=figPath, size = 1.5,color = "skyblue")
wordcloud2(res$nb.words, color = "green", backgroundColor = "white")

# names(FRO)
cont <- as.data.frame(res$cont.table)
for(d in descriptT){     ###Permet d'ajouter toutes les colonnes des descripteurs qui ne sont pas dans les tables --> permet la fusion des descripteurs sans soucis
    if(d %in% names(cont)){
    } else {
      col <- data.frame(d=rep(0,nrow(cont)))
      names(col) <- d
      cont <- cbind(cont,col)
    }
}
res$nb.words
sum(res$nb.words$nb.list)
dim(res$cont.table)
names(cont)
descfreq(cont[,c(1:62,64:86)],proba=0.2)$"35"

#####Fusion des colonnes#####
#Global#
# cont$rustique <- cont$agrume + cont$rustique
# cont$cereal <- cont$cereal + cont$farine + cont$ble +cont$muesli+ cont$avoine
# cont$vert <- cont$vert + cont$vegetal + cont$haricot + cont$herbec + cont$herbe +cont$gazon+cont$fleurcoupe
# cont$citron <- cont$citron + cont$citrique
# cont$malt <- cont$malt + cont$malte
# cont$epais <- cont$epais + cont$ferme + cont$gel + cont$beaugel + cont$solide
# cont$synerese <- cont$synerese + cont$`synerese+`
# cont$carton <- cont$carton + cont$papier + cont$tapisserie + cont$vieux
# cont$doux <- cont$doux + cont$tresdoux
# cont$amer <- cont$amer
# cont$floral <- cont$floral + cont$fleur + cont$lys + cont$jacinthe + cont$muguet +cont$rose+cont$lila+cont$geranium+cont$violette
# cont$fromage <- cont$emmental + cont$fromage + cont$gruyere
# cont$aigre <- cont$aigre + cont$vinaigre + cont$vinaigrebalsamique + cont$vinaigredevin + cont$vinaigrette + cont$acetique
# cont$levure <- cont$levure + cont$levain
# cont$grumeau <- cont$grumeau + cont$flocule +cont$caille + cont$destructure + cont$bulles + cont$mousse + cont$foisonne + cont$granules
# cont$pain <- cont$pain + cont$paindemie
# cont$plastique <- cont$plastique + cont$caoutchou + cont$latex
# cont$fluide <- cont$semil + cont$fluide
# cont$huile <- cont$huile + cont$gras
# cont$creme <- cont$creme + cont$cremefraiche + cont$cremeux
# cont$fruit <- cont$fruit + cont$fruite + cont$multifruit
# cont$moisi <- cont$moisi + cont$moisissure
# cont$neutre <- cont$nature + cont$neutre+ cont$rien
# cont$biere <- cont$biere + cont$brassage
# cont$chou <- cont$chou + cont$choux
# cont$pois <- cont$pois + cont$poischiche
# cont$lactique <- cont$lactique + cont$lait
# cont$jussoja <- cont$jus + cont$jussoja
# cont$beurre <- cont$beurre + cont$diacetyl +cont$margarine + cont$acetoine
# cont$fort <- cont$odeurforte + cont$fort
# cont$noix <- cont$fruitacoque + cont$noix + cont$noisette
# cont$cacahuete <- cont$cacahuete + cont$curly
# cont$bon <- cont$bon + cont$ok + cont$agreable + cont$gourmand
# cont$brule <- cont$brule + cont$crame
# cont$croquette <- cont$croquette + cont$viande
# cont$soufre <- cont$oeufpourri + cont$soufre
# cont$vin <- cont$vin + cont$mout
# cont$beany <- cont$graine + cont$beany
# cont$fruitpasmur <- cont$fruitpasmur + cont$mangueverte
# cont$foin <- cont$foin + cont$beuh
# cont$filant <- cont$filant + cont$polysaccharides + cont$textureparticuliere
# cont$synerese <- cont$synerese + cont$`synerese-` + cont$dephase + cont$phases
# cont$grille <- cont$grille + cont$toaste + cont$torrefie
# cont$bizarre <- cont$bizarre + cont$odeurbizarre
# cont$pate <- cont$pate + cont$pateapain
# cont$cracotte <- cont$cracotte + cont$gateauapero + cont$acetaldehyde + cont$orge + cont$oxyde + cont$souple + cont$gateau + cont$rizsouffle
# cont$piquant <- cont$piquant + cont$epice
# cont <- select(cont,-lila,-beuh,-odeurbizarre, -viande,-acetoine,-gateauapero,-caille,-torrefie,-toaste,-gourmand,-acetaldehyde,-souple,-orge,-oxyde,-oeufpourri,-mangueverte,-gateau,-rizsouffle, -mout,-rien,-graine,-gruyere,-agrume,-acetique,-jus,-na,-crame,-tapisserie,-vinaigrette,-fruitacoque,-violette,-margarine,-noisette,-farine,-ble,-lait,-muguet,-brassage,-moisissure, -vegetal,-herbe,-herbec,-haricot,-citrique,-gazon,-ferme,-gel,-`synerese+`,-papier,-beaugel,-tresdoux,-fleurcoupe, -fleur,-emmental,-geranium)
# cont <- select(cont,-cremeux,-epice,-curly, -pateapain,-agreable,-textureparticuliere,-polysaccharides,-ok,-phases,-`synerese-` ,-diacetyl,-vieux,-odeurforte,-cremefraiche,-lys,-granules,-jacinthe,-rose,-poischiche,-foisonne,-avoine,-choux,-malte,-solide, -vinaigre,-nature, -multifruit,-fruite, -vinaigrebalsamique, -vinaigredevin,-mousse,-bulles, -semil, -levain,-gras, -flocule,-muesli,-paindemie,-latex,-caoutchou,-destructure)
# as.data.frame(cont)   #,-dephase supprimé au dessus
# FRO <- FRO[order(FRO$RefCIRM),]   ###Pour trier les mots

###Creation Excel###
# write.table(cont1, "StatsCIRMsuptable1.xls")
# as.data.frame(A)
# #####
# # A<- descfreq(res$cont.table,proba=0.18)
# #####???Choisir si filter ou non : 
# row.names(cont1) <- cont1$RefCIRM  ###Afin de savoir qui sont chaque individue sur l'analyse des Correspondances
# row.names(cont2) <- c(cont2$Espece)   ###A modifier selon group Genre/Espece/Tcroissance
# cont1<- cont1[,-c((length(cont1)-6):length(cont1))]    ###Afin d'enlever les colonnes non quanti (RefCIRM, Especes...)
# for(i in 1:length(cont1)){
#   if(sum(cont1[,i])==0){
#     col <- cbind(col,i)
#   }
# }
# cont1 <- cont1[,-col]
# for(i in 1:length(cont2)){
#   if(sum(cont2[,i])==0){
#     col <- cbind(col,i)
#   }
# }
# cont2 <- cont2[,-col]
# 
# A <- descfreq(cont1, proba=0.1)
# B <- descfreq(cont2, proba=0.1)
# 
# # CAmes30A<- descfreq(contmes30, proba=0.18)
# # res.ca = CA(contmes30,graph=F, ncp = 8)
# write.table(A, "descfreq1019D.txt", sep="\t",row.names=F,na="")
# #ellipse#################
# res.ca = CA(cont1,graph=F, ncp = 6)
# plot.CA(res.ca, cex=0.8, axes = c(3,4), autoLab = "yes")
# 
# plot.CA(res.ca,cex=1.5,selectCol = "contrib 50", selectRow = "contrib 30",autoLab = "yes")
# plot.CA(res.ca, cex=1.5, axes = c(3,4),selectCol = "contrib 50", selectRow = "contrib 20", autoLab = "yes")
# plot.CA(res.ca, cex=1.5, axes = c(5,6),selectCol = "contrib 50", selectRow = "contrib 20", autoLab = "yes")
# plot.CA(res.ca, cex=0.8, axes = c(7,8),selectCol = "contrib 50", selectRow = "contrib 20", autoLab = "yes")
# 
# ###
# ellipseCA(res.ca,ellipse="row")
# ?ellipseCA
# ellipseCA(res.ca,ellipse="row",invisible="col")
# ###with colors
# color<-c("green", "orange", rep("black",2), "red")
# ellipseCA(res.ca,ellipse="row",invisible="col",col.row=color)

# #######Représentation des pHs#######
# select(fr,RefCIRM,Espece,`pH 10h`,`pHa T1`,`pHa T2`,`pHb T1`,`pHb T2`) %>%
#   gather(key = mpH, val = pH,-Espece,-RefCIRM)%>%
#   ggplot(aes(y=`pH`, x=RefCIRM)) +
#   geom_boxplot()+
#   theme_tufte()+
#   theme(axis.text.x = element_text(hjust = 1,vjust = 0.5, angle = 90))







###################################
###############TUKEY###############
###################################
mypara1 <- names(meltFRO)[11:16]
mymodel1 <- lapply(mypara1,function(k) {  #D?finition de la fonction test lm
  lm(eval(substitute(j ~ Produit,list(j=as.name(k)))),data = meltFRO[,c(1:5,11:16)] )  #Modifier data selon GC Ou Suc
}
)
tabgather1 <- meltFRO[,c(1:5,11:16)] %>%        #gather des Notes values
  gather(key = parametres, val = value, -`Response ID`,-`Date submitted`, -commentaire, -variable,-Produit)
tabgather1$value <- as.numeric(tabgather1$value)
tabgather1 <- filter(tabgather1, !is.na(value))
####Premier test uniquement sur note globale car certains produits avec une unique note dans une des catégories"
mypara1 <- mypara1[6]
tabgather1 <- filter(tabgather1,parametres =="Noteglobale")
unique(tabgather1$Produit)
####???Début tukey#### Première boucle print numéro de la variable de mypara1 traitée
variable_list1 <- list()
aov1 <- list()
tukey1 <- list()
tukey_groups1 <- list()
for(i in 1:length(mypara1)) { 
  variable_list1[[i]] <- tabgather1[tabgather1$parametres == mypara1[i],] 
  print(i)
  for (i in 1:length(variable_list1)) {
    aov1[[i]] <- (aov(value~Produit,data=variable_list1[[i]]))
    for (i in 1:length(aov1)) {
      tukey1[[i]] <- HSD.test(aov1[[i]],trt="Produit",group=FALSE, alpha = 0.1)
      tukey_groups1[[i]] <- HSD.test(aov1[[i]], trt = "Produit", alpha = 0.1)
    }
  }
}
####données statistique de l'anova
stats_data1 <- list()
for (i in 1:length(variable_list1)) {
  stats_data1[[i]] <- summarySE(variable_list1[[i]], measurevar = "value", groupvars = c("Produit"),na.rm = TRUE)
}
comp1 <- NA
for( i in 1 : length(mypara1)){
  comp1.i <- as.data.frame(tukey1[[i]]$comparison)
  comp1 <- rbind(comp1,comp1.i)
}
names(variable_list1) <- mypara1
names(tukey1) <- mypara1
names(tukey_groups1) <- mypara1
names(aov1) <- mypara1
names(stats_data1) <- mypara1

# pdf("residu.pdf")
my.pvalues1 <- NULL
par(mfrow=c(2,2))
for(i in 1:length(aov1)){
  my.pvalues1 <- rbind.data.frame(my.pvalues1, data.frame(para = names(aov1)[i], p.value = summary(aov1[[i]])[[1]][["Pr(>F)"]][1]))
  plot(aov1[[i]])
  title(names(aov1)[i])
}
# dev.off()
# write.table(my.pvalues1, "pval1_comp.csv", row.names=T, sep="\t",dec=".", na="")
# write.table(df_tukey_letters1, "tukeyiso.csv", row.names=T, sep="\t",dec=".", na="")

##Tukey test with letter groups
##Loop to extract data groups, mypara, and x & y positions for plotting
tukey_letters1<-list() 
i<-1
j<-1
for (i in 1:length(tukey_groups1)) {     
  k <- 1
  tukey_groups1[[i]]$groups$trt <- rownames(tukey_groups1[[i]]$groups)
  levels(tukey_groups1[[i]]$groups$trt) <- rownames(tukey_groups1[[i]]$means)     #tukey_groups[[i]]$groups$trt
  tukey_letters1[[i]] <- tukey_groups1[[i]]$groups
  stats_data1[[i]]$ymax <- stats_data1[[i]]$value+stats_data1[[i]]$se
  stats_data1[[i]]$ymin <- stats_data1[[i]]$value-stats_data1[[i]]$se
  xpos.t1 <- numeric()
  for(j in 1:nrow(tukey_letters1[[i]])) {
    xpos.t1[k] <- as.integer(rownames(stats_data1[[i]][stats_data1[[i]]$Produit==(as.vector(row.names(tukey_letters1[[i]])[j])),]))
    k<-k+1
    tukey_letters1[[i]]$Especes[j] <- filter(variable_list1[[1]], variable == tukey_letters1[[i]]$trt[j])[1,3]
    tukey_letters1[[i]]$Names[j] <- filter(variable_list1[[1]], variable == tukey_letters1[[i]]$trt[j])[1,2]
  }
  ymax1<-numeric()
  ymin1<-numeric()
  for(m in 1:length(xpos.t1))
  {
    ymax1[m]<-stats_data1[[i]][xpos.t1[m],"ymax"]
    ymin1[m]<-stats_data1[[i]][xpos.t1[m],"ymin"]
  }
  tukey_letters1[[i]]<-cbind(tukey_letters1[[i]],data.frame(xpos.t1,ymax1,ymin1))
}
##Naming list levels
names(tukey_letters1) <- mypara1
tukey_letters1
##We add new column in  tukey_letters that correspond with the name of the level list. This is very important for "facets"  
## function in ggplot.
for (i in 1:length(tukey_letters1)){  
  group<-numeric()
  group<-rep(names(tukey_letters1[i]),length.out=nrow(tukey_letters1[[i]]))
  tukey_letters1[[i]]<-cbind(tukey_letters1[[i]],group)
}
#We eliminate the levels of tukey_letters so we can set the postion and labels of each 
# tukey letters group automatically in our ggplot.   
df_tukey_letters1<-data.frame()
for(i in 1:length(tukey_letters1)){ 
  df_tukey_letters1<-rbind(df_tukey_letters1,tukey_letters1[[i]])
}





for(i in mypara1){     
  print(filter(df_tukey_letters1,group == i) %>%
          ggplot(aes(x= trt, y= value, ymax = ymax1, ymin=ymin1, fill= trt)) +
          geom_bar(stat = "identity")+
          # fillScale +
          #facet_wrap(~group, scales = "free", ncol = 5) +
          geom_errorbar() +
          theme_tufte() +
          # scale_y_continuous(limits = c(0, 0.1)) +
          theme (axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.3,size=15),title = element_text(size=15),axis.text.y = element_text(size=15),legend.text = element_text(size=15,face = "italic"))+   #,legend.position = "bottom" 
          geom_text(data =filter(df_tukey_letters1,group == i), aes(x=trt, y=ymax1,label=groups),vjust=0.2,hjust=-0.5,angle=90, size = 4,position=position_dodge(.5))+
          #geom_text(aes(label=round(value, digits=2)), vjust=3, color="black",position = position_dodge(0.9), size=2) +
          labs(x = "variable", y = "Concentration (en mg/L)", title =i))
}
