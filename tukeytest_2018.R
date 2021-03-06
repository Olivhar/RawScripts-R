rm(list = ls()) #vider l'environnement
#chargement des packages
# library(gplots)  ###Pour heatmaps.2
# library(pheatmap)
library(readxl)
library(writexl)
library(FactoMineR)
library(factoextra)
library(tidyverse) 
library(ggthemes)
library(agricolae) #pour t.test par paires --> HSD.test()
library(Rmisc)     #pour la fonction summarySE()
library(reshape)   #pour la fonction melt()
setwd("C:/Users/a0148148/Documents/Dossiers Olivier HARLE/R") #retrouve les fichiers - definition de l'espace de travail

##############################################
##############################################
##################TUKEY TEST##################
##############################################
##############################################

# Importation du jeu de donn�es
Suc <-read_xlsx("~/Dossiers Olivier HARLE/INRA/.xlsx", sheet=3, col_names = T) 

# D�finition des parametres a comparer (selon le nom des colonnes)
mypara1 <- names(Suc[c(1:50)])  ###


# TRAITEMENT DES DONNEES (SI NECESSAIRE)
#######CALCUL DES LOG10 (Avec si value=0 alors value=0.01 (car sinon bug: log(0)=-Inf)#####
for(i in mypara1[1:37]){      ##Pr�ciser de 1:37 pour les donn�es GC
  for(j in 1: length(Suc[,i])){
    if(Suc[i][j,]==0){
    Suc[i][j,] <- 0.01
    }
  }
}
Suc[c(1:37)] <- data.frame(lapply(Suc[c(1:37)], function(x) log10(x)))


###Creation d'un second tableau avec une colonne valeur et un colonne parametres (si inverse, faire spread au lieu de gather pour obtenir ma table Suc)
# Suc$Name <- Suc$Espece ####Faire modif si anova par especes ou d�finir les groupes si autres segment
tabgather1 <- Suc[,-c(53,55:56,59)] %>%        #POUR LES SUCRES
  gather(key = parametres, val = value, -RefCIRM,-Name, -Espece)

mymodel1 <- lapply(mypara1,function(k) {  #D?finition de la fonction test lm
  lm(eval(substitute(j ~ RefCIRM,list(j=as.name(k)))),data = Suc )  #Modifier data selon GC Ou Suc et modifier RefCIRM selon la nom de la colonne des individus a comparer
}
)
tabgather1 <- filter(tabgather1, !is.na(value))
tabgather1$value <- as.numeric(tabgather1$value)


####???D�but tukey#### Premi�re boucle print num�ro de la variable de mypara1 trait�e
variable_list1 <- list()
aov1 <- list()
tukey1 <- list()
tukey_groups1 <- list()
for(i in 1:length(mypara1)) { 
  variable_list1[[i]] <- tabgather1[tabgather1$parametres == mypara1[i],] 
  print(i)
  for (i in 1:length(variable_list1)) {
    aov1[[i]] <- (aov(value~RefCIRM,data=variable_list1[[i]]))
    for (i in 1:length(aov1)) {
      tukey1[[i]] <- HSD.test(aov1[[i]],trt="RefCIRM",group=FALSE, alpha = 0.1)
      tukey_groups1[[i]] <- HSD.test(aov1[[i]], trt = "RefCIRM", alpha = 0.1)
    }
  }
}

####donn�es statistique de l'anova
stats_data1 <- list()
for (i in 1:length(variable_list1)) {
  stats_data1[[i]] <- summarySE(variable_list1[[i]], measurevar = "value", groupvars = c("RefCIRM"),na.rm = TRUE)
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
write.table(my.pvalues1, "pval_comp.csv", row.names=T, sep="\t",dec=".", na="")

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
    xpos.t1[k] <- as.integer(rownames(stats_data1[[i]][stats_data1[[i]]$RefCIRM==(as.vector(row.names(tukey_letters1[[i]])[j])),]))
    k<-k+1
    tukey_letters1[[i]]$Especes[j] <- filter(variable_list1[[1]], RefCIRM == tukey_letters1[[i]]$trt[j])[1,3]
    tukey_letters1[[i]]$Names[j] <- filter(variable_list1[[1]], RefCIRM == tukey_letters1[[i]]$trt[j])[1,2]
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

#Create a custom color scale
library(scales)
scales::hue_pal()(14)
myColors <- c("#F8766D", "#E18A00", "#BE9C00", "#8CAB00", "#24B700", "#00BE70", "#00C1AB", "#00BBDA", "#00ACFC", "#8B93FF", "#D575FE", "#F962DD", "#FF65AC")
df_tukey_letters2 <- cbind(df_tukey_letters1[1:9])
df_tukey_letters2$Especes <- factor(c(df_tukey_letters2$Especes), levels = c("a","b","c","d","e", "f", "g","h", "i"))
names(myColors) <- levels(df_tukey_letters2$Especes)
fillScale <- scale_fill_manual(name = "Especes",values = myColors)

# factor(df_tukey_letters2$Especes,levels = c("acidophilus","casei","coryniformis","cvurvatus","delbrueckii","johnsonii","lactis","lactiscremoris","paracasei","pentosus","plantarum","t","thermophilus"))


####BOUCLES POUR TRACER TOUT LES GRAPHIQUES####
for(i in mypara1){     
print(filter(df_tukey_letters2,group == i)[1:53,] %>%
  ggplot(aes(x= trt, y= value, ymax = ymax1, ymin=ymin1, fill= Especes)) +
  geom_bar(stat = "identity")+
  fillScale +
  #facet_wrap(~group, scales = "free", ncol = 5) +
  geom_errorbar() +
  theme_tufte() +
  # scale_y_continuous(limits = c(0, 0.1)) +
  theme (axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.3,size=15),title = element_text(size=15),axis.text.y = element_text(size=15),legend.text = element_text(size=15,face = "italic"))+   #,legend.position = "bottom" 
  geom_text(data =filter(df_tukey_letters2,group == i), aes(x=trt, y=ymax1,label=groups),vjust=0.2,hjust=0.2,angle=90, size = 4,position=position_dodge(.5))+
  #geom_text(aes(label=round(value, digits=2)), vjust=3, color="black",position = position_dodge(0.9), size=2) +
  labs(x = "RefCIRM", y = "Value", title =i))
}


######FIN#####


