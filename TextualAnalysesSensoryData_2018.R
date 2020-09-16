#Script qui m'a surtout servi à assembler les descripteurs et à les transformer en fréquence
#via la fonction textual(). Voir le tableau TestConsortiums_062018.xlsx pour connaître le contenu
#des colonnes. Jusqu'au 31/05 le travail réalisé s'est contenté d'assemble tous les descripteurs afin
# de réaliser une analyse globale sur les descripteurs.
#Il faudra maintenant étudier: 
# L'effet souche 
# L'effet espèce
# L'effet pH
# L'effet Température
# 
# 

rm(list = ls()) #vider l'environnement
#chargement des packages
library(readxl)
library(writexl)
library(FactoMineR)
library(factoextra)
library(tidyverse) 
library(ggthemes) ##permet d'utiliser %>%
?textual
?descfreq
?

setwd("C:/Users/a0148148/Documents/Dossiers Olivier HARLE/R") #retrouve les fichiers - definition de l'espace de travail


#____________ lecture du fichier resultat et verifications
fr <-read_xlsx("~/Dossiers Olivier HARLE/INRA/Tabtravail/TestConsortiums_092018.xlsx", sheet=4, col_names = T) 
names(fr) <- fr[1,] 
fr <- fr[-1,]
fr$Especes1 <- factor(fr$Especes1, levels = c("acidophilus","casei","coryniformis","cvurvatus","delbrueckii","johnsonii","lactis","lactiscremoris","paracasei","pentosus","plantarum","t","thermophilus","paraplantarum"))
fr$Especes2 <- factor(fr$Especes2, levels = c("acidophilus","casei","coryniformis","cvurvatus","delbrueckii","johnsonii","lactis","lactiscremoris","paracasei","pentosus","plantarum","t","thermophilus","paraplantarum"))
fr$seance <- as.factor(fr$seance)
fr$Essai <- as.factor(fr$Essai)
fr$pH10h <- as.numeric(fr$pH10h)
fr$DegustINRA<- as.factor(fr$DegustINRA)
fr$DegustTrib<- as.factor(fr$DegustTrib)

FRO <- fr[,c(7:8,16:20,grep(pattern = "Autre",names(fr)))]
FRO <- filter(FRO,seance == 13|seance == 14|seance == 15|seance == 16|seance == 17|seance == 18|seance == 19)

FRO$seance
# for(i in 22:length(names(fr))){
#   fr[,i] <- as.numeric(unlist(fr[,i]))
# }
####Treatment OVER#####
names(FRO)
FRO$Ctout <- paste(FRO$AutreSD,FRO$AutreHF,FRO$AutreSP,FRO$AutreOI,FRO$AutreAT,FRO$AutreEG,FRO$AutreJN,FRO$AutreCT,FRO$AutreME,FRO$AutreAP,FRO$AutreOT,FRO$AutreTG,FRO$AutreCC,FRO$AutreCC,FRO$AutreML, sep=";")
FRO$Essai <- droplevels(FRO$Essai)
# FRO <- filter(FRO, nchar(FRO$Ctout) != 44)
summary(FRO)

##############Analyses par espèces##############
###GLOBAL###
res <-  textual(as.data.frame(FRO), num.text=22,contingence.by=1,sep.word =";")
lname <- rownames(res$cont.table)
cont <- as.data.frame(res$cont.table)
# rownames(cont) <- lname
# cont[gacide,]
# l <- !is.na(cont$splat)
# cont <- cont[l,]
# cont <- filter(cont, !is.na(splat))
# write.table(cont, "dcont.txt", sep="\t",row.names=T,na="")
# cont[38,]
# cont1 <- descfreq(cont,proba=0.05)
# rownames(cont)
# cont1[1]
# cont1$`2103-t-JC-43`
# cont$sepais

descriptT <- row.names(res$nb.words)
for(d in descriptT){     ###Permet d'ajouter toutes les colonnes des descripteurs qui ne sont pas dans les tables --> permet la fusion des descripteurs sans soucis
  if(d %in% names(cont)){
  } else {
    col <- data.frame(d=rep(0,nrow(cont)))
    names(col) <- d
    cont <- cbind(cont,col)
  }
}
#####Fusion des colonnes#####
#Global#
cont$srustique <- cont$sagrume + cont$srustique
cont$scereal <- cont$scereal + cont$smuesli+ cont$savoine
cont$svert <- cont$svert + cont$svegetal
cont$ocitron <- cont$ocitron + cont$ocitrique
cont$ofromage <- cont$ofromage + cont$oemmental
cont$scitron <- cont$scitron + cont$scitrique
cont$sfloral <- cont$sfloral +cont$sglycine
cont$saigre <- cont$saigre + cont$svinaigre
cont$shuile <- cont$shuile + cont$sgras
cont$screme <- cont$screme + cont$scremeux
cont$schou <- cont$schou + cont$schoux
cont$snoix <- cont$snoix + cont$snoisette
cont$smoisi <- cont$spourri + cont$smoisi
cont <- select(cont,-sepice,-oemmental,-ocitrique,-snoisette,-schoux, -sagrume,-smuesli,-savoine,-svegetal,-scitrique,-sglycine,-svinaigre,-sgras,-scremeux)


###Creation Excel###
write.table(cont, "cont1218.xls")
#levels(fro$V1) = c("C","Si","Sp")
#ou
#levels(fro$V1) = c("C1","C2","C3","Si1","Si2","Si3","Sp1","Sp2","Sp3")    
#Ou ajout de cont à la suite de fr

Senso <- read_xlsx("~/Dossiers Olivier HARLE/R/cont1218.xlsx", sheet=1, col_names = T) 
Senso <- as.data.frame(Senso)
row.names(Senso) <- as.data.frame(Senso)[,1]
summary(Senso)
res <- MFA (Senso[,c(1:7)],group=c(1,1,5), name.group = c("REFCIRM","pH","Sens"),type = c("n","s","s"),ncp=6, graph = TRUE)
res <- MFA (Senso,group=c(1,1,5,147), name.group = c("REFCIRM","pH","Sens","Desc"),type = c("n","s","s","f"),ncp=6,num.group.sup = 1, graph = TRUE)




###Analyse Var
quanti.var <- get_mfa_var(res)
head(quanti.var$cos2)
fviz_mfa_var (res,labelsize=5,geom = c("arrow","text"))  +theme_tufte()
fviz_mfa_var (res,cos2 =0.3,axes=c(3,4),labelsize=5,select.var = list(name = NULL, cos2 = 0.1, contrib = NULL),geom = c("arrow","text"))  +theme_tufte()
fviz_mfa_var (res,cos2 =0.3,axes=c(5,6),labelsize=5,select.var = list(name = NULL, cos2 = 0.1, contrib = NULL),geom = c("arrow","text"))  +theme_tufte()

fviz_contrib (res, choice = "quanti.var", axes = 1,top=20)

###Analyse Ind
ind <-get_mfa_ind (res)
head(ind$cos2)

fviz_mfa_ind (res,repel=T, invisible =("ind")) +theme_tufte()
fviz_mfa_ind (res,repel=T, invisible =("ind"),axes=c(5,6)) +theme_tufte()
fviz_mfa_ind (res,invisible =("quali.var"),geom = c("text")) +theme_tufte()
fviz_mfa_ind (res, habillage=2,invisible =("quali.var"),axes=c(3,4),geom = c("point")) +theme_tufte()
fviz_mfa_ind (res, habillage=2,invisible =("quali.var"),axes=c(5,6),geom = c("point")) +theme_tufte()














names(Senso)
#####
A<- descfreq(res$cont.table,proba=0.18)
A<- descfreq(cont, proba=0.18)
A<- descfreq(contI, proba=0.18)
A<- descfreq(contT, proba=0.18)
CAmes30A<- descfreq(contmes30, proba=0.18)
res.ca = CA(contmes30,graph=F, ncp = 8)
CAmes37A<- descfreq(contmes37, proba=0.18)
res.ca = CA(contmes37,graph=F, ncp = 8)
CAthe43A<- descfreq(contthe43, proba=0.18)
res.ca = CA(contthe43,graph=F, ncp = 8)
CAthe37A<- descfreq(contthe37, proba=0.18)
res.ca = CA(contthe37,graph=F, ncp = 8)

CA(cont1,graph=F, ncp = 8)

res=list()
  res[[1]]=A$C1;
  res[[3]]=A$C2;
  res[[2]]=A$C3;
write.table(res, "descfreq1.txt", sep="\t",row.names=F,na="")


#ellipse#################
res.ca = CA(cont,graph=F, ncp = 6)
plot.CA(res.ca, cex=0.8)
plot.CA(res.ca, cex=0.8, axes = c(3,4), autoLab = "yes")

plot.CA(res.ca,cex=0.8,selectCol = "contrib 50", selectRow = "contrib 30",autoLab = "yes")
plot.CA(res.ca, cex=0.8, axes = c(3,4),selectCol = "contrib 50", selectRow = "contrib 20", autoLab = "yes")
plot.CA(res.ca, cex=0.8, axes = c(5,6),selectCol = "contrib 50", selectRow = "contrib 20", autoLab = "yes")
plot.CA(res.ca, cex=0.8, axes = c(7,8),selectCol = "contrib 50", selectRow = "contrib 20", autoLab = "yes")


###
ellipseCA(res.ca,ellipse="row")
?ellipseCA
ellipseCA(res.ca,ellipse="row",invisible="col")
  
###with colors
color<-c("green", "orange", rep("black",2), "red")
ellipseCA(res.ca,ellipse="row",invisible="col",col.row=color)




#######Représentation des pHs#######
select(fr,RefCIRM,Espece,`pH 10h`,`pHa T1`,`pHa T2`,`pHb T1`,`pHb T2`) %>%
  gather(key = mpH, val = pH,-Espece,-RefCIRM)%>%
  ggplot(aes(y=`pH`, x=RefCIRM)) +
  geom_boxplot()+
  theme_tufte()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 0.5, angle = 90))

select(fr,RefCIRM,Espece,`pHa`,`pHb`,`pHa 37`,`pHb 37`,Tcroissance) %>%
  gather(key = mpHT, val = `30`,-Espece,-Tcroissance,-RefCIRM,-`pHa 37`,-`pHb 37`)%>%
  gather(key = mpH37, val = `37`,-Espece,-Tcroissance,-RefCIRM,-mpHT,-`30`) %>%
  gather(key =Temperature, val =pH,-Tcroissance,-Espece,-RefCIRM,-mpH37,-mpHT)%>%
  filter(Tcroissance==30)  %>%
  ggplot(aes(y=`pH`, x=Temperature,color=Espece)) +
  geom_boxplot()+
  theme_base() +
  scale_y_continuous(limits=c(4.1,7))+
  theme(axis.text.x = element_text(hjust = 1,vjust = 0.5, angle = 90))+
  facet_wrap(~RefCIRM, ncol=9)

select(fr,RefCIRM,Espece,`pHa`,`pHb`,`pHa 37`,`pHb 37`,Tcroissance) %>%
  gather(key = mpHT, val = `43`,-Espece,-Tcroissance,-RefCIRM,-`pHa 37`,-`pHb 37`)%>%
  gather(key = mpH37, val = `37`,-Espece,-Tcroissance,-RefCIRM,-mpHT,-`43`) %>%
  gather(key =Temperature, val =pH,-Tcroissance,-Espece,-RefCIRM,-mpH37,-mpHT)%>%
  filter(Tcroissance==43)  %>%
  ggplot(aes(y=`pH`, x=Temperature,color=Espece)) +
  geom_boxplot()+
  theme_base() +
  scale_y_continuous(limits=c(4.1,7))+
  theme(axis.text.x = element_text(hjust = 1,vjust = 0.5, angle = 90))+
  facet_wrap(~RefCIRM, ncol=9)

select(fr,RefCIRM,Espece,`pHa`,`pHb`,`pHa 37`,`pHb 37`,Tcroissance) %>%
  gather(key = mpHT, val = `43`,-Espece,-Tcroissance,-RefCIRM,-`pHa 37`,-`pHb 37`)%>%
  gather(key = mpH37, val = `37`,-Espece,-Tcroissance,-RefCIRM,-mpHT,-`43`) %>%
  gather(key =Temperature, val =pH,-Tcroissance,-Espece,-RefCIRM,-mpH37,-mpHT)%>%
  filter(Tcroissance==37)  %>%
  ggplot(aes(y=`pH`, x=Temperature,color=Espece)) +
  geom_boxplot()+
  theme_base() +
  scale_y_continuous(limits=c(4.1,7))+
  theme(axis.text.x = element_text(hjust = 1,vjust = 0.5, angle = 90))+
  facet_wrap(~RefCIRM, ncol=9)

# select(fr,RefCIRM,Espece,`pH 10h`,`pHa T1`,`pHa T2`,`pHb T1`,`pHb T2`)%>%
#   filter()
#  
# select(fr,RefCIRM,Espece,`pHa 37`,`pHb 37`) %>%
#   gather(key = mpH, val = pH,-Espece,-RefCIRM)%>%
#   ggplot(aes(y=`pH`, x=RefCIRM,color=Espece)) +
#   geom_boxplot()+
#   theme_tufte()+
#   theme(axis.text.x = element_text(hjust = 1,vjust = 0.5, angle = 90))

#######Représentation des notes#######
select(fr,RefCIRM,Espece,Tcroissance,NoteT1,NoteT2,NoteT3,NoteT4,NoteT5,NoteT6,NoteT7,Note371,Note372,Note373,Note374,Note375,Note376,Note377) %>%
  gather(key =numnotes37, val =37,-Tcroissance,-Espece,-RefCIRM,-c(NoteT1,NoteT2,NoteT3,NoteT4,NoteT5,NoteT6,NoteT7))%>%
  gather(key =numnotesT, val =43,-Tcroissance,-Espece,-RefCIRM,-`37`,-numnotes37)%>%
  gather(key =Tnotes, val =Notes,-Tcroissance,-Espece,-RefCIRM,-numnotes37,-numnotesT)%>%
  filter(Tcroissance == 43) %>%
  ggplot(aes(y=`Notes`, x=Tnotes,color=Espece)) +
  geom_boxplot()+
  theme_base()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 0.5, angle = 90))+
    facet_wrap(~RefCIRM, ncol=9)

select(fr,RefCIRM,Espece,Tcroissance,NoteT1,NoteT2,NoteT3,NoteT4,NoteT5,NoteT6,NoteT7,Note371,Note372,Note373,Note374,Note375,Note376,Note377) %>%
  gather(key =numnotes37, val =37,-Tcroissance,-Espece,-RefCIRM,-c(NoteT1,NoteT2,NoteT3,NoteT4,NoteT5,NoteT6,NoteT7))%>%
  gather(key =numnotesT, val =43,-Tcroissance,-Espece,-RefCIRM,-`37`,-numnotes37)%>%
  gather(key =Tnotes, val =Notes,-Tcroissance,-Espece,-RefCIRM,-numnotes37,-numnotesT)%>%
  filter(Tcroissance == 30|37) %>%
  ggplot(aes(y=`Notes`, x=Tnotes,color=Espece)) +
  geom_boxplot()+
  theme_base()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 0.5, angle = 90))+
  facet_wrap(~RefCIRM, ncol=9)

select(fr,RefCIRM,Espece,Tcroissance,NoteT1,NoteT2,NoteT3,NoteT4,NoteT5,NoteT6,NoteT7,Note371,Note372,Note373,Note374,Note375,Note376,Note377) %>%
  gather(key =numnotes37, val =37,-Tcroissance,-Espece,-RefCIRM,-c(NoteT1,NoteT2,NoteT3,NoteT4,NoteT5,NoteT6,NoteT7))%>%
  gather(key =numnotesT, val =43,-Tcroissance,-Espece,-RefCIRM,-`37`,-numnotes37)%>%
  gather(key =Tnotes, val =Notes,-Tcroissance,-Espece,-RefCIRM,-numnotes37,-numnotesT)%>%
  filter(Tcroissance == 37) %>%
  ggplot(aes(y=`Notes`, x=Tnotes,color=Espece)) +
  geom_boxplot()+
  theme_base()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 0.5, angle = 90))+
  facet_wrap(~RefCIRM, ncol=9)

# select(fr,RefCIRM,Espece,Tcroissance,NoteT1,NoteT2,NoteT3,NoteT4,NoteT5,NoteT6,NoteT7)%>%
#   filter(Tcroissance==30)  %>%
#   gather(key =numnote, val =Notes,-Espece,-RefCIRM,-Tcroissance)%>%
#   ggplot(aes(y=`Notes`, x=RefCIRM,color=Espece)) +
#   geom_boxplot()+
#   theme_tufte()+
#   theme(axis.text.x = element_text(hjust = 1,vjust = 0.5, angle = 90))
# select(fr,RefCIRM,Espece,Note371,Note372,Note373,Note374,Note375,Note376,Note377) %>%
#   gather(key =numnote, val =Notes,-Espece,-RefCIRM)%>%
#   ggplot(aes(y=`Notes`, x=RefCIRM,color=Espece)) +
#   geom_boxplot()+
#   theme_tufte()+
#   theme(axis.text.x = element_text(hjust = 1,vjust = 0.5, angle = 90))
# select(fr,RefCIRM,Espece,Tcroissance,NoteT1,NoteT2,NoteT3,NoteT4,NoteT5,NoteT6,NoteT7)%>%
#   filter(Tcroissance==43)  %>%
#   gather(key =numnote, val =Notes,-Espece,-RefCIRM,-Tcroissance)%>%
#   ggplot(aes(y=`Notes`, x=RefCIRM,color=Espece)) +
#   geom_boxplot()+
#   theme_tufte()+
#   theme(axis.text.x = element_text(hjust = 1,vjust = 0.5, angle = 90))
