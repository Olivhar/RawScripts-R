rm(list = ls()) #vider l'environnement
#chargement des packages
##################################################
library(readxl)
# library(writexl)
library(ggthemes)
#library(gsheet)
# library(grofit)
# library(deSolve)
# library(stringr)
# library(gridExtra)
#library(ggplot2) #Dejà dans tidyverse
#library(segmented)
# library(FactoMineR)
library(agricolae) #pour t.test par paires --> HSD.test()
library(Rmisc) #pour la fonction summarySE() ET FAIT BUGGER group_by() summarize()
#library(factoextra)#Permet de customiser le graph des ACM
library(reshape)   #pour la fonction melt()
library(data.table) ###Pour la fonction grep() et pour la fonction melt() --> package utilisé pour melt par patterns (utilisé setDT si table importé d'excel)
library(tidyverse)
# library(rmarkdown) #Permet de mieux décrire les scripts et de les imprimer en pdf ou word
# library(VennDiagram)
##################################################
setwd("C:/Users/a0148148/Documents/Dossiers Olivier HARLE/R") #retrouve les fichiers - definition de l'espace de travail
#Definition des fonctions
erreur.std <- function(x)
{ x <- na.omit(x)
ET <- sd(x)/sqrt(length(x))
return(ET)
}

##################################################
#Chargement du jeu de donnees data et mise en forme. 
##################################################
consor <- read_xlsx("~/Dossiers Olivier HARLE/INRA/EtudeConsortia/suivi pH_denum/essaicinetiquesoja_050819.xlsx", sheet=7) %>%
  tbl_df()
names(consor)<- consor[1,]
consor <- consor[-1,]
for(i in 4:length(consor)){
  consor[,i]<- as.numeric(unlist(consor[,i]))
}
summary(consor)
consor<-filter(consor,Souche!="T")
################################################
################################################
###################TUKEY #######################
################################################
################################################
#####POUR PASSER EN MOL/L#####
# consorS <- mutate(consorS, sucrosepH4.5=sucrosepH4.5/342.962,raffinosepH4.5=raffinosepH4.5/504.42,verbascosepH4.5=verbascosepH4.5/828.7183,stacchyosepH4.5=stacchyosepH4.5/666.5777,glucosepH4.5=glucosepH4.5/180.1559,fructosepH4.5=fructosepH4.5/180.1559,galactosepH4.5=galactosepH4.5/180.1559)
##############################
mypara1 <- names(consor)[4:length(names(consor))]
mypara1<-mypara1[1:(length(mypara1)-2)]
mymodel1 <- lapply(mypara1,function(k) {  #D?finition de la fonction test lm
  lm(eval(substitute(j ~ Consortia,list(j=as.name(k)))),data = consor )  
}
)
tabgather1 <- consor[,-c(1,25:26)] %>%        
  gather(key = parametres, val = value,-Consortia,-Replicat)
# tabgather1$Consortia <- paste(tabgather1$Consortia,tabgather1$parametres)
####???Début tukey#### Première boucle print numéro de la variable de mypara1 traitée
variable_list1 <- list()
aov1 <- list()
tukey1 <- list()
tukey_groups1 <- list()
for(i in 1:length(mypara1)) { 
  variable_list1[[i]] <- tabgather1[tabgather1$parametres == mypara1[i],] 
  print(i)
  for (i in 1:length(variable_list1)) {
    aov1[[i]] <- (aov(value~Consortia,data=variable_list1[[i]]))
    for (i in 1:length(aov1)) {
      tukey1[[i]] <- HSD.test(aov1[[i]],trt="Consortia",group=FALSE, alpha = 0.05)
      tukey_groups1[[i]] <- HSD.test(aov1[[i]], trt = "Consortia", alpha = 0.05)
    }
  }
}
####données statistique de l'anova
stats_data1 <- list()
for (i in 1:length(variable_list1)) {
  stats_data1[[i]] <- summarySE(variable_list1[[i]], measurevar = "value", groupvars = c("Consortia"),na.rm = TRUE)
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
    xpos.t1[k] <- as.integer(rownames(stats_data1[[i]][stats_data1[[i]]$Consortia==(as.vector(row.names(tukey_letters1[[i]])[j])),]))
    k<-k+1
    tukey_letters1[[i]]$Especes[j] <- filter(variable_list1[[1]], Consortia == tukey_letters1[[i]]$trt[j])[1,3]
    tukey_letters1[[i]]$Names[j] <- filter(variable_list1[[1]], Consortia == tukey_letters1[[i]]$trt[j])[1,2]
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
df_tukey_letters1

# writexl::write_xlsx(df_tukey_letters1,path="tukey_Dnumarticle.xlsx")
#Create a custom color scale
myColors <- c("gold2","gold2","red3","red3","black")
df_tukey_letters1$trt <- factor(df_tukey_letters1$trt, levels = c("B","B(BS)","S","S(BS)","T"))
names(myColors)<- levels(df_tukey_letters1$trt)
fillScale <- scale_fill_manual(name = "trt",values = myColors)

####BOUCLES POUR TRACER TOUT LES GRAPHIQUES####
for(i in mypara1){     
  print(filter(df_tukey_letters1,group == i) %>%
          ggplot(aes(x= trt, y= value, ymax = ymax1, ymin=ymin1, fill=trt)) +
          geom_bar(stat = "identity")+
          fillScale +
          #facet_wrap(~group, scales = "free", ncol = 5) +
          geom_errorbar() +
          theme_tufte() +
          scale_y_continuous(limits = c(0, 0.001)) +
          theme (axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.3,size=15),title = element_text(size=15),axis.text.y = element_text(size=15),legend.text = element_text(size=15,face = "italic"))+   #,legend.position = "bottom" 
          geom_text(data =filter(df_tukey_letters1,group == i), aes(x=trt, y=ymax1,label=groups),vjust=0.2,hjust=0.2,angle=90, size = 4,position=position_dodge(.5))+
          #geom_text(aes(label=round(value, digits=2)), vjust=3, color="black",position = position_dodge(0.9), size=2) +
          labs(x = "Consortia", y = "en mol/L", title =i))
}
# ,,ymax1,ymin1,groups
# df_tukey_letters1$trt<- row.names(df_tukey_letters1)
df_tukey_letters1 <- filter(df_tukey_letters1,trt!="T")


# DATABLE values
spreadTukey <-select(df_tukey_letters1,trt,value,group)%>%
spread(group,value)
row.names(spreadTukey)<-spreadTukey$trt
meltTukey <- melt.data.table(setDT(select(spreadTukey,-trt)), measure = patterns("^t","^pH","^D"),value.name = c("t","pH","D"))
meltTukey <- as_tibble(meltTukey)
meltTukey$variable<- rep(spreadTukey$trt,6)
# DATABLE Ymin
spreadymin<-select(df_tukey_letters1,trt,ymin1,group)%>%
  spread(group,ymin1)
row.names(spreadymin)<-spreadymin$trt
meltymin <- melt.data.table(setDT(select(spreadymin,-trt)), measure = patterns("^t","^pH","^D"),value.name = c("mint","minpH","minD"))
# DATABLE Ymax
spreadymax<-select(df_tukey_letters1,trt,ymax1,group)%>%
  spread(group,ymax1)
row.names(spreadymax)<-spreadymax$trt
meltymax <- melt.data.table(setDT(select(spreadymax,-trt)), measure = patterns("^t","^pH","^D"),value.name = c("maxt","maxpH","maxD"))
# DATABLE group
spreadgroup<-select(df_tukey_letters1,trt,groups,group)%>%
  spread(group,groups)
row.names(spreadgroup)<-spreadgroup$trt
meltgroup <- melt.data.table(setDT(select(spreadgroup,-trt)), measure = patterns("^t","^pH","^D"),value.name = c("groupt","grouppH","groupD"))
# CBIND
meltTukey<- cbind(meltTukey,meltymin[,-1],meltymax[,-1],meltgroup[,-1])

####GGPLOT ACIDIFICATION en fonction du TEMPS
ggplot(meltTukey,aes(x=t,y=pH,color=variable))+
  geom_line() +
  # geom_smooth(se=F,span=0.7)+
  geom_point(size=1)+
  geom_pointrange(aes(ymin=minpH,ymax=maxpH),fatten=1)+
  geom_errorbarh(aes(xmax=mint,xmin =maxt, height = 0))+
  theme_tufte()+
  # scale_y_continuous(limits=c(4,7.5))+
  # scale_x_continuous(limits=c(0,30))+
  scale_color_manual(values = c("gold2","darkorange2","red2","darkorange2"))
meltTukey$variable <- rep(c("B","B(BS)","S","S(BS)"),7) 
####GGPLOT NUM en fonction du pH
ggplot(meltTukey,aes(x=pH,y=D,color=variable,linetype=variable))+
  geom_line() +
  # geom_smooth(se=F,span=0.7)+
  geom_point(size=1)+
  geom_pointrange(aes(ymin=minD,ymax=maxD),fatten=1)+
  geom_errorbarh(aes(xmax=minpH,xmin =maxpH, height = 0))+
  theme_tufte()+
  scale_y_continuous(limits=c(10^4,10^9),trans='log10')+
  # scale_x_continuous(limits=c(0,30))+
  scale_color_manual(values = c("gold2","gold2","red2","red2"))+
  scale_linetype_manual(values = c(1,2,1,2))+
  geom_text(data =meltTukey, aes(x=pH, y=maxD,label=groupD),vjust=0.2,hjust=0.2,angle=0, size = 4,position=position_dodge(.5))+
  scale_x_reverse()
####GGPLOT NUM en fonction du t
ggplot(meltTukey,aes(x=t,y=D,color=variable,linetype=variable))+
  geom_line() +
  # geom_smooth(se=F,span=0.7)+
  geom_point(size=1)+
  geom_pointrange(aes(ymin=minD,ymax=maxD),fatten=1)+
  geom_errorbarh(aes(xmax=mint,xmin =maxt, height = 0))+
  theme_tufte()+
  scale_y_continuous(limits=c(10^4,10^9),trans='log10')+
  scale_x_continuous(limits=c(0,30))+
  scale_color_manual(values = c("gold2","gold2","red2","red2"))+
  scale_linetype_manual(values = c(1,2,1,2))+
  geom_text(data =meltTukey, aes(x=pH, y=maxD,label=groupD),vjust=0.2,hjust=0.2,angle=0, size = 4,position=position_dodge(.5))











##################################################
###CALCUL des moyennes et ERREURS STANDARD########
##################################################
m.consor<-consor%>%
  group_by(Souche,Consortia)%>%
  dplyr::summarise(m.t0 = mean (t0),
            m.t1 = mean (t1,na.rm=T),
            m.t2 = mean (t2,na.rm=T),
            m.t3 = mean (t3,na.rm=T),
            m.t4 = mean (t4,na.rm=T),
            m.t5 = mean (t5,na.rm=T),
            m.t6 = mean (t6,na.rm=T),
            m.pH0 = mean (pH0,na.rm=T),
            m.pH1 = mean (pH1,na.rm=T),
            m.pH2 = mean (pH2,na.rm=T),
            m.pH3 = mean (pH3,na.rm=T),
            m.pH4 = mean (pH4,na.rm=T),
            m.pH5 = mean (pH5,na.rm=T),
            m.pH6 = mean (pH6,na.rm=T),
            m.Da0 = mean(c(Da0,Db0),na.rm=T),
            m.D1 = mean(c(Da1,Db1),na.rm=T),
            m.D2 = mean(c(Da2,Db2),na.rm=T),
            m.D3 = mean(c(Da3,Db3),na.rm=T),
            m.D4 = mean(c(Da4,Db4),na.rm=T),
            m.D5 = mean(c(Da5,Db5),na.rm=T),
            m.D6 = mean(c(Da6,Db6),na.rm=T),
            es.t0 = erreur.std (t0),
            es.t1 = erreur.std (t1),
            es.t2 = erreur.std (t2),
            es.t3 = erreur.std (t3),
            es.t4 = erreur.std (t4),
            es.t5 = erreur.std (t5),
            es.t6 = erreur.std (t6),
            es.pH0 = erreur.std (pH0),
            es.pH1 = erreur.std (pH1),
            es.pH2 = erreur.std (pH2),
            es.pH3 = erreur.std (pH3),
            es.pH4 = erreur.std (pH4),
            es.pH5 = erreur.std (pH5),
            es.pH6 = erreur.std (pH6),
            es.D0 = erreur.std(c(Da0,Db0)),
            es.D1 = erreur.std(c(Da1,Db1)),
            es.D2 = erreur.std(c(Da2,Db2)),
            es.D3 = erreur.std(c(Da3,Db3)),
            es.D4 = erreur.std(c(Da4,Db4)),
            es.D5 = erreur.std(c(Da5,Db5)),
            es.D6 = erreur.std(c(Da6,Db6)))
meltconsor <- melt(setDT(as.data.frame(m.consor)), measure = patterns("^m.t","^m.pH","^m.D","^es.t","^es.pH","^es.D"),
                  value.name = c("m.t","m.pH","m.D","es.t","es.pH","es.D"))  #???[,variable:= NULL][order(form.ID)]
meltconsor$variable <- as.character(meltconsor$variable)
meltconsor$point <- as.character(meltconsor$variable)
meltconsor <- mutate(meltconsor,sum.D=m.D,esum.D=es.D)
for ( i in 1:length(meltconsor$Consortia)){
  j <- meltconsor$Consortia[i]
  p <- meltconsor$point[i]
  if(as.character(j) == "BS"){
    meltconsor$variable[i] <- "Melange"
    B <- filter(meltconsor,Consortia=="BS",Souche=="B",point== p)$m.D
    S <- filter(meltconsor,Consortia=="BS",Souche=="S",point== p)$m.D
    es.B <- filter(meltconsor,Consortia=="BS",Souche=="B",point== p)$es.D
    es.S <- filter(meltconsor,Consortia=="BS",Souche=="S",point== p)$es.D
    meltconsor$sum.D[i] <- S+B
    meltconsor$esum.D[i] <- es.B+es.S
      } else if(as.character(j) == "BC") {
    meltconsor$variable[i] <- "Melange"
    B <- filter(meltconsor,Consortia=="BC",Souche=="B",point== p)$m.D
    C <- filter(meltconsor,Consortia=="BC",Souche=="C",point== p)$m.D
    es.B <- filter(meltconsor,Consortia=="BC",Souche=="B",point== p)$es.D
    es.C <- filter(meltconsor,Consortia=="BC",Souche=="C",point== p)$es.D
    meltconsor$sum.D[i] <- C+B
    meltconsor$esum.D[i] <- es.B+es.C
      } else {
    meltconsor$variable[i] <- "Isole"
  }
}
meltconsor$Consortia <- as.factor(meltconsor$Consortia)
meltconsor <- mutate(meltconsor,Essai=paste(Souche, Consortia))
meltconsor <- mutate(meltconsor,ISB=Souche)
for ( i in 1:length(meltconsor$Consortia)){
  j <- meltconsor$ISB[i]
  if(meltconsor$ISB[i] == "B"){
    meltconsor$ISB[i] <- "B"
  } else {
    meltconsor$ISB[i] <- "Autre"
  }
}
##################################################
###GRAPH pH en fonction du temps###
##################################################
meltconsor %>% filter(Essai =="B B"| Essai =="B BS"|Essai =="S S"|Essai =="T T") %>%
  ggplot(aes(x=m.t,y=m.pH,color=Consortia))+
  geom_line() +
  # geom_smooth(se=F,span=0.7)+
  geom_point(size=1)+
  geom_pointrange(aes(ymin=m.pH+es.pH,ymax=m.pH+es.pH),fatten=1)+
  geom_errorbarh(aes(xmax=m.t+es.t,xmin =m.t-es.t, height = 0))+
  theme_tufte() +
  # scale_y_continuous(limits=c(4,7.5))+
  scale_x_continuous(limits=c(0,30))+
  scale_color_manual(values = c("gold2","darkorange3","red3","black"))
meltconsor %>% filter(Essai =="B B"| Essai =="B BC"|Essai =="C C"|Essai =="T T") %>%
  ggplot(aes(x=m.t,y=m.pH,color=Essai))+
  geom_line() +
  # geom_smooth(se=F,span=0.7)+
  geom_point(size=1)+
  geom_pointrange(aes(ymin=m.pH+es.pH,ymax=m.pH+es.pH),fatten=1)+
  geom_errorbarh(aes(xmax=m.t+es.t,xmin =m.t-es.t, height = 0))+
  theme_tufte() +
  scale_y_continuous(limits=c(4,7.5))+
  scale_x_continuous(limits=c(0,30))+
  scale_color_manual(values = c("gold2","green4","blue2","black"))
meltconsor %>% filter(Essai =="B B"| Essai =="B BC"|Essai =="C C"|Essai =="T T"| Essai =="B BS"|Essai =="S S") %>%
  ggplot(aes(x=m.t,y=m.pH,color=Essai))+
  geom_line() +
  # geom_smooth(se=F,span=0.7)+
  geom_point(size=1)+
  geom_pointrange(aes(ymin=m.pH+es.pH,ymax=m.pH+es.pH),fatten=1)+
  geom_errorbarh(aes(xmax=m.t+es.t,xmin =m.t-es.t, height = 0))+
  theme_tufte() +
  scale_y_continuous(limits=c(4,7.5))+
  scale_x_continuous(limits=c(0,30))+
  scale_color_manual(values = c("gold2","green4","darkorange3","blue2","red3","black"))
##################################################
### GRAPH Num en fonction du temps ###
##################################################
plot4 <- meltconsor %>% filter(Consortia!= "T") %>%
  filter(Essai == "B B"|Essai =="S S"|Essai =="B BS"|Essai =="S BS") %>%
  ggplot(aes(x=m.t,y=m.D,color=Essai,linetype= variable, shape=Souche))+
  geom_line()+
  geom_point(size=1)+
  # scale_y_continuous(limits=c(10^4,10^9),trans='log10',position = "left")+
  # scale_x_continuous(limits=c(0,30))+
  geom_pointrange(aes(ymin=m.D-es.D,ymax=m.D+es.D),fatten=1)+
  geom_errorbarh(aes(xmax=m.t+es.t,xmin =m.t-es.t, height = 0))+
  theme_tufte() +
  theme(legend.position = "none")
  scale_shape_manual(values=c(19,19,20))+
  scale_color_manual(values = c("gold2","gold2","red3","red3")) +
  # scale_linetype_manual( values = c(1,3,1))
  # facet_wrap(~Souche, ncol=1)
meltconsor %>% filter(Consortia!= "T") %>%
  filter(Essai == "B B"|Essai =="C C"|Essai =="B BC"|Essai =="C BC") %>%
  ggplot(aes(x=m.t,y=m.D,color=Consortia,linetype= variable, shape=Souche))+
  geom_line()+
  geom_point(size=1)+
  scale_y_continuous(limits=c(10^5,10^8),trans='log10')+
  scale_x_continuous(limits=c(0,30))+
  geom_pointrange(aes(ymin=m.D-es.D,ymax=m.D+es.D),fatten=1)+
  geom_errorbarh(aes(xmax=m.t+es.t,xmin =m.t-es.t, height = 0))+
  theme_tufte() +
  scale_shape_manual(values=c(19,19,20))+
  scale_color_manual(values = c("gold2","green4","blue2")) +
  scale_linetype_manual( values = c(1,3,1))+
  facet_wrap(~Souche, ncol=1)
meltconsor %>% filter(Essai == "B B"|Essai =="C C"|Essai =="B BC"|Essai =="C BC") %>%
  ggplot(aes(x=m.t,y=sum.D,color=Consortia))+
  geom_line()+
  geom_point(size=1)+
  scale_y_continuous(limits=c(10^3,10^9),trans='log10')+
  scale_x_continuous(limits=c(0,30))+
  geom_pointrange(aes(ymin=sum.D-esum.D,ymax=sum.D+esum.D),fatten=1)+
  geom_errorbarh(aes(xmax=m.t+es.t,xmin =m.t-es.t, height = 0))+
  theme_tufte() +
  scale_shape_manual(values=c(19,19))+
  scale_color_manual(values = c("gold2","green4","blue2"))
meltconsor %>% filter(Consortia!= "T") %>%
  ggplot(aes(x=m.t,y=sum.D,color=Consortia))+
  geom_line()+
  geom_point(size=1)+
  scale_y_continuous(limits=c(10^4,10^9),trans='log10')+
  scale_x_continuous(limits=c(0,30))+
  geom_pointrange(aes(ymin=sum.D-esum.D,ymax=sum.D+esum.D),fatten=1)+
  geom_errorbarh(aes(xmax=m.t+es.t,xmin =m.t-es.t, height = 0))+
  theme_tufte() +
  scale_shape_manual(values=c(19,19))+
  scale_color_manual(values = c("gold2","green4","darkorange3","blue2","red3"))
##################################################
### GRAPH Num en fonction du pH ###
##################################################
plot3 <- meltconsor %>% filter(Consortia!= "T") %>%
  filter(Essai == "B B"|Essai =="S S"|Essai =="B BS"|Essai =="S BS") %>%
  ggplot(aes(x=m.pH,y=m.D,color=Essai,linetype= variable, shape=Souche))+
  geom_line()+
  geom_point(size=1)+
  scale_y_continuous(limits=c(10^4,10^9),trans='log10')+
  scale_x_continuous(limits=c(4,7))+
  geom_pointrange(aes(ymin=m.D-es.D,ymax=m.D+es.D),fatten=1)+
  geom_errorbarh(aes(xmax=m.pH+es.pH,xmin =m.pH-es.pH, height = 0))+
  theme_tufte() +
  scale_shape_manual(values=c(19,19))+
  scale_color_manual(values = c("gold2","gold2","red3","red3")) +
  scale_linetype_manual( values = c(1,3,1))+
  # facet_wrap(~Souche, ncol=1)+
  scale_x_reverse()
meltconsor %>% filter(Consortia!= "T") %>%
  filter(Essai == "B B"|Essai =="C C"|Essai =="B BC"|Essai =="C BC") %>%
  ggplot(aes(x=m.pH,y=m.D,color=Consortia,linetype= variable, shape=Souche))+
  geom_line()+
  geom_point(size=1)+
  scale_y_continuous(limits=c(10^5,10^8),trans='log10')+
  scale_x_continuous(limits=c(4,7))+
  geom_pointrange(aes(ymin=m.D-es.D,ymax=m.D+es.D),fatten=1)+
  geom_errorbarh(aes(xmax=m.pH+es.pH,xmin =m.pH-es.pH, height = 0))+
  theme_tufte() +
  scale_shape_manual(values=c(19,20))+
  scale_color_manual(values = c("gold2","green4","blue2")) +
  scale_linetype_manual( values = c(1,3,1))+
  facet_wrap(~Souche, ncol=1)+
  scale_x_reverse()
meltconsor %>% filter(Essai == "B B"|Essai =="C C"|Essai =="B BC"|Essai =="C BC") %>%
  ggplot(aes(x=m.pH,y=sum.D,color=Consortia))+
  geom_line()+
  geom_point(size=1)+
  scale_y_continuous(limits=c(10^3,10^9),trans='log10')+
  scale_x_continuous(limits=c(4,7))+
  geom_pointrange(aes(ymin=sum.D-esum.D,ymax=sum.D+esum.D),fatten=1)+
  geom_errorbarh(aes(xmax=m.pH+es.pH,xmin =m.pH-es.pH, height = 0))+
  theme_tufte() +
  scale_shape_manual(values=c(19,19))+
  scale_color_manual(values = c("gold2","green4","blue2")) +
  scale_x_reverse()
meltconsor %>% filter(Consortia!= "T") %>%
  ggplot(aes(x=m.pH,y=sum.D,color=Consortia))+
  geom_line()+
  geom_point(size=1)+
  scale_y_continuous(limits=c(10^3,10^9),trans='log10')+
  scale_x_continuous(limits=c(4,7))+
  geom_pointrange(aes(ymin=sum.D-esum.D,ymax=sum.D+esum.D),fatten=1)+
  geom_errorbarh(aes(xmax=m.pH+es.pH,xmin =m.pH-es.pH, height = 0))+
  theme_tufte() +
  scale_shape_manual(values=c(19,19))+
  scale_color_manual(values = c("gold2","green4","darkorange3","blue2","red3")) +
  scale_x_reverse()
##################################################
### GRAPH Num en fonction du pH + tukey ###
##################################################
meltconsor %>% filter(Consortia!= "T") %>%
  filter(Essai == "B B"|Essai =="S S"|Essai =="B BS"|Essai =="S BS") %>%
  ggplot(aes(x=m.pH,y=m.D,color=Consortia,linetype= variable, shape=Souche))+
  geom_line()+
  geom_point(size=1)+
  scale_y_continuous(limits=c(10^3,10^9),trans='log10')+
  geom_pointrange(aes(ymin=m.D-es.D,ymax=m.D+es.D),fatten=1)+
  geom_errorbarh(aes(xmax=m.pH+es.pH,xmin =m.pH-es.pH, height = 0))+
  theme_tufte() +
  scale_shape_manual(values=c(19,19))+
  scale_color_manual(values = c("gold2","darkorange3","red3")) +
  scale_linetype_manual( values = c(1,3,1))+
  facet_wrap(~Souche, ncol=1)+
  geom_text(data =filter(as.data.frame(df_tukey_letters2),Essai=="D1"|"D2"|"D3"|"D4"|"D5")), aes(x=trt, y=ymax1,label=groups),vjust=0.2,hjust=0.2,angle=90, size = 4,position=position_dodge(.5))+
  scale_x_reverse()
meltconsor %>% filter(Consortia!= "T") %>%
  ggplot(aes(x=m.pH,y=sum.D,color=Consortia))+
  geom_line()+
  geom_point(size=1)+
  scale_y_continuous(limits=c(10^3,10^9),trans='log10')+
  geom_pointrange(aes(ymin=sum.D-esum.D,ymax=sum.D+esum.D),fatten=1)+
  geom_errorbarh(aes(xmax=m.pH+es.pH,xmin =m.pH-es.pH, height = 0))+
  theme_tufte() +
  scale_shape_manual(values=c(19,19))+
  scale_color_manual(values = c("gold2","darkorange3","red3")) +
  scale_x_reverse()
##################################################
### Extraction fit -> stock.para ###
##################################################
stock.para <- NULL
for (i in unique(meltconsor$Essai)) {
  sub.data <- meltconsor %>% filter(Essai == i)
  # sub.data <- select(sub.data, time = Time,intensity = pH)
  # sub.co2 <- myco2 %>% filter(Especes == i)
  # for (j in unique(sub.data$Essai)) {
  #     sub.sub <- sub.data %>% filter(Essai == j)
  #     ss.co2 <- sub.co2 %>% filter(Essai == j)
  # if (mean(ss.sub$PopV.Sacc) != 0 & mean(ss.sub$PopV.NonSacc) == 0) {
  #       res.sacc <- gcFitModel(ss.sub$Temps,ss.sub$PopV.Sacc)
  #       mypara <- data.frame(Especes = i, 
  #                            # Essai = j, 
  #                            mu.sacc = res.sacc$parameters$mu[1],
  #                            K.sacc = res.sacc$parameters$A[1],
  #                            latency.sacc = res.sacc$parameters$lambda[1],
  #                            mu.nonsacc = NA,
  #                            K.nonsacc = NA,
  #                            latency.nonsacc = NA)
  #       mypara$tScdomin <- 0
  #   }
  # if (mean(ss.sub$PopV.Sacc) == 0 & mean(ss.sub$PopV.NonSacc) != 0) {
  #   res.nonsacc <- gcFitModel(ss.sub$Temps,ss.sub$PopV.NonSacc)
  #   
  #   mypara <- data.frame(Especes = i, 
  #                        Essai = j, 
  #                        Fermenteur = k, 
  #                        mu.sacc = NA,
  #                        K.sacc = NA,
  #                        latency.sacc = NA,
  #                        mu.nonsacc = res.nonsacc$parameters$mu[1],
  #                        K.nonsacc = res.nonsacc$parameters$A[1],
  #                        latency.nonsacc = res.nonsacc$parameters$lambda[1])
  #   mypara$tScdomin <- max(ss.sub$Temps)       
  # }
  # if (mean(ss.sub$PopV.Sacc) != 0 & mean(ss.sub$PopV.NonSacc) != 0) {
  paste("res",i) <- gcFitModel(filter(sub.data, t!=0)$t,filter(sub.data, t!=0)$D)
  normlized_sub.data <- normalizeData(sub.data[1:450,])
  parameterVector <- multipleFitFunction(normlized_sub.data, model = "sigmoidal")
  intensityTheoretical <- sigmoidalFitFormula(sub.data$time[1:450],
                                              maximum = parameterVector$maximum_Estimate,
                                              slopeParam = parameterVector$slopeParam_Estimate,
                                              midPoint = parameterVector$midPoint_Estimate)
  comparisonData <- cbind(sub.data[1:450,], intensityTheoretical)
  ggplot(comparisonData)+
    geom_point(aes(x = time, y = intensity)) +
    geom_line(aes(x = time, y = intensityTheoretical), color = "orange") +
    expand_limits(x = 0, y = 0)
  # res.nonsacc <- gcFitModel(ss.sub$Temps,ss.sub$PopV.NonSacc)
  mypara <- data.frame(Variable = i, 
                       # Essai = j, 
                       mu.pH = paste("res",i)$parameters$mu[1],
                       K.sacc = paste("res",i)$parameters$A[1],
                       latency.sacc = paste("res",i)$parameters$lambda[1],
                       # mu.nonsacc = res.nonsacc$parameters$mu[1],
                       # K.nonsacc = res.nonsacc$parameters$A[1],
                       # latency.nonsacc = res.nonsacc$parameters$lambda[1])
                       # mypara$tScdomin <- min(filter(ss.sub, Prop_Scvie == ss.sub$Prop_Scvie[which.min(abs(ss.sub$Prop_Scvie - 0.5))])$Temps)
}

mypara$tVmax <- NA
mypara$t80 <- NA
mypara$mu.pop <- NA
mypara$K.pop <- NA
mypara$latency.pop <- NA
mypara$mu.co2 <- NA
mypara$K.co2 <- NA
mypara$latency.co2 <- NA
res.co2 <- NA
res.pop <- NA
mypara$tfin.co2 <- NA
mypara$co2max <- NA
mypara$popmax <- NA
mypara$tVmax <- mean(filter(sss.co2, debitCO2 == max(sss.co2$debitCO2))$Temps)
mypara$t80 <- mean(filter(sss.co2, cumul == sss.co2$cumul[which.min(abs(sss.co2$cumul - (0.8*max(sss.co2$cumul))))],sss.co2$cumul)$Temps)
res.pop <- gcFitModel(ss.sub$Temps,ss.sub$PopV)
mypara$mu.pop <- res.pop$parameters$mu[1]
mypara$K.pop <- res.pop$parameters$A[1]
mypara$latency.pop <- res.pop$parameters$lambda[1]
res.co2 <- gcFitModel(sss.co2$Temps, sss.co2$cumul)
mypara$mu.co2 <- res.co2$parameters$mu[1]
mypara$K.co2 <- res.co2$parameters$A[1]
mypara$latency.co2 <- res.co2$parameters$lambda[1]
mypara$tfin.co2 <- mean(filter(sss.co2, cumul == sss.co2$cumul[which.min(abs(sss.co2$cumul-75))])$Temps)
mypara$co2max <- max(sss.co2$cumul)
mypara$popmax <- max(ss.sub$Pop_totale)

stock.para <- rbind.data.frame(stock.para,mypara)
}

}
}


m.mydata %>%
  select(m.PopV.Souche,es.PopV.Souche, Temps, Souche, Essai, isole, Especes,color) %>%
  ggplot(aes(x= Temps, y= m.PopV.Souche, color = color,  linetype= isole )) +
  geom_line()+
  geom_point() +
  geom_pointrange(aes(ymin= m.PopV.Souche-es.PopV.Souche, ymax=m.PopV.Souche+es.PopV.Souche), fatten=1) +
  scale_color_manual(values = c("darkorange1", "red","green3","turquoise","black","darkorange3","red3","green4","turquoise3","darkorchid4", "darkorchid1")) +
  scale_linetype_manual( values = c(1,3)) +
  theme_tufte()+
  facet_wrap(Souche~Essai, ncol=5)

summary(as.factor(Criblage$specie))

sum(summary(as.factor(filter(Criblage,`NS_48`>(`Sac_48`+0.3))$specie)))


summary(as.factor(Criblage$specie))[-c(7,8,12,18,21,22,24)]

filter(Criblage, specie == "L. curvatus")$`Raf_48`
filter(Criblage, `NS_48`>(`Raf_48`+0.3))

filter(Criblage, `NS_10`>=(`Sta_10`+0.3),`NS_10`>(`Raf_10`+0.3))

i <- filter(Criblage, RefCIRM == 775)

filter(Criblage, specie == "L. johnsonii"| specie=="L. acidophilus")%>%
  select(specie,NS_10, `Sta_10`, `550_10`,`NS_48`,`Sta_48`,`550_48`)


filter(Criblage, `550_10` < 6, specie =="S. thermophilus")

write_xlsx(ph, path ="~/Dossiers Olivier HARLE/INRA/Souches/TabR/Souchespositivesacidi.xlsx")

is.na(Criblage$`550_48`)
tab <- select(Criblage, RefCIRM, specie,`550_48`) %>% group_by(specie) %>%
  summarize(`550`= mean(`550_48`))
print(as.data.frame(tab),topn=25)

Criblage$specie <- as.factor(Criblage$specie)
# levels(Criblage$specie) <- c(levels(Criblage$specie)[-c(1,22)],levels(Criblage$specie)[c(1,22)])

fig1 <- filter( Criblage, specie != "T") %>%
  ggplot(aes(x=specie, y=`550_10`, color=specie )) +
  geom_boxplot()+
  theme_tufte()+
  scale_y_continuous(limits=c(4,7)) +
  theme(legend.position ="None",axis.text.x = element_text(size=16,hjust = 1,vjust=0.1, angle = 90, face = "italic")) +
  # geom_hline(yintercept = 6 , linetype = 1) +
  labs(x = "Especes", y = "pH", title = "Jus 550_10h")

Criblage$specie <- factor(c(Criblage$specie), levels = c("L. acidophilus","L. amylovorus", "L. casei","L. coryniformis","L. curvatus","L. delbrueckii","L. diolivorans","L. helveticus","L. johnsonii","L. kunkeei","L. lactis","L. mali","L. mesenteroides","L. paracasei","L. paraplantarum","L. pentosus","L. plantarum","L. pontis","L. rhamnosus","L. sanfranciscensis","L.  xiangfangensis","L. zeae","T","S. thermophilus"))

jpeg("Fig1b0219.jpeg", width = 12, height = 6, units = 'in', res = 1000)
filter( Criblage, specie != "T") %>%
  ggplot(aes(x=specie, y=`550_48`, color=specie )) +
  geom_boxplot()+
  theme_tufte()+
  scale_y_continuous(limits=c(4,7)) +
  theme(legend.position ="None",axis.text.x = element_text(size=16,hjust = 1,vjust=0.1, angle = 90, face = "italic")) +
  # geom_hline(yintercept = 6 , linetype = 1) +
  labs(x = "", y = "", title = "Jus 550_48h")

dev.off()

################################################
################################################
################TUKEY SUCRES####################
################################################
################################################
consorS <- select(consor,Souche,Consortia,Replicat,names(consor)[(length(names(consor))-6):length(names(consor))]) %>%
  filter(Consortia=="B"|Consortia=="S"|Consortia=="BS"|Consortia=="T")
consorS <-filter(consorS, Consortia!="BS"|Souche!="S")
consorS <-filter(consorS, Replicat == "a"|Replicat == "b"|Replicat == "c")

#####POUR PASSER EN MOL/L#####
consorS <- mutate(consorS, sucrosepH4.5=sucrosepH4.5/342.962,raffinosepH4.5=raffinosepH4.5/504.42,verbascosepH4.5=verbascosepH4.5/828.7183,stacchyosepH4.5=stacchyosepH4.5/666.5777,glucosepH4.5=glucosepH4.5/180.1559,fructosepH4.5=fructosepH4.5/180.1559,galactosepH4.5=galactosepH4.5/180.1559)
##############################
mypara1 <- names(consorS)[4:length(names(consorS))]

mymodel1 <- lapply(mypara1,function(k) {  #D?finition de la fonction test lm
  lm(eval(substitute(j ~ Consortia,list(j=as.name(k)))),data = consorS )  
}
)
tabgather1 <- consorS[,-c(1)] %>%        #POUR LES SUCRES
  gather(key = parametres, val = value,-Consortia,-Replicat)

####???Début tukey#### Première boucle print numéro de la variable de mypara1 traitée
variable_list1 <- list()
aov1 <- list()
tukey1 <- list()
tukey_groups1 <- list()
for(i in 1:length(mypara1)) { 
  variable_list1[[i]] <- tabgather1[tabgather1$parametres == mypara1[i],] 
  print(i)
  for (i in 1:length(variable_list1)) {
    aov1[[i]] <- (aov(value~Consortia,data=variable_list1[[i]]))
    for (i in 1:length(aov1)) {
      tukey1[[i]] <- HSD.test(aov1[[i]],trt="Consortia",group=FALSE, alpha = 0.1)
      tukey_groups1[[i]] <- HSD.test(aov1[[i]], trt = "Consortia", alpha = 0.1)
    }
  }
}
####données statistique de l'anova
stats_data1 <- list()
for (i in 1:length(variable_list1)) {
  stats_data1[[i]] <- summarySE(variable_list1[[i]], measurevar = "value", groupvars = c("Consortia"),na.rm = TRUE)
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
    xpos.t1[k] <- as.integer(rownames(stats_data1[[i]][stats_data1[[i]]$Consortia==(as.vector(row.names(tukey_letters1[[i]])[j])),]))
    k<-k+1
    tukey_letters1[[i]]$Especes[j] <- filter(variable_list1[[1]], Consortia == tukey_letters1[[i]]$trt[j])[1,3]
    tukey_letters1[[i]]$Names[j] <- filter(variable_list1[[1]], Consortia == tukey_letters1[[i]]$trt[j])[1,2]
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
df_tukey_letters1
#Create a custom color scale
myColors <- c("gold2","darkorange3","red3","black")
df_tukey_letters1$trt <- factor(df_tukey_letters1$trt, levels = c("B","BS","S","T"))
names(myColors) <- levels(df_tukey_letters1$trt)
fillScale <- scale_fill_manual(name = "trt",values = myColors)

####BOUCLES POUR TRACER TOUT LES GRAPHIQUES####
for(i in mypara1){     
  print(filter(df_tukey_letters1,group == i) %>%
          ggplot(aes(x= trt, y= value, ymax = ymax1, ymin=ymin1, fill=trt)) +
          geom_bar(stat = "identity")+
          fillScale +
          #facet_wrap(~group, scales = "free", ncol = 5) +
          geom_errorbar() +
          theme_tufte() +
          scale_y_continuous(limits = c(0, 0.001)) +
          theme (axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.3,size=15),title = element_text(size=15),axis.text.y = element_text(size=15),legend.text = element_text(size=15,face = "italic"))+   #,legend.position = "bottom" 
          geom_text(data =filter(df_tukey_letters1,group == i), aes(x=trt, y=ymax1,label=groups),vjust=0.2,hjust=0.2,angle=90, size = 4,position=position_dodge(.5))+
          #geom_text(aes(label=round(value, digits=2)), vjust=3, color="black",position = position_dodge(0.9), size=2) +
          labs(x = "Consortia", y = "en mol/L", title =i))
}
