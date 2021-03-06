rm(list = ls()) 
library(tidyverse)
library(gridExtra)
library(caret)
library(ggthemes)
library(stringr)

setwd("C:/Users/harleoli/Google Drive/Stage Olivier/Cyto_SceVSTde")
setwd("/Users/thibaultnidelet_boulot/Google Drive/Stage Olivier/Cyto_SceVSTde")



names <- list.files("C:/Users/harleoli/Google Drive/Stage Olivier/Cyto_SceVSTde")[-1]
names <- list.files("/Users/thibaultnidelet_boulot/Google Drive/Stage Olivier/Cyto_SceVSTde")[-1]



all.temps <- NULL
stock.viabilite <- NULL
for (i in names) {
    mytemps <- str_split(i,"_")[[1]][2]
    # mytemps <- str_split(mytemps,"h")[[1]][1]
        
    all.temps <- c(all.temps, mytemps)
}
all.temps <- unique(all.temps)

#i = "24"

#all.temps <- sample(all.temps,5)
# all.temps <- all.temps[3]

pdf("prediction_machine_learningtd2.pdf")
for (i in all.temps) {
    names.sub <- names[str_detect(names,as.character(i))]
    
    names.sc <- names.sub[str_detect(names.sub,"_F21\\.") | str_detect(names.sub,"_F31\\.")]
    data.sc <- NULL
    for (j in names.sc) {
      data.sc <- rbind(data.sc, read.csv(file = j))
    }
    data.sc$souche <- "sc"
    data.sc$vivante <- "Oui"
    data.sc$vivante[data.sc$FL3.A > 10000] <- "Non"
    data.sc <- data.sc %>% 
        select(-Time, -Width) %>%
        filter(FSC.A > quantile(data.sc$FSC.A,0.01) & FSC.A < quantile(data.sc$FSC.A,0.99)) %>%
        filter(SSC.A > quantile(data.sc$SSC.A,0.01) & SSC.A < quantile(data.sc$SSC.A,0.99)) %>%
        filter(FL1.A > quantile(data.sc$FL1.A,0.01) & FL1.A < quantile(data.sc$FL1.A,0.99)) %>%
        filter(FL2.A > quantile(data.sc$FL2.A,0.01) & FL2.A < quantile(data.sc$FL2.A,0.99)) %>%
        filter(FL3.A > quantile(data.sc$FL3.A,0.01) & FL3.A < quantile(data.sc$FL3.A,0.99)) %>%
        filter(FL4.A > quantile(data.sc$FL4.A,0.01) & FL4.A < quantile(data.sc$FL4.A,0.99)) %>%
        filter(FSC.H > quantile(data.sc$FSC.H,0.01) & FSC.H < quantile(data.sc$FSC.H,0.99)) %>%
        filter(SSC.H > quantile(data.sc$SSC.H,0.01) & SSC.H < quantile(data.sc$SSC.H,0.99)) %>%
        filter(FL1.H > quantile(data.sc$FL1.H,0.01) & FL1.H < quantile(data.sc$FL1.H,0.99)) %>%
        filter(FL2.H > quantile(data.sc$FL2.H,0.01) & FL2.H < quantile(data.sc$FL2.H,0.99)) %>%
        filter(FL3.H > quantile(data.sc$FL3.H,0.01) & FL3.H < quantile(data.sc$FL3.H,0.99)) %>%
        filter(FL4.H > quantile(data.sc$FL4.H,0.01) & FL4.H < quantile(data.sc$FL4.H,0.99))
    
    names.td <- names.sub[str_detect(names.sub,"_F25\\.") | str_detect(names.sub,"_F26\\.") | str_detect(names.sub,"_F32\\.")]
    data.td <- NULL
    for (j in names.td) {
        data.td <- rbind(data.td, read.csv(file = j))
    }
    data.td$souche <- "td"
    data.td$vivante <- "Oui"
    data.td$vivante[data.td$FL3.A > 10000] <- "Non"
    data.td <- data.td %>% 
        select(-Time, -Width) %>%
        filter(FSC.A > quantile(data.td$FSC.A,0.01) & FSC.A < quantile(data.td$FSC.A,0.99)) %>%
        filter(SSC.A > quantile(data.td$SSC.A,0.01) & SSC.A < quantile(data.td$SSC.A,0.99)) %>%
        filter(FL1.A > quantile(data.td$FL1.A,0.01) & FL1.A < quantile(data.td$FL1.A,0.99)) %>%
        filter(FL2.A > quantile(data.td$FL2.A,0.01) & FL2.A < quantile(data.td$FL2.A,0.99)) %>%
        filter(FL3.A > quantile(data.td$FL3.A,0.01) & FL3.A < quantile(data.td$FL3.A,0.99)) %>%
        filter(FL4.A > quantile(data.td$FL4.A,0.01) & FL4.A < quantile(data.td$FL4.A,0.99)) %>%
        filter(FSC.H > quantile(data.td$FSC.H,0.01) & FSC.H < quantile(data.td$FSC.H,0.99)) %>%
        filter(SSC.H > quantile(data.td$SSC.H,0.01) & SSC.H < quantile(data.td$SSC.H,0.99)) %>%
        filter(FL1.H > quantile(data.td$FL1.H,0.01) & FL1.H < quantile(data.td$FL1.H,0.99)) %>%
        filter(FL2.H > quantile(data.td$FL2.H,0.01) & FL2.H < quantile(data.td$FL2.H,0.99)) %>%
        filter(FL3.H > quantile(data.td$FL3.H,0.01) & FL3.H < quantile(data.td$FL3.H,0.99)) %>%
        filter(FL4.H > quantile(data.td$FL4.H,0.01) & FL4.H < quantile(data.td$FL4.H,0.99))
    
    data.both <- rbind(data.sc,data.td)
    
    G.both <- ggplot(data.both, aes(x = FL1.H, y = FSC.A, color = souche)) + 
        scale_y_log10() +
        scale_x_log10() +
        theme_tufte() +
        geom_point(alpha = 0.2) +
        labs(title = "Culture isole")

    ggplot(data.sc, aes(x = FL1.H, y = FSC.A)) +
        scale_y_log10() +
        scale_x_log10() +
        theme_tufte() +
        geom_point(alpha = 0.2) +
        labs(title = "Culture isole")

    ggplot(data.td, aes(x = FL1.H, y = FSC.A)) +
        scale_y_log10() +
        scale_x_log10() +
        theme_tufte() +
        geom_point(alpha = 0.2) +
        labs(title = "Culture isole")

    
    inTrain <- createDataPartition(y = data.both$souche,p = 0.7, list = FALSE)
    training <- data.both[inTrain,]
    testing <- data.both[-inTrain,]

    #modFit <- train(souche ~ .,data = training,method = "rpart")
    #modFit <- train(souche ~ ., method = "gbm",data = training,verbose = FALSE)
    
    modFit <- train(souche ~ ., method = "gbm",data = data.both[-14],verbose = FALSE)
    my_pred <- predict(modFit,newdata = testing)
    testing$predRight <- my_pred == testing$souche
    print(sum(testing$predRight)/nrow(testing)*100)
    data.both$souche <- as.factor(data.both$souche)
    
    names.comp <- names.sub[str_detect(names.sub,"_F22\\.") | str_detect(names.sub,"_F23\\.") | str_detect(names.sub,"_F24\\.")]
    data.comp <- NULL
    for (j in names.comp) {
        data.int <- read.csv(file = j)
        data.int$vivante <- "Oui"
        data.int$vivante[data.int$FL3.A > 10000] <- "Non"
        
        
        
        data.int <- data.int %>% 
            select(-Time, -Width) %>%
            filter(FSC.A > quantile(data.int$FSC.A,0.01) & FSC.A < quantile(data.int$FSC.A,0.99)) %>%
            filter(SSC.A > quantile(data.int$SSC.A,0.01) & SSC.A < quantile(data.int$SSC.A,0.99)) %>%
            filter(FL1.A > quantile(data.int$FL1.A,0.01) & FL1.A < quantile(data.int$FL1.A,0.99)) %>%
            filter(FL2.A > quantile(data.int$FL2.A,0.01) & FL2.A < quantile(data.int$FL2.A,0.99)) %>%
            filter(FL3.A > quantile(data.int$FL3.A,0.01) & FL3.A < quantile(data.int$FL3.A,0.99)) %>%
            filter(FL4.A > quantile(data.int$FL4.A,0.01) & FL4.A < quantile(data.int$FL4.A,0.99)) %>%
            filter(FSC.H > quantile(data.int$FSC.H,0.01) & FSC.H < quantile(data.int$FSC.H,0.99)) %>%
            filter(SSC.H > quantile(data.int$SSC.H,0.01) & SSC.H < quantile(data.int$SSC.H,0.99)) %>%
            filter(FL1.H > quantile(data.int$FL1.H,0.01) & FL1.H < quantile(data.int$FL1.H,0.99)) %>%
            filter(FL2.H > quantile(data.int$FL2.H,0.01) & FL2.H < quantile(data.int$FL2.H,0.99)) %>%
            filter(FL3.H > quantile(data.int$FL3.H,0.01) & FL3.H < quantile(data.int$FL3.H,0.99)) %>%
            filter(FL4.H > quantile(data.int$FL4.H,0.01) & FL4.H < quantile(data.int$FL4.H,0.99))
        
        mypred <- predict(modFit,newdata = data.int)
        data.int$souche <- mypred
        data.int$both <- paste(data.int$vivante,data.int$souche)
        data.comp <- rbind(data.comp, data.int)
        
        fermenteur <- str_split(j,"_")[[1]][3]
        fermenteur <- str_split(fermenteur,".csv")[[1]][1]
        
        viabilite <- data.frame(temps = i, rep = fermenteur,
                                via.sc = sum(data.sc$vivante == "Oui")/nrow(data.sc), 
                                via.td = sum(data.td$vivante == "Oui")/nrow(data.td), 
                                via.comp.sc = sum(data.int$vivante == "Oui" & data.int$souche == "sc")/nrow(data.int[data.int$souche == "sc",]),
                                via.comp.td = sum(data.int$vivante == "Oui" & data.int$souche == "td")/nrow(data.int[data.int$souche == "td",]),
                                freq.sc = sum(data.int$souche == "sc")/nrow(data.int))
        stock.viabilite <- rbind.data.frame(stock.viabilite,viabilite)
        
        
       

        G.comp <- ggplot(data.int, aes(x = FL1.H, y = FSC.A, color = souche)) +
            scale_y_log10() +
            scale_x_log10() +
            theme_tufte() +
            geom_point(alpha = 0.2) +
            labs(title = "En competition")
        
        grid.arrange(G.both, G.comp,ncol = 2, top = paste(i,", fermenteur:",fermenteur,sep = ""))
        
    }
}
dev.off()

stock.viabilite

tab.taux.prediction <- data.frame(temps = all.temps,taux.prediction = c(99.5818, 100,99.99172,99.8333,98.97172,100,93.95126,100,100,99.99636,100,99.97897,99.88365))

write.csv(stock.viabilite, file = "tab_viabilite.csv")
write.csv(tab.taux.prediction, file = "taux_prediction.csv")

for (i in names) {
    mydata <- read.csv(file = paste("C:/Users/harleoli/Google Drive/Stage Olivier/Cyto_SceVSTde",i,sep = "/"))
    
    mydata$vivante <- "Oui"
    mydata$vivante[mydata$FL3.A > 10000] <- "Non"
    viabilite <- sum(mydata$vivante == "Oui")/nrow(mydata)
    
    # mydata <- mydata %>%
    #     filter(FL3.A != 0)
    
    # graph <- mydata %>%
    #     ggplot(aes(x = FL3.A, y = FL1.A, color = vivante)) +
    #     scale_y_log10() +
    #     scale_x_log10() +
    #     theme_tufte() +
    #     geom_vline(xintercept = 10000,col = 2) +
    #     geom_point(alpha = 0.2) +
    #     labs(title = i)
    # print(graph)
    
    
}

global.data <- NULL
sub_names <- names[str_detect(names, "F21B")]
for(i in sub_names) {
    newdata <- read.csv(file = paste("C:/Users/harleoli/Google Drive/Stage Olivier/Cyto_SceVSTde",i,sep = "/"))
    newdata$souche <- "Sc"
    global.data <- rbind(global.data, newdata)
}

sub_names <- c(names[str_detect(names, "F22B")],names[str_detect(names, "F23B")],names[str_detect(names, "F24B")])
for(i in sub_names) {
    newdata <- read.csv(file = paste("C:/Users/harleoli/Google Drive/Stage Olivier/Cyto_SceVSTde",i,sep = "/"))
    newdata$souche <- "Mp"
    global.data <- rbind(global.data, newdata)
}

global.data.vivante <- global.data %>%
    filter(FL3.A < 10000)



inTrain <- createDataPartition(y=global.data.vivante$souche,p = 0.7, list = FALSE)
training <- global.data.vivante[inTrain,]
testing <- global.data.vivante[-inTrain,]
        
model <- train(souche ~ .,data = training,method = "rpart")
my_pred <- predict(model,newdata=testing)
testing$predRight <- my_pred==testing$souche

table(my_pred,testing$souche)
sum(testing$predRight)/nrow(testing)*100 

plot(modFit$finalModel, uniform=TRUE,main="Classification Tree")
text(modFit$finalModel, use.n=TRUE, all=TRUE, cex=.8)
        

stock.taux <- NULL
sub_names <- c(names[str_detect(names, "F25B")],names[str_detect(names, "F26B")],names[str_detect(names, "F27B")])
for(i in sub_names) {
    newdata <- read.csv(file = paste("C:/Users/harleoli/Google Drive/Stage Olivier/Cyto_SceVSTde",i,sep = "/"))
    my_pred <- predict(modFit,newdata = newdata)
    taux.Sc <- sum(my_pred == "Sc")/length(my_pred)
    stock.taux <- rbind.data.frame(stock.taux, data.frame(i, taux.Sc))
}





#autre methode 
# modFit <- train(souche ~ ., method = "gbm",data=training,verbose=FALSE)
# print(modFit)
# 
# modFit <- train(Espece~ .,data=training,method="rf",prox=TRUE)
# modFit
# 
# 
# 
# 
# 
# pred <- predict(modFit,testing)
# testing$predRight <- pred == testing$souche
# table(pred,testing$souche)
# 
# 
# 
# 
# 
# 
# 
# fruit <- c("apple", "banana", "pear", "pinapple")
# str_subset(fruit, "a")
# str_which(fruit, "a")
# 
# str_subset(fruit, "^a")
# str_subset(fruit, "a$")
# str_subset(fruit, "b")
# str_subset(fruit, "[aeiou]")
# 
# # Missings never match
# str_subset(c("a", NA, "b"), ".")
# str_which(c("a", NA, "b"), ".")
# 
# 
# 
# plot(model$finalModel, uniform=TRUE,main="Classification Tree")
# text(model$finalModel, use.n=TRUE, all=TRUE, cex=.8)
# 
# 
# ls(model)
# 
# 
# i = "22.8"
