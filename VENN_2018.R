rm(list = ls()) #vider l'environnement


  
  ######DIAGRAMME DE VENN#######
library(VennDiagram)


###Sur toutes les souches
venn.plot <- draw.quad.venn(
  area1 =32,
  area2 =13,
  area3 =10,
  area4 =18,
  n12 = 8,
  n13 = 7,
  n14 = 8 ,
  n23 = 6,
  n24 = 5,
  n34 = 3,
  n123 = 3,
  n124 = 2,
  n134 = 2,
  n234 = 2,
  n1234=1 ,
  category = c("Acide", "Notes T", "Notes 37", "Gel+"),
  fill = c("orange", "red", "green", "blue"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("orange", "red", "green", "blue")
)

area1= 15 + 3 + 5 + 1 +2 +1 +1
area2= 4 + 3 + 1 + 1 + 2 + 1 + 1
area3= 3 + 3  + 1 + 1 + 1 + 1
area4= 7 + 5 + 1 + 2 + 1 + 1 + 1
n12= 1 + 2 + 1
n13= 3 + 1 + 1 + 1
n14= 5 + 2 + 1 + 1 + 1
n23= 3 + 1 + 1 + 1
n24= 1 + 2 + 1 + 1
n34= 1 + 1 + 1

###Sur les souches thermophiles
venn.plot <- draw.quad.venn(
  area1 =31,
  area2 =11,
  area3 =9,
  area4 =8,
  n12 = 8,
  n13 = 7,
  n14 = 7 ,
  n23 = 5,
  n24 = 3,
  n34 = 2,
  n123 = 3,
  n124 = 2,
  n134 = 2 ,
  n234 = 1,
  n1234=1 ,
  category = c("Acide", "Notes T", "Notes 37", "Gel+"),
  fill = c("orange", "red", "green", "blue"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("orange", "red", "green", "blue")
)


###Sur les souches mesophiles
venn.plot <- draw.quad.venn(
  area1 =1,
  area2 =2,
  area3 =1,
  area4 =10,
  n12 = 0,
  n13 = 0,
  n14 = 1 ,
  n23 = 1,
  n24 = 2,
  n34 = 1,
  n123 = 0,
  n124 = 0,
  n134 = 0 ,
  n234 = 1,
  n1234=0 ,
  category = c("Acide", "Notes T", "Notes 37", "Gel+"),
  fill = c("orange", "red", "green", "blue"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("orange", "red", "green", "blue")
)


venn.plot <- draw.quad.venn(
  area1 =56,
  area2 =154,
  area3 =30,
  area4 =6,
  n12 = 56,
  n13 = 12,
  n14 = 2 ,
  n23 = 23,
  n24 = 5,
  n34 = 6,
  n123 = 12 ,
  n124 = 2,
  n134 = 2,
  n234 = 5,
  n1234= 2,
  category = c("pH<4.8 en jus de soja (56)", "Sac + (154)", "Raf + (30)", "Sta + (6)"),
  fill = c("orange", "red", "green", "blue"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("orange", "red", "green", "blue")
)

jpeg("Fig1a2020.jpeg", width =9, height = 6, units = 'in', res = 1000)
venn.plot <- 
  draw.quad.venn(
  area1 =sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sac_48`+0.3))$specie))),
  area2 =sum(summary(as.factor(filter(Criblage, `NS_48`>(`Raf_48`+0.3))$specie))),
  area3 =sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sta_48`+0.3))$specie))),
  area4 =sum(summary(as.factor(filter(Criblage, `550_48`<6)$specie))),
  n12 = sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sac_48`+0.3),`NS_48`>(`Raf_48`+0.3))$specie))),
  n13 = sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sac_48`+0.3),`NS_48`>(`Sta_48`+0.3))$specie))),
  n14 = sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sac_48`+0.3),`550_48`<6)$specie))),
  n23 = sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sta_48`+0.3),`NS_48`>(`Raf_48`+0.3))$specie))),
  n24 = sum(summary(as.factor(filter(Criblage, `NS_48`>(`Raf_48`+0.3),`550_48`<6)$specie))),
  n34 = sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sta_48`+0.3),`550_48`<6)$specie))),
  n123 =sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sac_48`+0.3),`NS_48`>(`Raf_48`+0.3),`NS_48`>(`Sta_48`+0.3))$specie))),
  n124 =sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sac_48`+0.3),`NS_48`>(`Raf_48`+0.3),`550_48`<6)$specie))),
  n134 =sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sac_48`+0.3),`NS_48`>(`Sta_48`+0.3),`550_48`<6)$specie))),
  n234 =sum(summary(as.factor(filter(Criblage, `NS_48`>(`Raf_48`+0.3),`NS_48`>(`Sta_48`+0.3),`550_48`<6)$specie))),
  n1234=sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sac_48`+0.3)&`NS_48`>(`Raf_48`+0.3)&`NS_48`>(`Sta_48`+0.3)&`550_48`<6)$specie))),
  # category = c("pH<4.7 en jus de soja ()", "Sac + ()", "Raf + ()", "Sta + (29)"),
  fill = c("orange", "red", "green", "blue"),
  lty = "dashed",
  cex.cat = 2,
  cex = c(1.5,1.5,2,1.5,1.5,2,1.5,1.5,1.5,3,2.5,1.5,2,2,2),
  cat.col = c("orange", "red", "green", "blue")
)
dev.off()

draw.triple.venn(
  area1 =sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sac_48`+0.3))$specie))),
  area2 =sum(summary(as.factor(filter(Criblage, `NS_48`>(`Raf_48`+0.3))$specie))),
  area3 =sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sta_48`+0.3))$specie))),
  n12 = sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sac_48`+0.3),`NS_48`>(`Raf_48`+0.3))$specie))),
  n13 = sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sac_48`+0.3),`NS_48`>(`Sta_48`+0.3))$specie))),
  n23 = sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sta_48`+0.3),`NS_48`>(`Raf_48`+0.3))$specie))),
  n123 =sum(summary(as.factor(filter(Criblage, `NS_48`>(`Sac_48`+0.3),`NS_48`>(`Raf_48`+0.3),`NS_48`>(`Sta_48`+0.3))$specie))),
  fill = c("orange", "red", "green"),
  lty = "dashed",
  cex.cat = 2,
  cex = c(3,2.5,1.5,1.5,2,1.5,1.5),
  cat.col = c("orange", "red", "green")
)

