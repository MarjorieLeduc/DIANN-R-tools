MonDossier<-choose.dir(default = "", caption = "Select_folder")  
setwd(MonDossier)

#Import data
MyReport<-read.table("report.tsv",header=TRUE,sep="\t")

#import library
if (!require("reshape2", character.only = TRUE)) {
  install.packages("reshape2", dependencies = TRUE)
  library("reshape2", character.only = TRUE)
}
library(reshape2)

#calcul nb pep per protein group
x<-dcast(MyReport, Protein.Group ~ File.Name,value.var="PG.MaxLFQ",fun.aggregate=length)
colnames(x)<-c("Protein.Group",paste("Nb precurseur",colnames(x)[2:ncol(x)]))

#export matrice
write.table(x, file = paste0("report.pg_matrix_NbPrecurseur.tsv"),sep="\t",row.names = FALSE)

