MonFichierReport<-choose.files(default = "", caption = "Selectionner le fichier 'report.tsv'",multi = FALSE) 

library("stringr")

#############Calcul taille Protéine from Fasta
#trouve fasta
Mylog<-read.delim(gsub("report.tsv", "report.log.txt", MonFichierReport),sep="\t")
MyFasta<-Mylog[str_detect(Mylog[,1],".fasta")==TRUE ,][1]
MyFasta<-str_sub(MyFasta,53,nchar(MyFasta))

#import fasta et retire les lignes NA
Fasta<-read.delim(MyFasta,header = FALSE)
Fasta<-Fasta[!is.na(Fasta),]
Fasta<-as.data.frame(Fasta)

if (!require("stringr")) {install.packages("stringr", dependencies = TRUE)}

ExtractAccession <- function(x){
Accession<-strsplit(x, split = "\\|") [[1]][2]
return(Accession)
}

#calcul length et sequence
MyListeAccession<-c()
MyLenght<-c()
MySequence<-c()
for(i in 1:nrow(Fasta)){
    if(substr(Fasta[i,1], 1, 1)==">"){
      MyListeAccession<-c(MyListeAccession,ExtractAccession(Fasta[i,1]))
      L<-0
      S<-""
    }else{
      L<-L+nchar(Fasta[i,1]) 
      S<-paste0(S,Fasta[i,1])
      if(substr(Fasta[i+1,1], 1, 1)==">"||i==nrow(Fasta)){
        MyLenght<-c(MyLenght,L)
        MySequence<-c(MySequence,S)
      }
    }
}

FastaProtLenght <- data.frame (Accession  = MyListeAccession,ProtLenght = MyLenght,ProtSequence=MySequence)

rm(Fasta,Mylog,MyFasta)


##########Calcul nombre d'Acide aminé détecté par prot
x<-read.table(MonFichierReport,
              # colClasses = c("character",rep("NULL",2),"character",rep("NULL",10),"character",rep("NULL",40)),      
              header = TRUE,
              sep ="\t")
x<-x[,c("File.Name","Protein.Ids","Stripped.Sequence")]


#Extrait première prot du Protein.Ids
ExtractFirstAccession <- function(x){
  FirstAccession<-strsplit(x, split = ";") [[1]][1]
  return(FirstAccession)
}

x$FirstAccession<-as.character(lapply(x$Protein.Ids,ExtractFirstAccession))

#supprime doublon
if (!require("dplyr")) {install.packages("dplyr", dependencies = TRUE)}
library(dplyr)
x <- distinct(x)

#Merge report et FastaProtLenght
x<-merge(x,FastaProtLenght,all.X=TRUE,all.y = FALSE, by.x="FirstAccession",by.y="Accession")


#liste des positions des AA de chaque peptide
PositionPep<-function(i){
  P<-str_locate(pattern = x$Stripped.Sequence[i], x$ProtSequence[i])
  listePosition<-paste(P[1,1]:P[1,2],collapse =";")
  return(listePosition)
}

x$PosPep<-as.character(lapply(1:nrow(x),PositionPep))


#regroupe les PosPep par firstProtienGroup+FileName
Mypaste<-function(x){
  paste(x,collapse =";")
}
x<-aggregate(data=x,PosPep~FirstAccession+File.Name,"Mypaste")

#élimine doublon dans PosPep
for(i in 1:nrow(x)){
  x$Cov[i]<-length(unique(strsplit(x$PosPep[i], split = ";")[[1]]))
}


#Calcul%Cov
x<-merge(x,FastaProtLenght,all.X=TRUE,all.y = FALSE, by.x="FirstAccession",by.y="Accession")
x$PcCov<-100*x$Cov/x$ProtLenght
x<-x[,c(1,2,7)]

if (!require("reshape2")) {install.packages("reshape2", dependencies = TRUE)}
library(reshape2)
x<-reshape(data=x,idvar="FirstAccession",v.names = "PcCov", timevar = "File.Name", direction="wide")

#Output
Myoutput<-read.delim(gsub("report", "report.pg_matrix", MonFichierReport)  )
Myoutput$FirstAccession<-as.character(lapply(Myoutput$Protein.Ids,ExtractFirstAccession))
Myoutput<-merge(Myoutput[,c(1:4,ncol(Myoutput))],x,all.X=TRUE,all.y = FALSE, by.x="FirstAccession",by.y="FirstAccession")
write.table(Myoutput,gsub("report", "report.pg_matrix_PcCoverage", MonFichierReport),row.names = FALSE,sep = "\t")

for(i in 6:ncol(Myoutput)){
  hist(Myoutput[,i],main=paste0("colonne ",i,", moyenne=",mean(Myoutput[,i],na.rm = TRUE)))
}

cat("Le fichier PcCoverage.tsv est arrivé à destination.\nMerci d'avoir choisi Proteom'IC R-line ,\net à bientot.")
