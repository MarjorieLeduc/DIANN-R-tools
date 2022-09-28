if (!require("stringr")) {install.packages("stringr", dependencies = TRUE)}
library(stringr)

MonFichier<-choose.files(default = "", caption = "Selectionner le fichier 'report.txt'",multi = FALSE) 

df <- read.table(MonFichier,
                 comment.char = "", header = T, sep = "\t", na.strings = "NA",stringsAsFactors = F, quote = "")

AAMW<-matrix(nrow=0,ncol=2,dimnames = list(NULL,c("Letter","Monoisotopique")))
AAMW<-rbind(AAMW,c("A" ,71.03711))
AAMW<-rbind(AAMW,c("R" ,156.10111))
AAMW<-rbind(AAMW,c("N" ,114.04293))
AAMW<-rbind(AAMW,c("D" ,115.02694))
AAMW<-rbind(AAMW,c("C" ,103.00919))
AAMW<-rbind(AAMW,c("E" ,129.04259))
AAMW<-rbind(AAMW,c("Q" ,128.05858))
AAMW<-rbind(AAMW,c("G" ,57.02146))
AAMW<-rbind(AAMW,c("H" ,137.05891))
AAMW<-rbind(AAMW,c("I" ,113.08406))
AAMW<-rbind(AAMW,c("L" ,113.08406))
AAMW<-rbind(AAMW,c("K" ,128.09496))
AAMW<-rbind(AAMW,c("M" ,131.04049))
AAMW<-rbind(AAMW,c("F" ,147.06841))
AAMW<-rbind(AAMW,c("P" ,97.05276))
AAMW<-rbind(AAMW,c("S" ,87.03203))
AAMW<-rbind(AAMW,c("T" ,101.04768))
AAMW<-rbind(AAMW,c("W" ,186.07931))
AAMW<-rbind(AAMW,c("Y" ,163.06333))
AAMW<-rbind(AAMW,c("V" ,99.06841))
AAMW<-rbind(AAMW,c("(UniMod:4)" ,57.021464))
AAMW<-rbind(AAMW,c("(UniMod:35)" ,15.994915))
AAMW<-as.data.frame(AAMW)

MyMW<-df$Modified.Sequence

for(i in 1:nrow(AAMW)){
  MyAA<-AAMW$Letter[i]
  MyAAMW<-AAMW$Monoisotopique[i]
  
  MyMW<-gsub(MyAA, paste0(MyAAMW,";"), MyMW)
}


MyMW<-str_split(MyMW,";")
MyMW<-lapply(MyMW,as.numeric)
MyMW<-lapply(MyMW,na.exclude)
MyMW<-lapply(MyMW, sum)
MyMW<-unlist(MyMW)

H<-1.007825
O<-15.994915
H2O<-H+H+O

MyMW<-MyMW+H2O
df$MW<-MyMW
df$mz<-(df$MW+df$Precursor.Charge)/df$Precursor.Charge

hist(df$mz,breaks=40)
min(df$mz)#184.4119
max(df$mz)#1199.651

write.table(df,paste0(dirname(MonFichier),"/report_mz.txt"),quote = FALSE,row.names = FALSE,sep = "\t")
