setwd("../data")

#Installer packages pr manip données
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler", version = "3.8")

#GOstats package
source("https://bioconductor.org/biocLite.R")
biocLite("GOstats")
library(GOstats)

###COG###
#Extraire COG
setwd("../data/COG")
library(readr)
categories <- read_csv("cog2cat42163_10-nov-2018.xls.csv")
db <- read_excel("GO_COG.xlsx", col_names = FALSE)
View(db)
counts=table(na.omit(db$COG))
COG=data.frame(counts)
colnames(COG)=c('COG ID','Freq')
#Merge les tableaux
merged=merge(COG,categories)
PT=aggregate(merged$Freq~merged$Catergory,merged,sum)
PT$legend=c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','Z')
colnames(PT)=c('levels','Frequency','legend')
#COG graph
library(ggplot2)
p <- ggplot(data=PT, aes(x=levels, y=Frequency,fill=legend))+geom_bar(stat="identity", position=position_dodge(), colour='Seashell')
p + coord_flip() + scale_x_discrete(limits = rev(PT$levels)) + guides (fill = guide_legend(ncol = 1))+ylab("Frequency")+xlab('Factor class')+ggtitle("COG Representation by functional category")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_text(aes(label=Frequency), position=position_dodge(width=0.9), hjust=-0.25)
#COG graph en %
merged$Percent=(merged$Freq/sum(merged$Freq)*100)
datp <- data.frame(FunctionClass = factor(merged$`Catergory Code`), levels=merged$`Catergory Code`,legend = merged$Catergory,Percentage=merged$Percent)
library(ggplot2)
pp <- ggplot(data=datp, aes(x=FunctionClass, y=Percentage, fill=legend))+geom_bar(stat="identity", position=position_dodge(), colour="seashell")
pp + coord_flip() + guides (fill = guide_legend(ncol = 1))+xlab("Frequency")+ylab('Factor class')+ggtitle("COG Percentage by functional category")
#Tableau de %
Pt=aggregate(Percent~`Catergory Code`,merged,sum)
write.csv(Pt,file='Percent table.csv')

#COG/traitement
library(readr)
categories <- read_csv("COG/BET/cog2cat42163_10-nov-2018.xls.csv")
library(readxl)
BET_COG <- read_excel("COG/BET/BET_COG.xlsx")
STD_COG <- read_excel("COG/STD/STD_COG.xlsx")
counts_BET=table(na.omit(BET_COG$`COG ID`))
counts_STD=table(na.omit(STD_COG$`COG ID`))
COG_BET=data.frame(counts_BET)
COG_STD=data.frame(counts_STD)
colnames(COG_BET)=c('COG ID','Freq')
colnames(COG_STD)=c('COG ID','Freq')
merged_BET=merge(COG_BET,categories)
merged_STD=merge(COG_STD,categories)
PT_BET=aggregate(merged_BET$Freq~merged_BET$Catergory,merged_BET,sum)
PT_STD=aggregate(merged_STD$Freq~merged_STD$Catergory,merged_STD,sum)
colnames(PT_BET)=c('levels','Frequency')
colnames(PT_STD)=c('levels','Frequency')
PT_BET$legend=c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','Z')
PT_STD$legend=c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','Z')
library(ggplot2)
p_BET <- ggplot(data=PT_BET, aes(x=levels, y=Frequency,fill=legend))+geom_bar(stat="identity", position=position_dodge(), colour='Seashell')
p_BET + coord_flip() + scale_x_discrete(limits = rev(PT_BET$levels)) + guides (fill = guide_legend(ncol = 1))+ylab("Frequency")+xlab('Factor class')+ggtitle("COG Representation by functional category for betuline treatment") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ geom_text(aes(label=Frequency), position=position_dodge(width=0.9), hjust=-0.25)
p_STD <- ggplot(data=PT_STD, aes(x=levels, y=Frequency,fill=legend))+geom_bar(stat="identity", position=position_dodge(), colour='Seashell')
p_STD + coord_flip() + scale_x_discrete(limits = rev(PT_STD$levels)) + guides (fill = guide_legend(ncol = 1))+ylab("Frequency")+xlab('Factor class')+ggtitle("COG Representation by functional category for standard treatment") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ geom_text(aes(label=Frequency), position=position_dodge(width=0.9), hjust=-0.25)
Pt_BET=aggregate(Freq~`Catergory Code`,merged_BET,sum)
Pt_STD=aggregate(Freq~`Catergory Code`,merged_STD,sum)
write.csv(Pt_BET,file='Percent table_BET.csv')
write.csv(Pt_STD,file='Percent table_STD.csv')

###Graph abondance pfam###
#Load data
GO_pfam <- read_excel("Total/GO_pfam.xlsx", col_names = FALSE)
categories <- read_csv("Total/pfamlist48269_20-nov-2018.xls.csv")
#Faire disparaitre les NA's et merger
counts=as.data.frame(table(na.omit(GO_pfam$GO_Pfam)))
colnames(counts)=c('Pfam ID','Freq')
merged=merge(counts,categories)
#Ordoner pfam les plus de frequents et isoler 20
GO_sort=merged[order(-merged$Freq),]
GO_short=as.data.frame(GO_sort$Freq[1:20])
#Pfam en rownames
names <- paste(GO_sort$`Protein name`[1:20], sep=":")
rownames(GO_short)=names
colnames(GO_short)=c('Freq')
dat <- data.frame(FunctionClass = factor(row.names(GO_short)),legend=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'), levels=row.names(GO_short),Frequency=GO_short$Freq)
#Plot graph
library(ggplot2)
p <- ggplot(data=dat, aes(x=reorder(levels,-Frequency), y=Frequency))+geom_bar(stat="identity", position=position_dodge(), colour='Seashell')
p + coord_flip() + guides (fill = guide_legend(ncol = 1))+ylab("Frequency")+xlab('Protein family name')+ggtitle("Protein families representation by functional category")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_text(aes(label=Frequency), position=position_dodge(width=0.9), hjust=-0.25)

#Pfam/traitement
setwd("../data/pfam")
library(readr)
categories <- read_csv("Total/pfamlist48269_20-nov-2018.xls.csv")
library(readxl)
STD_GO_pfam <- read_excel("STD/STD_pfam.xlsx")
BET_GO_pfam <- read_excel("BET/BET_pfam.xlsx")
BET_counts=as.data.frame(table(na.omit(BET_GO_pfam$X__2)))
STD_counts=as.data.frame(table(na.omit(STD_GO_pfam$X__2)))
colnames(BET_counts)=c('Pfam ID','Freq')
colnames(STD_counts)=c('Pfam ID','Freq')
BET_merged=merge(BET_counts,categories)
STD_merged=merge(STD_counts,categories)
BET_GO_sort=BET_merged[order(-BET_merged$Freq),]
STD_GO_sort=STD_merged[order(-STD_merged$Freq),]
BET_GO_short=as.data.frame(BET_GO_sort$Freq[1:20])
STD_GO_short=as.data.frame(STD_GO_sort$Freq[1:20])
BET_names <- paste(BET_GO_sort$`Pfam ID`[1:20], BET_GO_sort$`Pfam Name`[1:20], sep=":")
rownames(BET_GO_short)=BET_names
STD_names <- paste(STD_GO_sort$`Pfam ID`[1:20], STD_GO_sort$`Pfam Name`[1:20], sep=":")
rownames(STD_GO_short)=STD_names
colnames(BET_GO_short)=c('Freq')
colnames(STD_GO_short)=c('Freq')
dat_BET <- data.frame(FunctionClass = factor(row.names(BET_GO_short)),legend=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'), levels=row.names(BET_GO_short),Frequency=BET_GO_short$Freq)
dat_STD <- data.frame(FunctionClass = factor(row.names(STD_GO_short)),legend=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'), levels=row.names(STD_GO_short),Frequency=STD_GO_short$Freq)
library(ggplot2)
p <- ggplot(data=dat_BET, aes(x=reorder(levels,-Frequency), y=Frequency,fill=legend))+geom_bar(stat="identity", position=position_dodge(), colour='Seashell')
p + coord_flip() + guides (fill = guide_legend(ncol = 1))+ylab("Frequency")+xlab('Factor class')+ggtitle("Protein families representation by functional category for betuline treatment")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_text(aes(label=Frequency), position=position_dodge(width=0.9), hjust=-0.25)
p <- ggplot(data=dat_STD, aes(x=reorder(levels,-Frequency), y=Frequency,fill=legend))+geom_bar(stat="identity", position=position_dodge(), colour='Seashell')
p + coord_flip() + guides (fill = guide_legend(ncol = 1))+ylab("Frequency")+xlab('Factor class')+ggtitle("Protein families representation by functional category for standard treatment")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_text(aes(label=Frequency), position=position_dodge(width=0.9), hjust=-0.25)

###KEGG###
setwd("../data")
KEGG <- read_csv("KEGG/KEGG.csv")
categories <- read_csv("KEGG/KO_categories.xls.csv")
#Faire disparaitre les NA's et merger
library("dplyr")
KEGG_simple=na.omit(select(KEGG,"EC Codes","KEGG KO"))
counts_KEGG_EC=as.data.frame(table(KEGG_simple))

counts_EC=as.data.frame(table(na.omit(KEGG$"EC Codes")))
counts_KEGG=as.data.frame(table(na.omit(KEGG$"KEGG KO")))

colnames(counts)=c('KO','Freq')
colnames(categories)=c('KO','Protein name')
merged=merge(counts,categories)
#Ordoner pfam les plus de frequents et isoler 20
GO_sort=merged[order(-merged$Freq),]
GO_short=as.data.frame(GO_sort$Freq[1:20])
#Pfam en rownames
names <- paste(GO_sort$`Protein name`[1:20], sep=":")
rownames(GO_short)=names
colnames(GO_short)=c('Freq')
dat <- data.frame(FunctionClass = factor(row.names(GO_short)),legend=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'), levels=row.names(GO_short),Frequency=GO_short$Freq)
#Plot graph
library(ggplot2)
p <- ggplot(data=dat, aes(x=reorder(levels,-Frequency), y=Frequency))+geom_bar(stat="identity", position=position_dodge(), colour='Seashell')
p + coord_flip() + guides (fill = guide_legend(ncol = 1))+ylab("Frequency")+xlab('Protein name')+ggtitle("KEGG abundance representation by category")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_text(aes(label=Frequency), position=position_dodge(width=0.9), hjust=-0.25)
