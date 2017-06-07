# FOS Course 2017
# R script belonging to practical on Epigenomics, June 22nd 2017

# Packages required
library(ggplot2)

#set location of files
setwd("change to the location of the files")
#setwd("/Users/roderick/Documents/Onderwijs/FOS/2017/")

#############
# Q1

# Load data
TissueData <- read.table("Tissue data subset.csv", sep=";", header=T)

# Look at data
head(TissueData) # First 6 lines of the data
str(TissueData) # Characteristics of the data

# Q1a
unique(TissueData$Tissue) # What groups are in the vector?

# Q1b
ggplot(TissueData, aes(x=Tissue, y=cg21507095, fill=Tissue))+ # Assign variables to x, y, and fill
	geom_boxplot()+ # Make a boxplot
	scale_fill_manual(values=c("#F9A23F","#009AC7","#8B1A4F")) # Give custom colors to the boxes using HEX colors
  
#Mean and SD
# Data = data to calculate statistic over, INDICES = grouping variable, FUN = function to be used
by(data = TissueData$cg21507095, INDICES = TissueData$Tissue, FUN=mean)
by(data = TissueData$cg21507095, INDICES = TissueData$Tissue, FUN=sd)

# Q1c
#Run the anova for one CpG
anova(lm(cg21507095~Tissue, data=TissueData)) 
  
#

#############
# Q2
TissueData22 <- read.table("Tissue data chromosome 22.csv", sep=";", header=T)
head(TissueData22)

# Q2a
nrow(TissueData22) #Number of rows

# Q2b
ggplot(TissueData22, aes(x=Blood))+
  geom_histogram(bins = 200) # Histogram

# Q2d
ggplot(TissueData22, aes(x=Gene_centric, y=Blood, fill=Gene_centric))+
  geom_boxplot()+
  scale_fill_manual(values=c("#D12633","#F9A23F","#4B8E55","#009AC7","#8B1A4F"))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) # Rotate axis labels

# Q2e
ggplot(TissueData22, aes(x=CGI_centric, y=Blood, fill=CGI_centric))+
  geom_boxplot()+
  scale_fill_manual(values=c("#D12633","#F9A23F","#4B8E55","#009AC7","#8B1A4F"))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# Q2f
load("TissueData chr22.RData") # Contains two files samplesheet, TissueDataChr22
head(samplesheet)
dim(samplesheet) # Number of rows and number of columns
head(TissueDataChr22)
dim(TissueDataChr22)

# Q2g

#Function to return the Pvalue for each of the CpG sites
getPvalue <- function(row, data, samplesheet)
{
  fit <- anova(lm(data[row,]~samplesheet$Tissue)) # Run the model
  out <- data.frame(CG = row, Fvalue = fit$`F value`[1], Pvalue = fit$`Pr(>F)`[1]) # Make a data.frame for the output
  as.matrix(out)
}

#Run for one CpG site
getPvalue('cg00017461', data=TissueDataChr22, samplesheet=samplesheet)

# Q2h
#Which CpGs are there
cpgs <- rownames(TissueDataChr22)
head(cpgs)
# Now run the model for all CpGs on chromosome 22
res <- lapply(cpgs, getPvalue, data=TissueDataChr22, samplesheet=samplesheet) # Returns a list of results per CpG
head(res)
ResultsChromosome22 <- do.call(rbind, res) #Combine to a table
head(ResultsChromosome22)

#############
# Q3

# Q3d
table(TissueData22$Pval.corrected)

# Q3e
TissueData22sign <- read.table("Tissue data chromosome 22 significant.csv", sep=";", header=T)
head(TissueData22sign)
table3e <- table(TissueData22sign$Gene_centric, TissueData22sign$CGI_centric)
table3e

# Q3f
table3f <- (table3e / 1805)*100
round(table3f,1)

#############
# Q4
TissueData22signOrdered <- TissueData22sign[order(TissueData22sign$Pval),]
head(TissueData22signOrdered)



