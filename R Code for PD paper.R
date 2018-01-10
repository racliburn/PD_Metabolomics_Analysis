setwd("/Users/Farmer/Desktop")
#Importing the original, unfiltered data.  We won't actually use this data for analysis, but it's useful for some 
#linear models with the demographic table later 
MetabData <- read.table("Metabolomics Excel Data txt.txt", sep="\t", header=T)
#First, we want to do some basic formatting.  This makes our data prettier and compatible with MetaboAnalyst if
#need be.
#We merge the m/z and time columns together
MetabData[,1] <- paste0(MetabData$mz, "/", MetabData$time)
#We title this colum "Sample"
names(MetabData)[1] <- "Sample"
#"We delete the column with only time"
MetabData <- MetabData[,-2]
#we remove the .cdf from the column names.  This is useful because we need to write functions involving the names
#in this data frame and in others, so we need them to be in the same format (the names don't have .cdf in the other
#locations
#names(MetabData) <- gsub(".cdf", "", names(MetabData))
#We create a row to place clinical status in later on.
MetabData <- rbind(rep(NA, times = 38), MetabData)
#

# Mapping file
setwd("/Users/Farmer/Desktop")
#We import the mapping file, which contains file name (e.g. "DW_140904_013"	), sample ID (e.g. "CR0059_01"	),
#and batch (e.g. "1") 
PDmapping <- read.table("sampleID_mapping_pos.txt", sep="\t", header=T)
#We convert the file name to a character so it's in the same class as it is in other locations
PDmapping$FileName <- as.character(PDmapping$FileName)
#We convert the sample ID to a character so it's in the same class as it is in other locations 
#PDmapping$SampleID <- as.character(PDmapping$SampleID)
####
#Make sampleIDs identical in format as in PDdata, rename so the variable names are consistent
#First, we split the sample ID column into two columns, one containing the part of the ID before the underscore
#(e.g. CR0059) and the other containing the part of the ID after the underscore (e.g. 01) 
PDmapping$SampleID <- unlist(lapply(strsplit(PDmapping$SampleID, "_"), function(i) i[1]))
#We rename the sample ID column to CR00 so it's consistent with the demographic table
names(PDmapping)[2] <- "CR00"
##
#We import the table w demographic info
PDdata <- read.table("PD_summary_variables_converted_from_SPSS.txt", sep="\t", header=T)
#We create a subset of the table that excludes alzheimer's patients
Bothdata_subset <- subset(PDdata, caco == "PD" | caco == "control")
#We convert the sample ID to characters so it is in the same data class as everywhere else
Bothdata_subset$CR00 <- as.character(Bothdata_subset$CR00)
#This seems unnecessary but I'm going to leave it in just in case
Bothdata_subset$CR00 <- gsub(" ", "", Bothdata_subset$CR00)

# Merge diagnosis
#Merge the mapping file and the demographic data into one table 
PDmapping1 <- merge(PDmapping, Bothdata_subset[,c(3,5)], by="CR00")
#
# Subset of IDs in the metabolomics dataset
#Now that we've merged the mapping file and demographic table, we can create a demographic table that only has
#the patients in the metabolomics data
PDmapping1 <- PDmapping1[PDmapping1$FileName %in% names(MetabData),]
#Convert the clinical data to character for easier usage    
PDmapping1$caco <- as.character(PDmapping1$caco)
##
setwd("/Users/Farmer/Desktop")
#
#Here we are finally taking all of our organized data and determining clinical status for each patient in the
#metabolomics data
#first, we create a subset of the demographic table with only PD patients and one with only control patients
PDmapping1_subset <- subset(PDmapping1, caco == "PD")
PDmapping1_controlsubset <- subset(PDmapping1, caco == "control")
#We bind a row containing the word "Sample" to the table so that it will match up with the Metabolomic Data table
Fakerow <- list("nope", "Sample", "nope", "nope")
PDmapping1 <- rbind(Fakerow, PDmapping1)
names(MetabData)[1] <- "Sample"
#We eliminate any patients that aren't in the demographic table
MetabData <- MetabData[which(names(MetabData) %in% PDmapping1$FileName)]
#We assign case status to each patient in Metabolomics Data
MetabData[1,] <- ifelse(names(MetabData) %in% PDmapping1_subset$FileName, "PD", "Control")
#x
#Because of the above, the sample column was assigned ctrl, so we remedy that
MetabData[1,1] <- "Label"
########
setwd("/Users/Farmer/Desktop")
#Now we import the data we're actually using for t tests, which already has been filtered by CV, has zeroes
#imputed with half minimum value for each row, and has case status 
MetabDataCV <- read.csv("MetaboAnalyst feature table with mins replaced with half value per row.csv")
#We merge the mz and time columns
MetabDataCV[,1] <- paste0(MetabDataCV$mz, "/", MetabDataCV$time)
#We remove the time only column
MetabDataCV <- MetabDataCV[,-2]
#We add label to the top, so that this is in metaboanalyst format
MetabDataCV[1,1] <- "Label"
#We name the column "Sample"
names(MetabDataCV)[1] <- "Sample"
#Here's the part where we do the glog transformation.  We first create a separate data frame to be transformed.
#Since we want everything to be numeric, we remove the character info 
Metab2 <- as.data.frame(MetabDataCV[2:9072, 2:35])
#We make sure that the data is all numeric
Metab2[] <- lapply(Metab2, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
Metab2[] <- lapply(Metab2, as.numeric)

#Here we define the generalized log function 
glog <- function (i) {
  log2((i+sqrt((i^2)+1))/2)
}

#We apply the generalized log function to the data
Metab2[,1:34] <- lapply (Metab2[,1:34], glog)
#We re-bind the character data 
Metab2 <- cbind (MetabDataCV[-1,1], Metab2)
names(Metab2)[1] <- "Sample"
#We convert the data into character form
MetabDataCV[] <- lapply(MetabDataCV, function(x) {
  if(is.factor(x)) (as.character(x)) else x
})
MetabDataCV[1,] <- lapply(MetabDataCV[1,], as.character)
class(MetabDataCV[1,7])
Metab2 <- rbind(MetabDataCV[1,], Metab2)

#Now, we make our data the log transformed data
MetabDataCV <- Metab2

#We transpose the data, which makes t tests easier
tdata <- t(MetabDataCV)
#We delete the sample names so they don't interfere with the t tests    
tdata <- tdata[-1,]
#We make sure the tdata is expressed as a data frame
tdata <- as.data.frame(tdata)

#We give the column names the sample names so they are excluded from the t tests 
names(tdata) <- MetabDataCV$Sample

#We name the tdata column with Diagnosis "Diagnosis"
names(tdata)[1] <- "Diagnosis"

#We ensure tdata is numeric 
tdata[,-1] <- lapply(tdata[,-1], function(x) {
  if(is.factor(x)) as.numeric((as.character(x))) else x
})
tdata[,-1] <- lapply(tdata[,-1], as.numeric)

#We create a df w results
Ttestresult <- data.frame(MetabDataCV$Sample)
Ttestresult$tstat <- NA
Ttestresult$pval <- NA
Ttestresult$cor <- NA
Ttestresult$pcor <- NA
#We create a separate data frame with only pd patients to be used in the correlation loop 
pdata <- tdata[which(tdata$Diagnosis=="PD"),]

#We run the t-test and correlation loop!
for (i in 1:ncol(tdata)){if (i == 1){next}
  if (i > 1){
    ttest <- t.test(tdata[[i]]~tdata$Diagnosis,var.equal=T)
    Ttestresult$tstat[[i]] <- ttest$statistic
    Ttestresult$pval[[i]] <- ttest$p.value
    Ttestresult$cor[[i]] <- cor(tdata[2639], tdata[[i]], method = "spearman")
    Ttestresult$pcor[[i]] <- cor(pdata[2639], pdata[[i]], method = "spearman")
  }
}

tdata[2639]
pdata[2639]
#We rename the df with results 
MetabResults <- Ttestresult



MetabResults$foldchange <- NA
#We remove the empty first row
MetabResults <- MetabResults[-1,]
#write.csv(MetabResults, file  = "MetabCor.csv")
#write.csv(MetabResults, file = "MetabwithCor.csv")
write.csv(MetabResults, file = "MetabResultswithLogTransformation")
#MetabSignificant <- subset(MetabResults, MetabResults$foldchange > 2 & MetabResults$pval < 0.01) 








for (i in 1:ncol(tdata)){
  if (i == 1){
    next
  }
  if (i > 1){
    MetabResults$foldchange[i-1] <- (mean(tdata[which(tdata$Diagnosis=="PD"),i]))/(mean(tdata[which(tdata$Diagnosis=="Control"),i]))
  }
}

MetabResults$foldchangelog2 <- log2(MetabResults$foldchange)

#getwd ()

#write.csv(MetabResults, file = "MetabResultsCSV.csv")
#We add a column with adjusted p value 
MetabResults$padj <- p.adjust(MetabResults$pval, method = c("fdr"))

#We create a data frame with transposed data used in boxplot creation
boxplotdata <- t(MetabDataCV)
#We create a data frame of boxplotdata (originally it stores as a matrix)  
boxplotdata <- as.data.frame(boxplotdata)
#We give it the same names as tdata
names(boxplotdata) <- MetabDataCV$Sample
#We remove the row with the names 
boxplotdata <- boxplotdata[-1,]
#We give the column with Diagnosis the name of "Diagnosis"w
names(boxplotdata)[1] <- "Diagnosis"
#These libraries are necessary for the next steps
library(lattice)
library(Rcpp)
library(ggplot2)
#We make sure the data that is supposed to be numeric is
boxplotdata[,-1] <- lapply(boxplotdata[,-1], as.character)
boxplotdata[,-1] <- lapply(boxplotdata[,-1], as.numeric)
class(boxplotdata[,5])
#We create a function to make jittery boxplots
jitterrrr  <- function(i) {
  boxplotdata[i] <- lapply(boxplotdata[i], as.character)
  boxplotdata[i] <- lapply(boxplotdata[i], as.numeric)
  boxplotdata$detectionlevel <- boxplotdata[[i]]
  ggplot(boxplotdata, aes(Diagnosis, detectionlevel)) + geom_boxplot(colour="grey50") +geom_jitter(width = 0.1, height = 0.1) +labs(title=names(boxplotdata)[i])
}
jitterrrr(2639)
#class(networktdata[[5]])


#Here we performed a PCA 
#We defined two data frames: one with the tdata numeric data and one with the diagnosis
p <-tdata[,-1] 
p.groups <- tdata[,1]
#We made sure the diagnosis was as factor
groupcol <- factor(p.groups, levels = c("Control", "PD"))
#We perforemd the PCA 
#C
pca <-prcomp(p, center = TRUE, scale. = TRUE)
#We viewed and plotted the PCA
print(pca)
summary(pca)
plot(pca, type = "lines")
#install.packages("devtools")
library("devtools")
#install_github("ggbiplot", "vqv")
library("ggbiplot")
#Using pQckage ggbiplot, we plotted PCA
g <- ggbiplot(pca, var.scale = 1, obs.scale = 1, groups = p.groups, circle = T, ellipse = T, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme (legend.direction = 'horizontal', legend.position = 'top')
print(g)
install.packages("pca3d")
library(pca3d)
#Using the package pca3d, we plotted a 3d pca
pca3d(pca, group=p.groups)
pca2d(pca, group=p.groups)
#We plotted an elbow plot. 
#First, we defined standard deviation
sd <- pca$sdev
#Then, we defined variance as sd^2
var <- sd ^ 2
#Then, we defined a variable as each variance divided by sum of all variances
propvar <- var/sum(var)
#Finally, we created the elbow plots
plot(propvar, xlab = "Principal components", ylab = "Proportion of Variance Explained", type = "b")
plot(cumsum(propvar), xlab ="Principal components", ylab ="Cumulative Proportion of Variance Explained", type = "b")


#PLSDA
#I tried a bunch of different methods to get this to work, but they all failed
#install.packages("muma")
#library("muma")
#oplsda("Metab2v4.csv")
#source("https://bioconductor.org/biocLite.R")
#biocLite("ropls")
#browseVignettes("ropls")
#library("ropls")
#opls1 <- opls(p, p.groups, predI = 1)
#install.packages("mdatools")
#library(mdatools)
l <- as.matrix(tdata)
o <- as.matrix(l[,1])
l <- as.matrix(l[,-1])
l <- as.data.frame(tdata)
o <- as.data.frame(l[,1])
l <- as.data.frame(l[,2:5])
for (i in 1:ncol(l)) {for (j in 1:nrow(l)) {l[j,i] <- as.numeric(l[j,i])}} 

#plsda(l, o)
#plsda1<-plsda(tdata[,2:9072], tdata[,1], ncomp = 1, coeffs.ci ='jk')
plsda1
ci <- plsda1$coeffs$ci
ci <-as.data.frame(ci)
ci <- plsda1$coeffs$p.values
ci <- as.data.frame(ci)
plsda1
xloading <- plsda1$xloadings
pvalss <- plsda1$coeffs$p.values
pvalss<- as.data.frame(pvalss)
plsda1$ncomp
#head(tdata)
#edit(tdata)
#names(tdata)[2] <- "Label",
#tdata2<-tdata[complete.cases(tdata),]
#edit(MetabData)
#MetabData[2,2:35] <- unlist(MetabData[2,2:35])
#t.test(unlist(MetabData[which(names(MetabData) %in% PDmapping1_subset$FileName)]) ~ unlist(MetabData[which(names(MetabData) %in% PDmapping1_controlsubset$FileName)]), data=MetabData[2,2:35], var.equal=TRUE, conf.level=0.95)

#MetabData[-1,-1] <- sapply(MetabData[-1,-1], as.numeric)
#class(MetabData[3, 10])
install.packages("xmsannotator")


sewtd("Users:Farmer:Downloads")
setwd("/Users/Farmer/Downloads")



#Here we are using the full demographic data to create linear models for different metabolites against some data
#First, we import the full demographic tabl3
FullPDData <- read.csv("Master PD Control Dataset.csv")

#Fulldata_subset <- subset(FullPDdata, caco == "PD" | caco == "control")
#We create a subset that excludes Alzheimer's patients
Fulldata_subset <- subset (FullPDData, caco == "PD"| caco=="control")
#We make sure it is all in the right form (as character, etc)
Fulldata_subset$CR00 <- as.character(Fulldata_subset$CR00)
Fulldata_subset$CR00 <- gsub(" ", "", Fulldata_subset$CR00)
# Merge diagnosis
#We merge the demographic table with the mapping as we saw earlier
PDmappingfull <- merge(PDmapping, Fulldata_subset[,c(2,3,5, 19:27)], by="CR00")

# Subset of IDs in the metabolomics dataset
#We make sure that we're only looking at patients w metabdata  
PDmappingfull<- PDmappingfull[(PDmappingfull$FileName %in% names(MetabData)),]
#We make sure the data is in the right form (as.character, remove .cdf)
PDmappingfull$caco <- as.character(PDmappingfull$caco)
gsubcdf <- function (i ) {gsub(".cdf","",i)} 
rownames(tdata) <- lapply(rownames(tdata),gsubcdf)

#We add a column with the metabolic data DOPAL
PDmappingfull <- merge(PDmappingfull, tdata[1170], by.x="FileName", by.y = "row.names")



#names(PDmappingfull[15]) <- "DOPAL"
#We name the column "DOPAL"
names(PDmappingfull)[15 ] <- "DOPAL"
#lm(PDmappingfull[15], PDmappingfull[5])
#PDmappingfull[15] <- lapply(PDmappingfull[15], as.numeric)
#PDmappingfull[5] <- lapply(PDmappingfull[5], as.character)
#PDmappingfull[5] <- lapply(PDmappingfull[6], as.numeric) 
#PDmappingfull[15] <- lapply(PDmappingfull[15], as.character)
#PDmappingfull[15] <- lapply(PDmappingfull[15], as.numeric)

#We create a linear model and examine the results
lmdopal <- lm(PDmappingfull$levodopa_eq ~ PDmappingfull$DOPAL)
summary(lmdopal)
#We calculate the correlation coefficient
cor(PDmappingfull[15], PDmappingfull[6])
#We create a scatterplot
plot(PDmappingfull$levodopa_eq, PDmappingfull$DOPAL)

#REPEAT for different metabolites  ( Harmalol, Genistein, Cotinine, D)

PDmappingfull <- merge(PDmappingfull, tdata[2865], by.x="FileName", by.y ="row.names" )
names(PDmappingfull)[16] <- "Harmalol" 
lmharmalol <- lm(PDmappingfull$levodopa_eq ~ PDmappingfull$Harmalol)
summary(lmharmalol)
plot(PDmappingfull$levodopa_eq, PDmappingfull$Harmalol)

PDmappingfull <- merge(PDmappingfull, tdata[2909], by.x = "FileName", by.y = "row.names")
names(PDmappingfull)[17] <- "Genistein"
lmgenistien <- lm(PDmappingfull$levodopa_eq ~ PDmappingfull$Genistein)
summary(lmgenistien)
plot(PDmappingfull$levodopa_eq, PDmappingfull$Genistein)

PDmappingfull <- merge (PDmappingfull, tdata[1833], by.x = "FileName", by.y = "row.names")
names(PDmappingfull)[18] <- "Cot"
lmcotharm <- lm(PDmappingfull$Harmalol ~ PDmappingfull$Cot)
summary(lmcotharm)
plot(PDmappingfull$Harmalol, PDmappingfull$Cot)

tdata[1170]
tdata[2865]
write.csv(MetabResults, file = "PDMetabolomicsResults.csv")

