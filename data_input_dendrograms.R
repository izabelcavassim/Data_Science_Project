######################################################################
# Project: Data Science in Bioinformatics - Visualization and Analysis 
# Student: Maria Izabel Cavassim Alves
# Professor: Thomas Mailund
# Data Input and cleaning
# Clustering dendrogram of samples based on their Euclidean distance
######################################################################

#install.packages("WGCNA", repos=c("http://rstudio.org/_WGCNA", "http://cran.rstudio.com"))
#Not in CRAN repository:
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("GO.db", "impute", "preprocessCore"))
require(WGCNA)

#Editing the phenotype data from the article: Lorwy et al., 2013
getwd()
setwd("C:/Users/Maria/Desktop/MASTER/DATA_SCIENCE/Data/Coexpression/")

pheno = read.table("probes_log2_normalized.txt", header=TRUE)
row.names(pheno) = pheno$Line.Marker
pheno = pheno[,-1]
#the genes ending with a even number are from dry environment, 
#ending with odd numbers are from the wet environment
#View(pheno)

#separating the wet sample
pheno_wet = pheno[,seq(1, ncol(pheno),by=2)] 

#separating the dry sample
pheno_dry = pheno[,seq(0, ncol(pheno),by=2)]

#changing the genotypes names for the same names_code written in genotype data
#They are listed in the same order as the phenotype data
RIL_names = read.csv("names_ril_cor.csv", header=FALSE, sep = ";" )
head(RIL_names)
RIL_names = unique(RIL_names[,2])
length(RIL_names) #it is right! checked

colnames(pheno_dry) = RIL_names
colnames(pheno_wet) = RIL_names
#View(pheno_dry) #checking if everything is ok!

pheno_dry = as.data.frame(t(pheno_dry))
pheno_wet = as.data.frame(t(pheno_wet))

#check for genes and samples with too many missing values:
gsg = goodSamplesGenes(pheno_wet, verbose = 3);
gsg$allOK

sampleTree1 = hclust(dist(pheno_wet), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
pdf(file = "Wet Clustering.pdf");
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree1, main = "Wet sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree1, cutHeight = 56, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = pheno_wet[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

save(datExpr, file = "Wet-dataInput.RData")
########################################################
# Analysing the pheno dry data
########################################################
gsg2 = goodSamplesGenes(pheno_dry, verbose = 3);
gsg2$allOK

sampleTree2 = hclust(dist(pheno_dry), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
pdf(file = "Dry Clustering.pdf");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree2, main = "Dry sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Determine cluster under the line
clust = cutreeStatic(sampleTree2, cutHeight = 55, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = pheno_dry[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

save(datExpr, file = "Dry-dataInput.RData")