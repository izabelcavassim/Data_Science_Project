######################################################################
# Project: Data Science in Bioinformatics - Visualization and Analysis 
# Student: Maria Izabel Cavassim Alves
# Professor: Thomas Mailund
# Data Input and cleaning
# Modules construction for wet and dry envyronment
######################################################################

#install.packages("WGCNA", repos=c("http://rstudio.org/_WGCNA", "http://cran.rstudio.com"))
#Not in CRAN repository:
#source("https://bioconductor.org/biocLite.R")
#biocLite(c("GO.db", "impute", "preprocessCore"))
require(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()

#lnames = load('C:/Users/Maria/Documents/GitHub/Data_Science_Project/Wet-dataInput.RData')
lnames = load('Wet-dataInput.RData')
lnames
#############################################################################################################################
#Constructing a weighted gene network entails the choice of the soft thresholding power ?? to which co-expression
#similarity is raised to calculate adjacency [1]. The authors of [1] have proposed to choose the soft thresholding power
#based on the criterion of approximate scale-free topology. We refer the reader to that work for more details; here
#we illustrate the use of the function pickSoftThreshold that performs the analysis of network topology and aids the
#user in choosing a proper soft-thresholding power.
#############################################################################################################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


bwnet_wet = blockwiseModules(datExpr, maxBlockSize = 26000,
                              power = 6, TOMType = "unsigned", minModuleSize = 30,
                              reassignThreshold = 0, mergeCutHeight = 0.25,
                              numericLabels = TRUE,
                              saveTOMs = TRUE,
                              saveTOMFileBase = "wetTOM-blockwise",
                              verbose = 3)

save(bwnet_wet, file='modules_wet.RData')