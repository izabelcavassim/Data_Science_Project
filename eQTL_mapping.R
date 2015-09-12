######################################################################
# Project: Data Science in Bioinformatics - Visualization and Analysis 
# Student: Maria Izabel Cavassim Alves
# Professor: Thomas Mailund
# Using the software r-qtl to eQTL mapping
######################################################################

#Instaling r/qtl package
#install.packages('qtl')
require(qtl)

#Editing the phenotype data from the article: Lorwy et al., 2013
getwd()
setwd("C:/Users/Maria/Desktop/MASTER/DATA_SCIENCE/Data/eQTL mapping")

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

#saving the data in a .csv file
write.csv(t(pheno_dry), file = "pheno_dry.csv")
write.csv(t(pheno_wet), file = "pheno_wet.csv")

#treating the genotype file
genotypes <-read.csv("tpc115352SupplementalDS5.csv", header=TRUE,sep=";") 
linkage_map <- genotypes[1:2,]

#Converting letters to numbers
#That is not necessary, just wanted to see how to change
genotypes[3:5,3:5]
#Converting A to 1
genotypes1 <-sapply(genotypes,as.character) 
genotypes1 <-data.frame(gsub("A","1", genotypes1))
#Converting H to 0
genotypes1 <-sapply(genotypes1,as.character) 
genotypes1 <-data.frame(gsub("H","0", genotypes1))

genotypes1[3:5,3:5] 

genotypes_match = genotypes[match(RIL_names, genotypes$Line.Marker),]

#putting the information of the linkage map:
genotypes_match = rbind(linkage_map, genotypes_match)

write.csv(genotypes_match, file = "genotypes_match.csv")

#############################################################################
# eQTL analysis for the WET environment separated using r/qtl package
#############################################################################
#Downloading the linkage map and genotype of KAS-1 TSU-1
#Loading the phenotype expression data (wet and dry)
require(qtl)
genotype_map = read.cross("csvs", dir=".", genfile = "genotypes_match.csv", 
                          phefile = "pheno_wet.csv", genotypes =c("A","H"), header=FALSE, sep=',', dec=',' )

summary(genotype_map)
nind(genotype_map) #numer of individuals
nchr(genotype_map) #number of chromossomes
totmar(genotype_map) #total markers
nmar(genotype_map) #number of markers per chromossome
nphe(genotype_map) #number of phenotypes

#plotting the summary of the the data
#pdf("Linkage map")
plot(genotype_map, pheno.col = 3) #Analizing just the first phenotype sample/column
#dev.off()

#plotting the missing values of linkage map
plotMissing(genotype_map) 
#The number of markers is relatively low, thats why there are so many 'missing values

#pdf("Genetic Map")
plotMap(genotype_map)
#dev.off()

#Estimating the recombination fraction. Uses the hidden Markov model method to calculate the 
#probabilities of the true underlying genotypes given the observed multipoint marker data,
#with possible allowance for genotyping errors.

#The Haley Knott (HK) regression method continues to be a popular approximation to 
# standard interval mapping (IM) of quantitative trait loci (QTL) in experimental crosses.
# Converting genetic distances into recombination fractions using Kosambi mapping function

genotype_map = calc.genoprob(genotype_map, step=1, map.function="kosambi", error.prob = 0.01, stepwidth= "fixed")
genotype_map = convert2riself(genotype_map)

pheno_wetLOD1=scanone(genotype_map, pheno.col=2:5000, model=c("normal"), method=c("hk"))
pheno_wetLOD2=scanone(genotype_map, pheno.col=5001:10000, model=c("normal"), method=c("hk"))
pheno_wetLOD3=scanone(genotype_map, pheno.col=10001:15000, model=c("normal"), method=c("hk"))
pheno_wetLOD4=scanone(genotype_map, pheno.col=15001:20000, model=c("normal"), method=c("hk"))
pheno_wetLOD5=scanone(genotype_map, pheno.col=20001:25663, model=c("normal"), method=c("hk"))
transcriptLODall=cbind(pheno_wetLOD1, pheno_wetLOD2, pheno_wetLOD3, pheno_wetLOD4, pheno_wetLOD5)


# function scanone is used to perform a permutation test to get a genome-wide LOD significance threshold. 
set.seed(123456)
maxlod1=scanone(genotype_map, pheno.col=2:5000, model=c("normal"), method=c("hk"), n.perm=1000)
tabbycol_5000<-summary(pheno_wetLOD1, format="tabByCol", pvalues=TRUE, perm=maxlod1, ci.function="lodint")
write.table(tabbycol_5000, file = "Tab_by_5000_Mckay_sum_SNP.txt", quote=FALSE, row.name=FALSE)


set.seed(123456)
maxlod2=scanone(genotype_map, pheno.col=5001:10000, model=c("normal"), method=c("hk"), n.perm=1000)
tabbycol_5001<-summary(pheno_wetLOD2, format="tabByCol", pvalues=TRUE, perm=maxlod2, ci.function="lodint")
write.table(tabbycol_5001, file = "Tab_by_5001_Mckay_sum_SNP.txt", quote=FALSE, row.name=FALSE)


set.seed(123456)
maxlod3=scanone(genotype_map, pheno.col=10001:15000, model=c("normal"), method=c("hk"), n.perm=1000)
tabbycol_10001<-summary(pheno_wetLOD3, format="tabByCol", pvalues=TRUE, perm=maxlod3, ci.function="lodint")
write.table(tabbycol_10001, file = "Tab_by_10001_Mckay_sum_SNP.txt", quote=FALSE, row.name=FALSE)


set.seed(123456)
maxlod4=scanone(genotype_map, pheno.col=15001:20000, model=c("normal"), method=c("hk"), n.perm=1000)
tabbycol_15001<-summary(pheno_wetLOD4, format="tabByCol", pvalues=TRUE, perm=maxlod4, ci.function="lodint")
write.table(tabbycol_10001, file = "Tab_by_15001_Mckay_sum_SNP.txt", quote=FALSE, row.name=FALSE)

set.seed(123456)
maxlod5=scanone(genotype_map, pheno.col=20001:25663, model=c("normal"), method=c("hk"), n.perm=1000)
tabbycol_20001<-summary(pheno_wetLOD5, format="tabByCol", pvalues=TRUE, perm=maxlod5, ci.function="lodint")
write.table(tabbycol_10001, file = "Tab_by_20001_Mckay_sum_SNP.txt", quote=FALSE, row.name=FALSE)

save.image(pheno_wetLOD5, pheno_wetLOD4, pheno_wetLOD3, pheno_wetLOD2, pheno_wetLOD1, file='EQTL_wet.RData')

###Analyzing the different statistical methods in scanone function for locating QTLs and estimating their effects. 

#EM method maximum likelihood method, EM algorithm (Dempster et al. 1977)
out.em = scanone(genotype_map, pheno.col=3:10, method = 'em')
out.em[1:5,1:5]

#Haley-Knott Regression method (regression of the phenotypes on the multipoint QTL genotype probabilities): faster
out.hk <- scanone(genotype_map, method="hk", pheno.col = 3:10)
out.hk[1:5,1:5]

plot(out.em, out.hk, col=c("blue", "red")) #plotting the first phenotype
legend(8,0.7, legend=c("EM","HK"), fill=c("blue","red"), bty="n")

#multiple imputation is used, as described by Sen and Churchill (2001)
out.imp <- sim.geno(genotype_map, step=1, n.draws=64)
out.imp <- scanone(out.imp, method="imp", pheno.col = 3:10)

plot(out.em, out.hk, out.imp, col=c("green","blue", "red"))

legend(8,0.7, legend=c("EM","HK","imp"), fill=c("blue","red","green"), bty="n")

#############################################################################
# eQTL analysis for the DRY environment separated using r/qtl package
#############################################################################

#Downloading the linkage map and genotype of KAS-1 TSU-1
#Analizing first the Wet data
genotype_map_dry = read.cross("csvs", dir=".", genfile = "genotypes_match.csv", 
                          phefile = "pheno_dry.csv", header=FALSE, sep=',', dec=',' )


summary(genotype_map_dry)
nind(genotype_map_dry) #numer of individuals
nchr(genotype_map_dry) #number of chromossomes
totmar(genotype_map_dry) #total markers
nmar(genotype_map_dry) #number of markers per chromossome
nphe(genotype_map_dry) #didn't input the phenotype yet. Default is 1

#estimate recombination fraction?
#plotting the summary of the the data
#pdf("summary_map_dry")
plot(genotype_map_dry, pheno.col = 3) #Analizing just the first phenotype sample/column
#dev.off()

#plotMissing(genotype_map_dry)
pdf("Genetic Map")
plotMap(genotype_map_dry)
#dev.off()

genotype_map = calc.genoprob(genotype_map_dry, step=1, map.function="kosambi", error.prob = 0.01)

pheno_dryLOD1=scanone(genotype_map, pheno.col=2:5000, model=c("normal"), method=c("hk"))
pheno_dryLOD2=scanone(genotype_map, pheno.col=5001:10000, model=c("normal"), method=c("hk"))
pheno_dryLOD3=scanone(genotype_map, pheno.col=10001:15000, model=c("normal"), method=c("hk"))
pheno_dryLOD4=scanone(genotype_map, pheno.col=15001:20000, model=c("normal"), method=c("hk"))
pheno_dryLOD5=scanone(genotype_map, pheno.col=20001:25663, model=c("normal"), method=c("hk"))
transcriptLODall=cbind(pheno_dryLOD1, pheno_dryLOD2, pheno_dryLOD3, pheno_dryLOD4, pheno_dryLOD5)


# summarizing the data by extracting the lowest genome-wide permutation corrected P value per each expression trait. 
set.seed(123456)
maxlod1=scanone(genotype_map, pheno.col=2:5000, model=c("normal"), method=c("hk"), n.perm=1000)
tabbycol_5000 <-summary(pheno_dryLOD1, format="tabByCol", pvalues=TRUE, perm=maxlod1, ci.function="lodint")
write.table(tabbycol_5000, file = "Tab_by_5000_Mckay_sum_SNP_dry.txt", quote=FALSE, row.name=FALSE)


set.seed(123456)
maxlod2=scanone(genotype_map, pheno.col=5001:10000, model=c("normal"), method=c("hk"), n.perm=1000)
tabbycol_5001<-summary(pheno_wetLOD2, format="tabByCol", pvalues=TRUE, perm=maxlod2, ci.function="lodint")
write.table(tabbycol_5001, file = "Tab_by_5001_Mckay_sum_SNP_dry.txt", quote=FALSE, row.name=FALSE)


set.seed(123456)
maxlod3=scanone(genotype_map, pheno.col=10001:15000, model=c("normal"), method=c("hk"), n.perm=1000)
tabbycol_10001<-summary(pheno_wetLOD3, format="tabByCol", pvalues=TRUE, perm=maxlod3, ci.function="lodint")
write.table(tabbycol_10001, file = "Tab_by_10001_Mckay_sum_SNP_dry.txt", quote=FALSE, row.name=FALSE)


set.seed(123456)
maxlod4=scanone(genotype_map, pheno.col=15001:20000, model=c("normal"), method=c("hk"), n.perm=1000)
tabbycol_15001<-summary(pheno_wetLOD4, format="tabByCol", pvalues=TRUE, perm=maxlod4, ci.function="lodint")
write.table(tabbycol_10001, file = "Tab_by_15001_Mckay_sum_SNP_dry.txt", quote=FALSE, row.name=FALSE)

set.seed(123456)
maxlod5=scanone(genotype_map, pheno.col=20001:25663, model=c("normal"), method=c("hk"), n.perm=1000)
tabbycol_20001<-summary(pheno_wetLOD5, format="tabByCol", pvalues=TRUE, perm=maxlod5, ci.function="lodint")
write.table(tabbycol_10001, file = "Tab_by_20001_Mckay_sum_SNP_dry.txt", quote=FALSE, row.name=FALSE)


save.image(pheno_dryLOD5, pheno_dryLOD4, pheno_dryLOD3, pheno_dryLOD2, pheno_dryLOD1, file='EQTL_dry.RData')