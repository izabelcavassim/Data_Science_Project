######################################################################
# Project: Data Science in Bioinformatics - Visualization and Analysis 
# Student: Maria Izabel Cavassim Alves
# Professor: Thomas Mailund
# Identifying the eQTLs from each module by taking their positions
######################################################################

# Analazing the qtl positions for each phenotype;
require('qtl')
setwd("C:/Users/Maria/Desktop/MASTER/DATA_SCIENCE/Data/eQTL mapping/")
load("EQTL_dry.RData")

##########################################
#Function for identification of the transcripts based on the maxlod
##########################################

unlist_select = function(pheno, max){
mysummary = summary(pheno, format='allpeaks', pvalues = TRUE, perm =maxlod1)
unlist = unlist(mysummary)
names(unlist[which(!is.na(sapply(names(unlist), function(x) pmatch("AT",x))))]) -> transcripts
as.numeric(unlist[which(!is.na(sapply(names(unlist), function(x) pmatch("pos",x))))]) -> positions
as.numeric(unlist[which(!is.na(sapply(names(unlist), function(x) pmatch("pval",x))))]) -> pvals

#taking the positions of p-values which are lower than 0.05
threshold_positions = which(as.numeric(unlist[which(!is.na(sapply(names(unlist), function(x) pmatch("pval",x))))]) < 0.05)
pvals_match = pvals[threshold_positions]
positions_match = positions[threshold_positions]
transcripts_match = transcripts[threshold_positions]

#taking the last character of the string indentified as chromossome
chromossomes_match = as.numeric(substring(transcripts_match[1:length(transcripts_match)], 10))
eqtl = mapply(rbind, positions_match, transcripts_match, chromossomes_match, pvals_match, SIMPLIFY=T)
return(eqtl)
}  
  
#it is necessary to load the qtl package before  
eqtl1 = unlist_select(pheno = pheno_dryLOD1, max = maxlod1)  
eqtl2 = unlist_select(pheno = pheno_dryLOD2, max = maxlod2)
eqtl3 = unlist_select(pheno = pheno_dryLOD3, max = maxlod3)
eqtl4 = unlist_select(pheno = pheno_dryLOD4, max = maxlod4)
eqtl5 = unlist_select(pheno = pheno_dryLOD5, max = maxlod5)
  
####################################################  
#putting all the identified eqtls in a unique vector
pvalues_total = as.numeric(c(eqtl_1[4,], eqtl_2[4,], eqtl_3[4,], eqtl_4[4,], eqtl_5[4,])) 
position_total = as.numeric(c(eqtl_1[1,], eqtl_2[1,], eqtl_3[1,], eqtl_4[1,], eqtl_5[1,])) 
transcripts_names_total = c(eqtl_1[2,], eqtl_2[2,], eqtl_3[2,], eqtl_4[2,], eqtl_5[2,])
chromossomes_total = c(eqtl_1[3,], eqtl_2[3,], eqtl_3[3,], eqtl_4[3,], eqtl_5[3,])

total_eqtl = mapply(rbind,transcripts_names_total, pvalues_total, position_total, chromossomes_total)
total_eqtl = t(total_eqtl)
save(total_eqtl, file='Total_eqtl_dry.RData')
load('Total_eqtl_dry.RData')

#the p-values are localized at 4th row from etql_[1:5] object
#install.packages('fdrtool')
require(fdrtool)

pvalues = as.numeric(total_eqtl[,2])
positions = as.numeric(total_eqtl[,3])
write.table(pvalues, file='eqtl_todos_pvalues_dry.txt')

#pval.estimate.eta0(pvalues, method="conservative")
#pval.estimate.eta0(pvalues, method="adaptive")
#pval.estimate.eta0(pvalues, method="bootstrap")
fdrtool(pvalues, statistic="pvalue")
fndr.cutoff(pvalues, statistic= "pvalue") # cutoff is 0.049
p_adjust = p.adjust(pvalues, method = 'BH', n = length(pvalues))

length(total_eqtl[,2])
position = which(p_adjust <= 0.01)
#total_eqtl_names = total_eqtl[position,1] # just the names of the trancripts which has p-value lower than 0.01
total_eqtl_pvalues = total_eqtl[position,1:4] #defining the threshold of 0.01 and extracting the position of the eQTL
length(total_eqtl_pvalues[,1])

z = load('C:/Users/Maria/Desktop/MASTER/DATA_SCIENCE/Data/Coexpression/Modules_transcriptname_dry.RData') # the data with the modules classified by their color in wet data
list(z)

# taking off the last element of the transcript name
x = rownames(total_eqtl_pvalues) # the rownames continue with the chromossome number
total_names = substr(x, 0,9) 

length(module_blue)
module_bluematch = total_eqtl_pvalues[match(module_blue,total_names),] # matching module blue names with the eqtl names significant
module_bluematch = na.omit(module_bluematch)
length(module_bluematch[,1])
View(module_bluematch)

#Defining if the eQTL is cis or trans
#Taking the third word of the transcript name and looking at the lineage group
# If lineage group is the same as the chromossome == 'cis'
# Else == 'trans'
############################################################
cis_trans = function(module){
x = module[,1]
chr = substr(x, 3,3) 
lineage_group = module[,4]
length(chr)
list= vector()
for(i in 1:length(chr)){
if(chr[i] == lineage_group[i]){
  list[i] = 'cis';
}else{
  list[i] = 'trans'
  }
}
#defining the col names
names = c('Transcript', 'P-value', 'eQTL position', 'Lineage group', 'Cis/Trans')
module = cbind(module,list)
colnames(module) = names
return(module)
}

module_bluematch = cis_trans(module_bluematch)
module_bluematch
############################################################
#Creating a function to detect the EQTL significant transcripts in each module

match_mod = function(modules, total){
  total_names = list()
  # taking off the last element of the transcript name
  x = rownames(total) # the rownames continue with the chromossome number
  total_names = substr(x, 0,9) 
  module_match = total[match(modules,total_names),]
  module_match = na.omit(module_match)
  module_match = cis_trans(module_match)
  if(length(module_match) > 5){
      module_match = module_match[order(module_match[,3]),];
    }
  return(module_match)
}
#############################################################  
module_redmatch = match_mod(modules = module_red, total = total_eqtl_pvalues) 
module_bluematch = match_mod(modules = module_blue, total = total_eqtl_pvalues) 
module_greenmatch = match_mod(modules = module_green, total = total_eqtl_pvalues)
module_yellowmatch = match_mod(modules = module_yellow, total = total_eqtl_pvalues)
module_blackmatch = match_mod(modules = module_black, total = total_eqtl_pvalues)
module_pinkmatch = match_mod(modules = module_pink, total = total_eqtl_pvalues)
module_purplematch = match_mod(modules = module_purple, total = total_eqtl_pvalues)
module_magentamatch = match_mod(modules = module_magenta, total = total_eqtl_pvalues)
module_turquoisematch = match_mod(modules = module_turquoise, total = total_eqtl_pvalues)
module_salmonmatch = match_mod(modules = module_salmon, total = total_eqtl_pvalues)
module_greymatch = match_mod(modules = module_grey, total = total_eqtl_pvalues)
module_brownmatch = match_mod(modules = module_brown, total = total_eqtl_pvalues)
module_brazilmatch = match_mod(modules = module_brazil, total = total_eqtl_pvalues)
module_tanmatch = match_mod(modules= module_tan, total = total_eqtl_pvalues)

save(module_bluematch, module_brownmatch, module_greenmatch, module_yellowmatch, module_blackmatch, module_tanmatch, module_pinkmatch, module_salmonmatch, module_turquoisematch, module_redmatch, module_brazilmatch, module_purplematch, module_greymatch, module_magentamatch, file="Modules_detected_dry.RData")
write.csv(module_bluematch, file='modulos_match_qtl_blue_dry.csv')
write.csv(module_yellowmatch, file='modulos_match_qtl_yellow_dry.csv')
write.csv(module_pinkmatch, file='modulos_match_qtl_pink_dry.csv')
write.csv(module_turquoisematch, file='modulos_match_qtl_turquoise_dry.csv')