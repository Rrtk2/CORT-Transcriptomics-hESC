#-----------------------------------------------------------------------------------------------------#
#							Call settings
#-----------------------------------------------------------------------------------------------------#
source("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Settings.R")



#-----------------------------------------------------------------------------------------------------#
#							Input
#-----------------------------------------------------------------------------------------------------#
load(file = paste0(s_OUT_dir,"Scorpius/Scorpius_obj.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/logcpm.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/Pheno.Rdata"))
load(file = paste0(s_OUT_dir,"LME/lme_table.Rdata"))

#-----------------------------------------------------------------------------------------------------#
#							Get relevant genes for pro-diffy
#-----------------------------------------------------------------------------------------------------#
genes_pro_diffy =  lme_table[which(lme_table$pro_diffy.p.value<0.05),"feature"]
genes_diffy_diffm =  lme_table[which(lme_table$diffy_diffm.p.value<0.05),"feature"]
allgenes = unique(c(genes_pro_diffy,genes_diffy_diffm))

#paste0(allgenes,collapse = "@")
# put it in excel?
# interface is going to be faster like this...

#-----------------------------------------------------------------------------------------------------#
#							Get eQTL data
#-----------------------------------------------------------------------------------------------------#
#https://www.metabrain.nl/
#https://download.metabrain.nl/files.html

# VIA ONLINE BROWSER....

#-----------------------------------------------------------------------------------------------------#
#							Get all grenerated csv files and make 1 big file data
#-----------------------------------------------------------------------------------------------------#
# basically the genes with any cis eQTL is included in this set of files, it does not fully overlap for all genes, some genes didnt have any. eQTL per gene varied between 1 and 800

allcsvfiles = list.files("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Anno/eQTL metabrain/eqtl/")

eQTLdata = list()
for(i in allcsvfiles){
	eQTLdata[[i]] = read.csv(paste0("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Anno/eQTL metabrain/eqtl/",i))
}

eQTLdata = do.call(rbind.data.frame, eQTLdata)
rownames(eQTLdata) = eQTLdata$ID

#-----------------------------------------------------------------------------------------------------#
#							output
#-----------------------------------------------------------------------------------------------------#
save(eQTLdata,file = paste0(s_OUT_dir,"eQTL/eQTLdata.Rdata"))

#-----------------------------------------------------------------------------------------------------#
#							out
#-----------------------------------------------------------------------------------------------------#

# get upper 10% -> 1500 genes?
#linkList <- getLinkList(weightMat,threshold = quantile(weightMat,0.9))
#dim(linkList)