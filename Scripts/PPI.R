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
#							Just get the relevant parts and save into file to be put in ppi manually
#-----------------------------------------------------------------------------------------------------#

py = lme_table$feature[which(lme_table$pro_diffy.p.value<0.05)]

ym = lme_table$feature[which(lme_table$diffy_diffm.p.value<0.05)]

#-----------------------------------------------------------------------------------------------------#
#							output
#-----------------------------------------------------------------------------------------------------#
write.table(format(py, digits=3),paste0(s_ROOT,s_out_folder,"PPI/Pro_diffy.tsv"),col.names = FALSE,row.names = FALSE,sep = "\t",quote = FALSE)
write.table(format(ym, digits=3),paste0(s_ROOT,s_out_folder,"PPI/Diffy_diffm.tsv"),col.names = FALSE,row.names = FALSE,sep = "\t",quote = FALSE)


#-----------------------------------------------------------------------------------------------------#
#							out
#-----------------------------------------------------------------------------------------------------#

# get upper 10% -> 1500 genes?
#linkList <- getLinkList(weightMat,threshold = quantile(weightMat,0.9))
#dim(linkList)