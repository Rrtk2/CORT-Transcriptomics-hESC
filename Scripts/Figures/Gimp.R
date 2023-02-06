
#-----------------------------------------------------------------------------------------------------#
# 							GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#
# File description:
#	Name
#		Gimp.R
#
#	Purpose 
# 		Dedicated script to generate Volcanoplot
#
# Author comment:
#	Rick A. Reijnders 
#	ra.reijnders@maastrichtuniversity.nl
#
#	Personal comment:
#
#
#-----------------------------------------------------------------------------------------------------#
#							Input
#-----------------------------------------------------------------------------------------------------#
load("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Data_QC/DE/DE_res_all.Rdata")
load(paste0(s_OUT_dir,"LME/LME_img.Rdata"))

#-----------------------------------------------------------------------------------------------------#
#							save
#-----------------------------------------------------------------------------------------------------#

write.table(format(data.frame(gimp), digits=3),paste0(s_ROOT,s_out_folder,".FnT/gimp.tsv"),col.names = TRUE,row.names = TRUE,sep = "\t",quote = FALSE)
