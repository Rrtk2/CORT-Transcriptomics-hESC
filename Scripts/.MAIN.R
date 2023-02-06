#-----------------------------------------------------------------------------------------------------#
#							Call settings
#-----------------------------------------------------------------------------------------------------#
source("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Settings.R")


#-----------------------------------------------------------------------------------------------------#
#							RUN MAIN ALGORITHM
#-----------------------------------------------------------------------------------------------------#
source(paste0(s_ROOTSCRIPTS,"Prep.R"))


source(paste0(s_ROOTSCRIPTS,"/DE.R")) # for extra


source(paste0(s_ROOTSCRIPTS,"/Scorpius.R"))


source(paste0(s_ROOTSCRIPTS,"/LME.R"))


source(paste0(s_ROOTSCRIPTS,"/GO.R"))


#source(paste0(s_ROOTSCRIPTS,"/GRN.R"))
