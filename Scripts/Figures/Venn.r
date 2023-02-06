#-----------------------------------------------------------------------------------------------------#
#							Call settings
#-----------------------------------------------------------------------------------------------------#
source("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Settings.R")

#-----------------------------------------------------------------------------------------------------#
#							Input
#-----------------------------------------------------------------------------------------------------#
load(file = paste0(s_OUT_dir,"LME/lme_table.Rdata"))

x = list( pro_diffy = lme_table$feature[lme_table$pro_diffy.p.value<0.05],
diffy_diffm = lme_table$feature[lme_table$diffy_diffm.p.value<0.05])


  
#devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF","#CD534CFF"),
  stroke_size = 0.8, set_name_size = 8
  )