#-----------------------------------------------------------------------------------------------------#
#							SEED
#-----------------------------------------------------------------------------------------------------#
# set seed for consistency
s_seed <<- 42

set.seed(s_seed)
#-----------------------------------------------------------------------------------------------------#
#							ROOT FOLDER
#-----------------------------------------------------------------------------------------------------#
# define ROOT here 
s_ROOT = "C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/"
s_ROOTSCRIPTS = paste0(s_ROOT,"Scripts/")

#-----------------------------------------------------------------------------------------------------#
#							DEFINE OUTPUT FOLDER
#-----------------------------------------------------------------------------------------------------#
# change this to generate a different run, stored in new location
s_out_folder <<- "Data_QC/test/" 

# define out loc
s_OUT_dir <<- paste0(s_ROOT,s_out_folder)

#-----------------------------------------------------------------------------------------------------#
#							Require libraries
#-----------------------------------------------------------------------------------------------------#
source(paste0(s_ROOTSCRIPTS,"Libraries.R"))
source(paste0(s_ROOTSCRIPTS,"Functions.R"))

#-----------------------------------------------------------------------------------------------------#
#							MAKE OUTPUT FOLDERS IF NEEDED
#-----------------------------------------------------------------------------------------------------#
temp_all_dirs = c(".FnT","Prep","DE","Scorpius","LME","GO","PPI","eQTL","Coloc")
if(!dir.exists(paste0(s_OUT_dir))){dir.create(file.path(paste0(s_OUT_dir)))}
lapply(temp_all_dirs,function(i){if(!dir.exists(paste0(s_OUT_dir,i))){dir.create(file.path(paste0(s_OUT_dir,i)))}})

#-----------------------------------------------------------------------------------------------------#
#							GENERAL SETTINGS
#-----------------------------------------------------------------------------------------------------#

s_fcTH  = 1.5
s_pvalTH = 0.05
#-----------------------------------------------------------------------------------------------------#
#							Make SETTINGS file (IF NOT EXISTS)
#-----------------------------------------------------------------------------------------------------#

if(!file.exists(paste0(s_OUT_dir,"SETTINGS.txt"))){
	temp_settings_names = ls()[grep(ls(),pattern = "^s_")]

	temp_set = data.frame(variable = temp_settings_names, value = NA)

	for(i in 1:nrow(temp_set)){
		temp_set[i,"value"] =  get(temp_set[i,"variable"] )
	}

	write.table("SETTINGS:",file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
	write.table(temp_set,file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE,append = TRUE)
	write.table("\nVERSION:",file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
	write.table(data.frame(name=names(version),value=unlist(version)),file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
	write.table("\nSYSTEM INFO:",file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
	write.table(data.frame(name=names(Sys.info()),value=unlist(Sys.info())),file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
	write.table("\nBASE PACKAGES:",file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
	write.table(sessionInfo()$basePkgs,file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
	write.table("\nOTHER PACKAGES:",file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)
	write.table(names(sessionInfo()$otherPkgs),file = paste0(s_OUT_dir,"SETTINGS.txt"),sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE,append = TRUE)

}