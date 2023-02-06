#-----------------------------------------------------------------------------------------------------#
#							Call settings
#-----------------------------------------------------------------------------------------------------#
source("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Settings.R")


#-----------------------------------------------------------------------------------------------------#
#							Libraries
#-----------------------------------------------------------------------------------------------------#
library("ggplot2")
library("ggfortify")

#-----------------------------------------------------------------------------------------------------#
#							Input
#-----------------------------------------------------------------------------------------------------#
#temp_dat = Data # load(C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Data_QC/Salmon/Data.Rdata")
#temp_pheno = Pheno # load(C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Data_QC/Salmon/Pheno.Rdata")
# needs DE_res
load(file = paste0(s_OUT_dir,"Prep/logcpm.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/Pheno.Rdata"))


temp_pheno = Pheno 
#-----------------------------------------------------------------------------------------------------#
#							main
#-----------------------------------------------------------------------------------------------------#
# source: https://github.com/rcannood/SCORPIUS
expression = as.matrix(t(logcpm ))
plot(prcomp(t(logcpm ))$x)
plot(prcomp((expression ))$x)
group_name = temp_pheno$stage

# Only select the DMSO samples
expression_DMSO = expression[temp_pheno$Corisol==FALSE,]
group_name_DMSO = group_name[temp_pheno$Corisol==FALSE]
pheno_DMSO = temp_pheno[temp_pheno$Corisol==FALSE,]

# Reduce dimention here, check how many components contain a lot of info, then push into infer algorithm
# PCA
pca = prcomp(expression_DMSO,scale=TRUE)

# make autoplot of PCA to show variance explained on first 2 comp
autoplot(pca)

# See effect of components
plot(round(pca$sdev,3))
number_of_components = which(!abs(diff(round(pca$sdev,3)))>mean(abs(diff(round(pca$sdev,3)))))[1] - 1
# this checks whcih is the last point above average change, then gets this point -1
print(paste0("Estimated no comp = ",number_of_components))

# reduce dimentions
#space_DMSO <- reduce_dimensionality(expression_DMSO, "euclidean", ndim = number_of_components)
space_DMSO = pca$x[,1:number_of_components]

# infer trajectory
traj_DMSO <- infer_trajectory(space=space_DMSO,
							k=1, # set to 1!
							smoother = "lowess",
							approx_points = 100,
							thresh = 0.001,
							maxit=1000)
							
# Show trajectory
draw_trajectory_plot(space_DMSO, group_name_DMSO, traj_DMSO$path, contour = FALSE)
round(as.numeric(traj_DMSO$time),2)

# Save traj plot
Rplot(insert={draw_trajectory_plot(space_DMSO, group_name_DMSO, traj_DMSO$path, contour = FALSE)
},title="SCORPIUS_DMSO_pca_trajectory",resolution = 350, width = 480, height = 480)

# get Gene importance based on traj
gimp <- gene_importances(expression_DMSO, traj_DMSO$path, ntree = 10000, num_permutations = 1000, num_threads = 8)

# get significant genes
Gene_set = gimp$gene[gimp$pvalue<0.05]
print(paste0("Amount of genes sign = ",length(Gene_set)))

# calculate modules
Scaled_exp_DMSO = scale_quantile(expression_DMSO)
modules <- extract_modules(Scaled_exp_DMSO[,Gene_set], traj_DMSO$time, verbose = TRUE) # 
draw_trajectory_heatmap(expression_DMSO[,Gene_set], traj_DMSO$time, interaction(pheno_DMSO $Corisol,pheno_DMSO $stage), modules)
# Modules make sense; 7 reuslted. When checking dist between all genes, around 7 clusters identified!
# plot(hclust(dist(t(expression_DMSO[,Gene_set]))))

# 

Scorpius_obj = list(trajectory = traj_DMSO, modules = modules, Exp = expression_DMSO, Pheno = pheno_DMSO)
#-----------------------------------------------------------------------------------------------------#
#							Plotting
#-----------------------------------------------------------------------------------------------------#
# make heatmap
Rplot(insert={draw_trajectory_heatmap(expression_DMSO[,Gene_set], traj_DMSO$time, interaction(pheno_DMSO $Corisol,pheno_DMSO $stage), modules)
},title="SCORPIUS_DMSO_stage_heatmap",resolution = 350, width = 480, height = 480)


Rplot(insert={autoplot(pca)},title="SCORPIUS_DMSO_PCA_autoplot",resolution = 350, width = 480, height = 480)



modules$module[match(Gene_set,modules$feature)]

 

#-----------------------------------------------------------------------------------------------------#
#							output
#-----------------------------------------------------------------------------------------------------#
save.image(file = paste0(s_OUT_dir,"Scorpius/Scorpius_img.Rdata"))
save(Scorpius_obj,file = paste0(s_OUT_dir,"Scorpius/Scorpius_obj.Rdata"))

Relative_exp = scale_quantile(expression_DMSO[,Gene_set])
Relative_exp =data.frame(t(Relative_exp))
Relative_exp$module = modules$module[match(Gene_set,modules$feature)]
write.table(format(Relative_exp, digits=3),paste0(s_ROOT,s_out_folder,".FnT/Relative_exp.tsv"),col.names = TRUE,row.names = TRUE,sep = "\t",quote = FALSE)




#-----------------------------------------------------------------------------------------------------#
#							old
#-----------------------------------------------------------------------------------------------------#
## show modules
#refmod = data.frame(gene =modules$feature,mod = modules$module)
#plot(y = rep(NA,length(traj_DMSO$time)),x = as.numeric(traj_DMSO$time), type="l",col="white",ylim=range(Scaled_exp_DMSO))
#for( current_mod in 1:length(unique(refmod$mod))){
#	temp_dat = Scaled_exp_DMSO[,refmod[refmod$mod==current_mod,"gene"]]
#	apply(temp_dat,2,function(x){lines(y=x,as.numeric(traj_DMSO$time),type="l",col=alpha(current_mod,0.05))})
#	apply(temp_dat,2,function(x){lines(y=x,as.numeric(traj_DMSO$time),type="l",col=alpha(current_mod,0.05))})
#	lines(y=apply(temp_dat,1,mean) ,x=as.numeric(traj_DMSO$time),col=current_mod,lwd=2,lty=2)
#}