
#-----------------------------------------------------------------------------------------------------#
# 							GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#
# File description:
#	Name
#		Traj_modules.R
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
##-----------------------------------------------------------------------------------------------------#
#							Call settings
#-----------------------------------------------------------------------------------------------------#

source("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Settings.R")
load(file = paste0(s_OUT_dir,"LME/lme_table.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/logcpm.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/Pheno.Rdata"))
load(file = paste0(s_OUT_dir,"Scorpius/Scorpius_obj.Rdata"))
load(file = paste0(s_OUT_dir,"LME/LME_img.Rdata"))


#-----------------------------------------------------------------------------------------------------#
#							Libraries
#-----------------------------------------------------------------------------------------------------#
library(tidyverse)
library(ggplot2)

#-----------------------------------------------------------------------------------------------------#
#							MAIN
#-----------------------------------------------------------------------------------------------------#
Gene_set = Scorpius_obj$modules$feature

scaled_data = t(scale_quantile(expression_DMSO[,Gene_set]))

Scorpius_obj$modules$module

modres = list()
for( i in 1:max(Scorpius_obj$modules$module)){
	temp_select = Scorpius_obj$modules$module == i

	modres[[paste0("Module ",i)]] = apply(scaled_data[temp_select,],2,mean)
}

modres = as.data.frame(modres)

modres_avg = rbind(pro = (apply(modres[1:3,],2,mean)),diffy = (apply(modres[4:6,],2,mean)),diffm = (apply(modres[7:9,],2,mean)))
Pseudotime_avg = c(mean(Scorpius_obj$trajectory$time[1:3]),mean(Scorpius_obj$trajectory$time[4:6]),mean(Scorpius_obj$trajectory$time[7:9]))

plot(x=Pseudotime_avg,modres_avg[,1],col=1,lty=1,type="l")
for(i in 2:dim(modres_avg)[2]){
	lines(Pseudotime_avg,y=modres_avg[,i],col=i)
}

#-----------------------------------------------------------------------------------------------------#
#							GGplot variant AVG
#-----------------------------------------------------------------------------------------------------#

library(ggplot2)
library(ggformula)
library(firatheme)
s_margin = 1

# make data in long
modres_avg2 = data.frame(modres_avg)
modres_avg2$sample = rownames(modres_avg2)



temp = pivot_longer(data = modres_avg2,cols = c(!sample),names_to = "Module",values_to = "Relative expression")
temp = data.frame(temp)
temp$pseudo =  case_when(
  temp$sample == "pro" ~ Pseudotime_avg[1],
  temp$sample == "diffy" ~ Pseudotime_avg[2],
  temp$sample == "diffm" ~ Pseudotime_avg[3])


# Basic line plot with points
g1 = ggplot(data=temp, aes(x=pseudo, y=Relative.expression, color=Module)) +
  geom_line(lwd=1.)+
  geom_point()+
  labs(x = "Pseudotime",y="Average relative expression",title= "Average relative expression of modules driving neuronal differentiation")+
  scale_x_continuous(limits = c(0, 1),expand = c(-s_margin, s_margin)) + scale_y_continuous(limits = c(0, 1),expand = c(-s_margin, s_margin)); g1
  
  
Rplot(g1,"Traj_modules",width = 480*2)

#-----------------------------------------------------------------------------------------------------#
#							GGplot variant ALL
#-----------------------------------------------------------------------------------------------------#

library(ggplot2)
library(ggformula)
library(firatheme)
s_margin = 1

# make data in long
modres3 = data.frame(modres)
modres3$sample = rownames(modres3)
modres3$pseudotime = Scorpius_obj$trajectory$time
modres3$sd = apply(modres3[,1:8],2,sd)

temp2 = pivot_longer(data = modres3,cols = c(1:8),names_to = "Module",values_to = "Relative expression")
temp2 = data.frame(temp2)



#temp$pseudo =  case_when(
#  temp$sample == "pro" ~ Pseudotime_avg[1],
#  temp$sample == "diffy" ~ Pseudotime_avg[2],
 # temp$sample == "diffm" ~ Pseudotime_avg[3])


# Basic line plot with points
g2 = ggplot() +
  geom_line(data=temp, aes(x=pseudo, y=Relative.expression, color=Module),lwd=1.5, alpha = 0.5)+
  geom_point(data=temp2, aes(x=pseudotime, y=Relative.expression, color=Module,group  = sample), size = 3, alpha = 0.2)+
  theme_minimal()+
  labs(x = "Pseudotime",y="Relative expression",title= "Relative expression of modules driving neuronal differentiation"); g2
  

Rplot2(g2,"Traj_modules",width = 480*2,s_figure_folder = "C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Data_QC/test/.FnT/")
