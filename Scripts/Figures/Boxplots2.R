
#-----------------------------------------------------------------------------------------------------#
# 							GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#
# File description:
#	Name
#		Volcanoplot.R
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
#source("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Settings.R")
source("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Settings.R")
load(file = paste0(s_OUT_dir,"LME/lme_table.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/logcpm.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/Pheno.Rdata"))
load(file = paste0(s_OUT_dir,"Scorpius/Scorpius_obj.Rdata"))
load(file = paste0(s_OUT_dir,"LME/LME_img.Rdata"))
library(ggpubr)   
library(ggplot2)
library(gridExtra)
library(cowplot)

#-----------------------------------------------------------------------------------------------------#
#							Functions
#-----------------------------------------------------------------------------------------------------#

Plotgene = function(datalong = pdy_data_long, selectedgene = pdy_genes$feature[1], genename = NA){
	plot_1 = ggplot(datalong[datalong$Gene==selectedgene,] ,aes(x=Condition, y=Expression, fill = Condition)) +
	geom_boxplot()+ #fill = Condition
	stat_compare_means(comparisons  = temp_stat_comp_list, method = "t.test", label = "p.signif")+
	ggtitle(paste0(selectedgene," (",genename,")"))+
	xlab("Stage")+
	ylab("Expression (log2 cpm)")+
	theme( axis.text.x=element_blank(),
	axis.ticks.x=element_blank())+
	#scale_x_discrete(breaks=seq(1,12,1),labels= c("a","b"))
	#facet_grid(cols =vars(Stage), scales = "free")+
	facet_wrap(~Stage, scales = "free_x")
   
  return(plot_1)
}


#-----------------------------------------------------------------------------------------------------#
#							pro-diffy
#-----------------------------------------------------------------------------------------------------#

# get genes for pro-diffy
pdy_genes =  lme_table[lme_table$pro_diffy.p.value<0.05,]
pdy_genes = pdy_genes[order(pdy_genes $pro_diffy.p.value),]

pdy_genes = pdy_genes[1:2,]

pdy_data = logcpm[match(pdy_genes$feature,rownames(logcpm)),]
pdy_data = t(pdy_data)
pdy_data = data.frame(pdy_data, data.frame(t(as.data.frame(strsplit(rownames(pdy_data),split = "_")))))
pdy_data$stage = substr(pdy_data$X2,1,nchar(pdy_data$X2)-1)
pdy_data_long = pivot_longer(pdy_data,c(1,2))
colnames(pdy_data_long) = c("Condition", "ID", "Stage", "Gene", "Expression")
pdy_data_long$Condition = factor(pdy_data_long$Condition,levels=c("D","C"))
levels(pdy_data_long$Condition) = c("DMSO","Cortisol")

pdy_data_long$Stage_cond = factor(interaction(pdy_data_long$Stage,pdy_data_long$Condition),levels=c("pro.DMSO","pro.Cortisol","diffy.DMSO","diffy.Cortisol","diffm.DMSO","diffm.Cortisol"))
pdy_data_long$Stage = factor(pdy_data_long$Stage,levels=c("pro","diffy","diffm"))

# Define things to be compared
temp_stat_comp_list = list(c(1,2))


Plotgene(datalong = pdy_data_long, selectedgene = pdy_genes$feature[1], genename = "LRRTM2")
Plotgene(datalong = pdy_data_long, selectedgene = pdy_genes$feature[2], genename = "TSPAN5")


Rplot2({Plotgene(datalong = pdy_data_long, selectedgene = pdy_genes$feature[1], genename = "LRRTM2")},	title = "box_LRRTM2_pdy")
Rplot2({Plotgene(datalong = pdy_data_long, selectedgene = pdy_genes$feature[2], genename = "TSPAN5")},	title = "box_TSPAN5_pdy")

#-----------------------------------------------------------------------------------------------------#
#							diffy-diffm
#-----------------------------------------------------------------------------------------------------#

# get genes for pro-diffy
dydm_genes =  lme_table[lme_table$diffy_diffm.p.value<0.05,]
dydm_genes = dydm_genes[order(dydm_genes $diffy_diffm.p.value),]

dydm_genes = dydm_genes[1:3,]

dydm_data = logcpm[match(dydm_genes$feature,rownames(logcpm)),]
dydm_data = t(dydm_data)
dydm_data = data.frame(dydm_data, data.frame(t(as.data.frame(strsplit(rownames(dydm_data),split = "_")))))
dydm_data$stage = substr(dydm_data$X2,1,nchar(dydm_data$X2)-1)
dydm_data_long = pivot_longer(dydm_data,c(1,2,3))
colnames(dydm_data_long) = c("Condition", "ID", "Stage", "Gene", "Expression")
dydm_data_long$Condition = factor(dydm_data_long$Condition,levels=c("D","C"))
levels(dydm_data_long$Condition) = c("DMSO","Cortisol")

dydm_data_long$Stage_cond = factor(interaction(dydm_data_long$Stage,dydm_data_long$Condition),levels=c("pro.DMSO","pro.Cortisol","diffy.DMSO","diffy.Cortisol","diffm.DMSO","diffm.Cortisol"))
dydm_data_long$Stage = factor(dydm_data_long$Stage,levels=c("pro","diffy","diffm"))

# Define things to be compared
temp_stat_comp_list = list(c(1,2))
	
Plotgene(datalong = dydm_data_long, selectedgene = dydm_genes$feature[1], genename = "KCND3")
Plotgene(datalong = dydm_data_long, selectedgene = dydm_genes$feature[2], genename = "KCNIP4")
Plotgene(datalong = dydm_data_long, selectedgene = dydm_genes$feature[3], genename = "GRIA3")
		  

Rplot2({Plotgene(datalong = dydm_data_long, selectedgene = dydm_genes$feature[1], genename = "KCND3")},	title = "box_KCND3_dydm")
Rplot2({Plotgene(datalong = dydm_data_long, selectedgene = dydm_genes$feature[2], genename = "KCNIP4")},	title = "box_KCNIP4_dydm")
Rplot2({Plotgene(datalong = dydm_data_long, selectedgene = dydm_genes$feature[3], genename = "GRIA3")},	title = "box_GRIA3_dydm")