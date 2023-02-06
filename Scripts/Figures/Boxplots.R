
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

# Define things to be compared
temp_stat_comp_list = list(c(1,2),c(3,4),c(5,6))
	
# Change box plot colors by groups
selectedgene = pdy_genes$feature[1]

Plotgene = function(selectedgene = pdy_genes$feature[1], genename = NA){
 plot_1 = ggplot(pdy_data_long[pdy_data_long$Gene==selectedgene,] , aes(x=Stage_cond, y=Expression, fill=Condition)) +
  geom_boxplot()+
  ggtitle(paste0(selectedgene," (",genename,")"))+
  stat_compare_means(comparisons = temp_stat_comp_list,method = "t.test",label = "p.signif")+
  xlab("Stage and Condition")+
  ylab("Expression (log2 cpm)")+
  theme( axis.text.x=element_blank(),
  axis.ticks.x=element_blank())+
  #scale_x_discrete(breaks=seq(1,12,1),labels= c("a","b"))
  geom_text(x=1.5, y=round(range(pdy_data_long[pdy_data_long$Gene==selectedgene,]$Expression)[1],0)+0.05, label="Pro",family = "sans",size=5)+
  geom_text(x=3.5, y=round(range(pdy_data_long[pdy_data_long$Gene==selectedgene,]$Expression)[1],0)+0.05, label="Diffy",family = "sans",size=5)+
  geom_text(x=5.5, y=round(range(pdy_data_long[pdy_data_long$Gene==selectedgene,]$Expression)[1],0)+0.05, label="Diffm",family = "sans",size=5)
  
  return(plot_1)
}


Plotgene(selectedgene = pdy_genes$feature[1], genename = "LRRTM2")
			  
# Change box plot colors by groups
ggplot(pdy_data_long[pdy_data_long$Stage=="diffy",] , aes(x=Gene, y=Expression, fill=Condition)) +
  geom_boxplot()

# Change box plot colors by groups
ggplot(pdy_data_long[pdy_data_long$Stage=="diffm",] , aes(x=Gene, y=Expression, fill=Condition)) +
  geom_boxplot()



lme_table$feature
