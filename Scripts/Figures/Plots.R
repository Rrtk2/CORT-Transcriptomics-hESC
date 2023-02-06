#-----------------------------------------------------------------------------------------------------#
#							Call settings
#-----------------------------------------------------------------------------------------------------#
source("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Settings.R")

#-----------------------------------------------------------------------------------------------------#
#							Input
#-----------------------------------------------------------------------------------------------------#
load(file = paste0(s_OUT_dir,"Prep/logcpm.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/Pheno.Rdata"))
load(file = paste0(s_OUT_dir,"LME/Genes_of_interest.Rdata"))
load(file = paste0(s_OUT_dir,"Scorpius/Scorpius_img.Rdata"))

#-----------------------------------------------------------------------------------------------------#
#							paired boxplot
#-----------------------------------------------------------------------------------------------------#
library("ggplot2")

plotggboxplot = function(gene){
	ggdata = cbind(Pheno,exp = logcpm[gene,])

	e <- ggplot(ggdata , aes(x = stage, y = exp, fill = Corisol))+
			geom_boxplot() + ggtitle(gene)
	e
}

mod = modules[match(Genes_of_interest[["intersect"]],modules$feature),]
mod2 = mod[mod$module==6,]


plotggboxplot(mod2$feature[1])
plotggboxplot(mod2$feature[2])
plotggboxplot(mod2$feature[3])

#-----------------------------------------------------------------------------------------------------#
#							output
#-----------------------------------------------------------------------------------------------------#
save.image(file = paste0(s_OUT_dir,"LME/LME_img.Rdata"))
save(lme_table,file = paste0(s_OUT_dir,"LME/lme_table.Rdata"))


