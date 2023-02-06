
#-----------------------------------------------------------------------------------------------------#
# 							GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#
# File description:
#	Name
#		LME_table.R
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

if(FALSE){
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
#							Heatmap
#-----------------------------------------------------------------------------------------------------#
geneset = lme_table$feature[lme_table$pro_diffy.p.value<0.05 & abs(lme_table$pro_diffy.Value)>0.5 | lme_table$diffy_diffm.p.value<0.05 & abs(lme_table$diffy_diffm.Value)>0.5]
genedata = logcpm[geneset,]

# DMSO scaled set
#a=scale(t(genedata[,Pheno$Corisol==F]))
a=scale(t(genedata[,]))
gene_scale = attributes(a)$`scaled:scale`
gene_center = attributes(a)$`scaled:center`

genedata2 = t(scale(t(genedata),scale = gene_scale, center = gene_center))

colnames(genedata2) = colnames(genedata)
genedata = data.frame(genedata2 )
genedata$Gene = rownames(genedata )


genedata_heatmap = genedata %>%
	pivot_longer(!Gene, names_to = "Sample", values_to = "Expression") %>%
	mutate(COND = factor(ifelse(grepl(Sample,pattern = "C_"),yes = "CORT",no = "DMSO"),levels=c("DMSO","CORT")))  %>%
	mutate(TIME = factor(gsub(gsub(Sample,pattern = "._",replacement = ""),pattern = ".$",replacement = ""),levels=c("pro","diffy","diffm")))  %>%
	mutate(REP = gsub(Sample,pattern = "[A-z]",replacement = ""))  %>%
	mutate(Sample = as.character(Sample))%>%
	mutate(MODULE = factor(data.frame(Scorpius_obj$modules[match(Gene,Scorpius_obj$modules$feature),"module"])[,1],levels=1:max(Scorpius_obj$modules$module)))#%>%
	#mutate(Gene = fct_reorder(Gene, order(MODULE)))%>%
	#mutate(Sample = fct_reorder(Sample, order(TIME,REP))) 

# reorder samples
genedata_heatmap =  genedata_heatmap[ order(genedata_heatmap$COND,genedata_heatmap$TIME,genedata_heatmap$REP ),]

ggplot(genedata_heatmap, aes(x = Sample, y = Gene, fill = Expression, group = MODULE)) + 
	facet_grid( MODULE ~ COND+TIME, scales = "free") + #, switch = "x", scales = "free_x", space = "free_x"
	#facet_grid(~ MODULE) +
	scale_fill_gradient2(low="navy", mid="white", high="red",midpoint = 0,aesthetics = "fill", na.value = NA)+ #limit = c(-round(max(abs(range(a)))*1.5,3),round(max(abs(range(a)))*1.5,3))
	geom_tile()


#-----------------------------------------------------------------------------------------------------#
#							Process table to LONG, so ggplot can understand
#-----------------------------------------------------------------------------------------------------#
# we want X=feature; Ymean = Value; Xsd = std error; Moment = moment 1 or 2

temp = pivot_longer(data = lme_table,cols = c("pro_diffy.Value","diffy_diffm.Value"),names_to = "Moment",values_to = "FC")
temp2 = pivot_longer(data = temp,cols = c("pro_diffy.Std.Error","diffy_diffm.Std.Error"),names_to = "stde",values_to = "stde2")

ggplot(data=temp[1:10,], aes(x=feature,y=FC,col=Moment))+geom_point()


ggplot(data=temp2[1:20,], aes(x=feature,y=FC,col=Moment, ymin = FC-stde2, ymax = FC+stde2))+
	geom_hline(yintercept = 0,lty=2)+
	geom_errorbar(width=0.2,position = position_dodge(0.3))+
	geom_point(size=2,position = position_dodge(0.3))+
	coord_flip()




#-----------------------------------------------------------------------------------------------------#
#							Main
#-----------------------------------------------------------------------------------------------------#
lme_table

}
cat("this script is outdated\n")