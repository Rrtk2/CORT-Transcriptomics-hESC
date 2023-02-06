#-----------------------------------------------------------------------------------------------------#
#							Call settings
#-----------------------------------------------------------------------------------------------------#
source("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Settings.R")

#-----------------------------------------------------------------------------------------------------#
#							Input
#-----------------------------------------------------------------------------------------------------#
load(file = paste0(s_OUT_dir,"Prep/logcpm.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/Pheno.Rdata"))
load(file = paste0(s_OUT_dir,"Scorpius/Scorpius_obj.Rdata"))





#-----------------------------------------------------------------------------------------------------#
#							main ALL + auto RELEVEL
#-----------------------------------------------------------------------------------------------------#
library(nlme)


lme_table = data.frame(feature=Scorpius_obj$modules$feature, pro_diffy.Value = NA, pro_diffy.Std.Error = NA, pro_diffy.DF = NA, pro_diffy.t.value = NA, pro_diffy.p.value = NA, diffy_diffm.Value = NA, diffy_diffm.Std.Error = NA, diffy_diffm.DF = NA, diffy_diffm.t.value = NA, diffy_diffm.p.value = NA)
	for(i in 1:dim(lme_table)[1]){
		o = logcpm[Scorpius_obj$modules$feature,][i,]
		model_pro_diffy <- lme(y ~ Cort*Stage, random=~ 1|ID, data = data.frame(y = o,
                                                     ID = gsub(pattern = "[A-z]",replacement = "",x = rownames(Pheno)),Cort = Pheno$Corisol, 
													 Stage = factor(Pheno$stage,levels = c("pro", "diffy", "diffm"))),control = lmeControl(opt = 'optim')) #to run the model
		model_diffy_diffm <- lme(y ~ Cort*Stage, random=~ 1|ID, data = data.frame(y = o,
                                                     ID = gsub(pattern = "[A-z]",replacement = "",x = rownames(Pheno)),Cort = Pheno$Corisol, 
													 Stage = factor(Pheno$stage,levels = c("diffy", "diffm", "pro"))),control = lmeControl(opt = 'optim')) #to run the model
													 
		temp_result = data.frame(pro_diffy = t(summary(model_pro_diffy)$tTable["CortTRUE:Stagediffy",]), diffy_diffm = t(summary(model_diffy_diffm)$tTable["CortTRUE:Stagediffm",]))
												 
													 
		lme_table[i,(2:dim(lme_table)[2])] = temp_result#anova(model)$'p-value'[4]
	}

#-----------------------------------------------------------------------------------------------------#
#							Print outputs
#-----------------------------------------------------------------------------------------------------#
Genes_of_interest = list()

# pro_diffy
table(lme_table$pro_diffy.p.value<0.05)
Genes_of_interest[["pro_diffy"]] = lme_table$feature [which(lme_table$pro_diffy.p.value<0.05)]

# diffy_diffm
table(lme_table$diffy_diffm.p.value<0.05)
Genes_of_interest[["diffy_diffm"]] = lme_table$feature [which(lme_table$diffy_diffm.p.value<0.05)]

# those who inersect
Genes_of_interest[["intersect"]]  =intersect(lme_table$feature[which(lme_table$pro_diffy.p.value<0.05)], lme_table$feature[which(lme_table$diffy_diffm.p.value<0.05)])


#-----------------------------------------------------------------------------------------------------#
#							redo heatmap pro_diffy
#-----------------------------------------------------------------------------------------------------#
load(file = paste0(s_OUT_dir,"Scorpius/Scorpius_img.Rdata"))
draw_trajectory_heatmap(expression_DMSO[,Genes_of_interest[["pro_diffy"]] ], traj_DMSO$time, interaction(pheno_DMSO $Corisol,pheno_DMSO $stage),modules[match(Genes_of_interest[["pro_diffy"]],modules$feature),])


# Reduce dimention here, check how many components contain a lot of info, then push into infer algorithm
# PCA
pca = prcomp(t(logcpm[Genes_of_interest[["pro_diffy"]],]),scale=TRUE)

# make autoplot of PCA to show variance explained on first 2 comp
autoplot(pca, col = as.numeric(Pheno$Corisol)+1)


#-----------------------------------------------------------------------------------------------------#
#							redo heatmap diffy_diffm
#-----------------------------------------------------------------------------------------------------#
draw_trajectory_heatmap(expression_DMSO[,Genes_of_interest[["diffy_diffm"]] ], traj_DMSO$time, interaction(pheno_DMSO $Corisol,pheno_DMSO $stage),modules[match(Genes_of_interest[["diffy_diffm"]],modules$feature),])


# Reduce dimention here, check how many components contain a lot of info, then push into infer algorithm
# PCA
pca = prcomp(t(logcpm[Genes_of_interest[["diffy_diffm"]],]),scale=TRUE)

# make autoplot of PCA to show variance explained on first 2 comp
autoplot(pca, col = as.numeric(Pheno$Corisol)+1)

#-----------------------------------------------------------------------------------------------------#
#							redo heatmap OVERLAP
#-----------------------------------------------------------------------------------------------------#

draw_trajectory_heatmap(expression_DMSO[,Genes_of_interest[["intersect"]] ], traj_DMSO$time, interaction(pheno_DMSO $Corisol,pheno_DMSO $stage),modules[match(Genes_of_interest[["intersect"]],modules$feature),])


# Reduce dimention here, check how many components contain a lot of info, then push into infer algorithm
# PCA
pca = prcomp(t(logcpm[Genes_of_interest[["intersect"]],]),scale=TRUE)

# make autoplot of PCA to show variance explained on first 2 comp
autoplot(pca, col = as.numeric(Pheno$Corisol)+1)

#-----------------------------------------------------------------------------------------------------#
#							Check the modules; see if genes enriched in some modules
#-----------------------------------------------------------------------------------------------------#
Genes_of_interest[["pro_diffy"]]
Genes_of_interest[["diffy_difm"]]

pro_diffy_mods = data.frame(Scorpius_obj$modules[match(Genes_of_interest[["pro_diffy"]],Scorpius_obj$modules$feature),])
diffy_difm_mods = data.frame(Scorpius_obj$modules[match(Genes_of_interest[["diffy_diffm"]],Scorpius_obj$modules$feature),])

tablesizes = table(Scorpius_obj$modules$module)

pro_diffy_mods_perc = round(table(pro_diffy_mods$module) / tablesizes,2)
#   1    2    3    4    (5)    6    7    8 
#0.27 0.13 0.14 0.07 0.29 0.20 0.11 0.29 
diffy_difm_mods_perc = round(table(diffy_difm_mods$module) / tablesizes,2)
#   1    2    3    4    (5)    6    7    8 
#0.13 0.17 0.07 0.18 0.71 0.36 0.16 0.25
# mod 5 is huge!

#-----------------------------------------------------------------------------------------------------#
#							generate tables
#-----------------------------------------------------------------------------------------------------#

write.table(format(lme_table, digits=3),paste0(s_ROOT,s_out_folder,".FnT/lme_table.tsv"),col.names = TRUE,row.names = TRUE,sep = "\t",quote = FALSE)


#-----------------------------------------------------------------------------------------------------#
#							output
#-----------------------------------------------------------------------------------------------------#
save.image(file = paste0(s_OUT_dir,"LME/LME_img.Rdata"))
save(lme_table,file = paste0(s_OUT_dir,"LME/lme_table.Rdata"))
save(Genes_of_interest,file = paste0(s_OUT_dir,"LME/Genes_of_interest.Rdata"))


#-----------------------------------------------------------------------------------------------------#
#							old
#-----------------------------------------------------------------------------------------------------#

#
##lme_table$p.adj = p.adjust(lme_table$pval)
#
## get genes significant in pro-diffy OR diffy-diffm
#sign_names = lme_table[which(lme_table$pro_diffy.p.value<0.05 | lme_table$diffy_diffm.p.value <0.05),1]
#length(sign_names) # 139 genes affected
#
#
## Get the FC and sign from the DE
#load(file = paste0(s_OUT_dir,"DE/DE_res_all.Rdata"))
#
## extract gene names following sign cutoff and FC cutoff
#genes_dif_cor_diffm = Extract_DEGs(DE_res_all$DIFFM[sign_names,]) # 2 hits; FOXJ1, RIPPLY2
#genes_dif_cor_diffy = Extract_DEGs(DE_res_all$DIFFY[sign_names,]) # 7 hits; CDK6-AS1, MT-TG, ABI3BP, TTC36-AS1, C1QTNF1, JAKMIP1-DT, RAB6C-AS1
#
#
#
#
#plot(t(logcpm[c(genes_dif_cor_diffm,genes_dif_cor_diffy),])~interaction(Pheno$stage,Pheno$Corisol))
#
#
#load(file = paste0(s_OUT_dir,"Scorpius/Scorpius_img.Rdata"))
#
#modules_diff = as.data.frame(modules)
#modules_diff =  modules_diff[modules_diff[,1]%in%sign_names,]
#SCORPIUS::draw_trajectory_heatmap(expression_DMSO[,sign_names], traj_DMSO$time, interaction(pheno_DMSO $Corisol,pheno_DMSO $stage), modules_diff)
#
#expression_DMSO[,sign_names]