#-----------------------------------------------------------------------------------------------------#
#							Call settings
#-----------------------------------------------------------------------------------------------------#
source("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Settings.R")


#-----------------------------------------------------------------------------------------------------#
#							Input
#-----------------------------------------------------------------------------------------------------#
load(file = paste0(s_OUT_dir,"Prep/Pheno.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/Data.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/logcpm.Rdata"))


#-----------------------------------------------------------------------------------------------------#
#							main
#-----------------------------------------------------------------------------------------------------#


# make design matrix
temp_design <- model.matrix(~ 0 + Corisol:stage,Pheno)
colnames(temp_design) = make.names(colnames(temp_design))

# Make LMmodel and mod t test
temp_fit <- lmFit(logcpm,temp_design)

# Show the coefs for first gene
data.frame(temp_fit$coef[1,])

temp_contrast_matrix <- makeContrasts(
CORT = (CorisolTRUE.stagediffy  + CorisolTRUE.stagepro + CorisolTRUE.stagediffm) - (CorisolFALSE.stagepro + CorisolFALSE.stagediffm + CorisolFALSE.stagediffy),
PRO = (CorisolTRUE.stagepro) - ( CorisolFALSE.stagepro ),
DIFFY = (CorisolTRUE.stagediffy) - ( CorisolFALSE.stagediffy ),
DIFFM = (CorisolTRUE.stagediffm) - ( CorisolFALSE.stagediffm ),
pro_diffy_dmso =  (CorisolFALSE.stagediffy - CorisolFALSE.stagepro),
diffy_diffm_dmso =  (CorisolFALSE.stagediffm - CorisolFALSE.stagediffy),
pro_diffm_dmso =  (CorisolFALSE.stagediffm - CorisolFALSE.stagepro),
pro_diffy_cort =  (CorisolTRUE.stagediffy - CorisolTRUE.stagepro),
diffy_diffm_cort =  (CorisolTRUE.stagediffm - CorisolTRUE.stagediffy),
pro_diffm_cort =  (CorisolTRUE.stagediffm - CorisolTRUE.stagepro),
pro_diffy_cortVSpro_diffy_dmso = (CorisolTRUE.stagediffy - CorisolTRUE.stagepro) - (CorisolFALSE.stagediffy - CorisolFALSE.stagepro),
diffy_diffm_cortVSdiffy_diffm_dmso = (CorisolTRUE.stagediffm - CorisolTRUE.stagediffy) - (CorisolFALSE.stagediffm - CorisolFALSE.stagediffy),
pro_diffm_cortVSpro_diffm_dmso = (CorisolTRUE.stagediffm - CorisolTRUE.stagepro) - (CorisolFALSE.stagediffm - CorisolFALSE.stagepro),
levels = colnames(temp_design))
 
 
 
fit2 <- contrasts.fit(temp_fit, temp_contrast_matrix)

fit2 <- eBayes(fit2,trend=TRUE)


# save all conte
DE_res_all = lapply(colnames(temp_contrast_matrix),function(i){topTable(fit2, coef = i, num=Inf, sort.by = "P", adjust="BH")})
names(DE_res_all ) = colnames(temp_contrast_matrix)

DE_res = lapply(colnames(temp_contrast_matrix),function(i){topTable(fit2, coef = i, num=Inf, sort.by = "P", adjust="BH",p.value = 0.05)})
names(DE_res ) = colnames(temp_contrast_matrix)


DE_res2 = lapply(colnames(temp_contrast_matrix),function(i){topTable(fit2, coef = i, num=Inf, sort.by = "P", adjust="BH",p.value = 0.00005)})
names(DE_res2 ) = colnames(temp_contrast_matrix)



DE_res2$Bestgenes = unique(c(rownames(DE_res2[["PRO"]]),rownames(DE_res2[["DIFFY"]]),rownames(DE_res2[["DIFFM"]])))

#-----------------------------------------------------------------------------------------------------#
#							Make plot
#-----------------------------------------------------------------------------------------------------#

Rvolcano(data = DE_res_all$PRO, title = "condition effects in PRO", s_pvalTH = s_pvalTH, s_fcTH = s_fcTH, showtopgenesnumber = 5)
Rvolcano(data = DE_res_all$DIFFY, title = "condition effects in DIFFY", s_pvalTH = s_pvalTH, s_fcTH = s_fcTH, showtopgenesnumber = 5)
Rvolcano(data = DE_res_all$DIFFM, title = "condition effects in DIFFM", s_pvalTH = s_pvalTH, s_fcTH = s_fcTH, showtopgenesnumber = 5)

Rplot(insert={Rvolcano2(data = DE_res_all$PRO, title = "condition effects in PRO", s_pvalTH = s_pvalTH, s_fcTH = s_fcTH)},title="DE_volcano_pro",resolution = 350, width = 480, height = 480)

Rplot(insert={Rvolcano2(data = DE_res_all$DIFFY, title = "condition effects in DIFFY", s_pvalTH = s_pvalTH, s_fcTH = s_fcTH)},title="DE_volcano_diffy",resolution = 350, width = 480, height = 480)

Rplot(insert={Rvolcano2(data = DE_res_all$DIFFM, title = "condition effects in DIFFM", s_pvalTH = s_pvalTH, s_fcTH = s_fcTH)},title="DE_volcano_diffm",resolution = 350, width = 480, height = 480)

#-----------------------------------------------------------------------------------------------------#
#							GO
#-----------------------------------------------------------------------------------------------------#
FuncAnno(DE_obj = DE_res_all$PRO, .s_pvalTH = s_pvalTH, .s_fcTH = s_fcTH, .s_ROOT_dir = s_ROOT,.s_out_folder = s_out_folder)

FuncAnno(DE_obj = DE_res_all$DIFFY, .s_pvalTH = s_pvalTH, .s_fcTH = s_fcTH, .s_ROOT_dir = s_ROOT,.s_out_folder = s_out_folder)

FuncAnno(DE_obj = DE_res_all$DIFFM, .s_pvalTH = s_pvalTH, .s_fcTH = s_fcTH, .s_ROOT_dir = s_ROOT,.s_out_folder = s_out_folder)


#-----------------------------------------------------------------------------------------------------#
#							generate tables
#-----------------------------------------------------------------------------------------------------#

write.table(format(DE_res_all$PRO, digits=3),paste0(s_ROOT,s_out_folder,".FnT/DE_res_all_PRO.tsv"),col.names = TRUE,row.names = TRUE,sep = "\t",quote = FALSE)

write.table(format(DE_res_all$DIFFY, digits=3),paste0(s_ROOT,s_out_folder,".FnT/DE_res_all_DIFFY.tsv"),col.names = TRUE,row.names = TRUE,sep = "\t",quote = FALSE)

write.table(format(DE_res_all$DIFFM, digits=3),paste0(s_ROOT,s_out_folder,".FnT/DE_res_all_DIFFM.tsv"),col.names = TRUE,row.names = TRUE,sep = "\t",quote = FALSE)

#-----------------------------------------------------------------------------------------------------#
#							For STRING
#-----------------------------------------------------------------------------------------------------#
paste0(Extract_DEGs(DE_res_all$PRO),collapse = "@")
paste0(Extract_DEGs(DE_res_all$DIFFY),collapse = "@")
paste0(Extract_DEGs(DE_res_all$DIFFM),collapse = "@")
# then in STRING 
#	- FULL string network
#	- confidence
#	- (minus) texmining
#	- highest confidence (0.9>)
#	- 10 extra in second shell
#	- Hide disconnected nodes

#-----------------------------------------------------------------------------------------------------#
#							PCA
#-----------------------------------------------------------------------------------------------------#
Rplot(insert={
	autoplot(prcomp(t(logcpm)),col=c("black","red3","green4")[as.numeric(Pheno$stage)],shape = c(19,17)[as.numeric(Pheno$Corisol)+1],size = 3)+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black"))

},title="DE_PCA",resolution = 350, width = 480, height = 480)


#-----------------------------------------------------------------------------------------------------#
#							Venn diagramme
#-----------------------------------------------------------------------------------------------------#
x = list( PRO = rownames(DE_res_all$PRO)[DE_res_all$PRO$adj.P.Val<0.05 & 2^DE_res_all$PRO$logFC>s_fcTH],
		Diffy = rownames(DE_res_all$DIFFY)[DE_res_all$DIFFY$adj.P.Val<0.05 & 2^DE_res_all$DIFFY$logFC>s_fcTH], 
		Diffm = rownames(DE_res_all$DIFFM)[DE_res_all$DIFFM$adj.P.Val<0.05 & 2^DE_res_all$DIFFM$logFC>s_fcTH] )



  
#devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
g = ggvenn(
  x, 
  fill_color = c("gray30","#771F1B","#3E7D22"),
  stroke_size = 0.8, set_name_size = 8, show_outside = "none"
  );g
  


#-----------------------------------------------------------------------------------------------------#
#							Save 
#-----------------------------------------------------------------------------------------------------#
save(DE_res,file = paste0(s_OUT_dir,"DE/DE_res.Rdata"))
save(DE_res_all,file = paste0(s_OUT_dir,"DE/DE_res_all.Rdata"))


