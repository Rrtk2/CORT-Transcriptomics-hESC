#-----------------------------------------------------------------------------------------------------#
#							Call settings
#-----------------------------------------------------------------------------------------------------#
source("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Settings.R")

#-----------------------------------------------------------------------------------------------------#
#							Input
#-----------------------------------------------------------------------------------------------------#
load(file = paste0(s_OUT_dir,"DE/logcpm.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/Pheno.Rdata"))
load(file = paste0(s_OUT_dir,"Scorpius/Scorpius_obj.Rdata"))



#-----------------------------------------------------------------------------------------------------#
#							main
#-----------------------------------------------------------------------------------------------------#
#install.packages("gam")
library(gam)
library(ggplot2)
library(ggfortify)

# Fit GAM for each gene using pseudotime as independent variable.
#
#gam.pval <- apply(t(Scorpius_obj$Exp), 1, function(z){
#  d <- data.frame(z=z, t=t)
#  tmp <- gam(z ~ lo(t), data=d)
#  p <- summary(tmp)[4][[1]][1,5]
#  p
#})

t <- round(Scorpius_obj$trajectory$time,2)
# run the correction over DMSO data
gam_fitted <- apply(t(Scorpius_obj$Exp[,Scorpius_obj$modules$feature]), 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  tmp$fitted.values
})
PCA = prcomp(gam_fitted)
autoplot(PCA,col=colorRampPalette(colors = c("red","green"))(length(as.numeric(t)))[rank(as.numeric(t))])

#-----------------------------------------------------------------------------------------------------#
#							Experimental interpoation
#-----------------------------------------------------------------------------------------------------#
n_interpolations = 100

t <- round(Scorpius_obj$trajectory$time,2)
# run the correction over DMSO data
gam_fitted_interpolated <- apply(t(Scorpius_obj$Exp[,Scorpius_obj$modules$feature]), 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  predict(tmp,data.frame(t=(1:n_interpolations/n_interpolations)))
})

# interpolated
#scaled_inter <- scale((gam_fitted_interpolated), center= PCA$center)
Projected_PCA_values <- scale(gam_fitted_interpolated, PCA$center, PCA$scale) %*% PCA$rotation 



colnames ( Projected_PCA_values ) = colnames(PCA$rotation)
rownames ( Projected_PCA_values ) = (1:n_interpolations/n_interpolations)

# fitteed (to show residuals)
Projected_PCA_values2 <- scale(gam_fitted, PCA$center, PCA$scale) %*% PCA$rotation 

colnames ( Projected_PCA_values2 ) = colnames(PCA$rotation)
rownames ( Projected_PCA_values2 ) = rownames(gam_fitted)


# actual! (to show residuals)
Projected_PCA_values3 <- scale(Scorpius_obj$Exp[,Scorpius_obj$modules$feature], PCA$center, PCA$scale) %*% PCA$rotation 

colnames ( Projected_PCA_values3 ) = colnames(PCA$rotation)
rownames ( Projected_PCA_values3 ) = rownames(Scorpius_obj$Exp[,Scorpius_obj$modules$feature])


# plotting RAW, FITTED and INTERPOLATED
ggplot()+
geom_point(PCA$x[,1:2],mapping = aes(x = PCA$x[,1],y = PCA$x[,2]),col=colorRampPalette(colors = c("red","green"))(length(as.numeric(t)))[rank(as.numeric(t))])+
geom_point(Projected_PCA_values[,1:2],mapping = aes(x = Projected_PCA_values[,1],y = Projected_PCA_values[,2]),col=colorRampPalette(colors = c("red","green"))(length(as.numeric(1:n_interpolations/n_interpolations)))[rank(as.numeric(1:n_interpolations/n_interpolations))])+
geom_point(Projected_PCA_values2[,1:2],mapping = aes(x = Projected_PCA_values2[,1],y = Projected_PCA_values2[,2]),col=colorRampPalette(colors = c("black","black"))(length(as.numeric(t)))[rank(as.numeric(t))])+
geom_point(Projected_PCA_values3[,1:2],mapping = aes(x = Projected_PCA_values3[,1],y = Projected_PCA_values3[,2]),col=colorRampPalette(colors = c("blue","blue"))(length(as.numeric(t)))[rank(as.numeric(t))])

# PLOTTING RAW and INTERPOLATED
ggplot()+
geom_point(PCA$x[,1:2],mapping = aes(x = PCA$x[,1],y = PCA$x[,2]),col=colorRampPalette(colors = c("red","green"))(length(as.numeric(t)))[rank(as.numeric(t))])+
geom_point(Projected_PCA_values[,1:2],mapping = aes(x = Projected_PCA_values[,1],y = Projected_PCA_values[,2]),col=colorRampPalette(colors = c("red","green"))(length(as.numeric(1:n_interpolations/n_interpolations)))[rank(as.numeric(1:n_interpolations/n_interpolations))])+
geom_point(Projected_PCA_values3[,1:2],mapping = aes(x = Projected_PCA_values3[,1],y = Projected_PCA_values3[,2]),col=colorRampPalette(colors = c("black","black"))(length(as.numeric(t)))[rank(as.numeric(t))],size=2)

# PLOTTING RAW and INTERPOLATED
ggplot()+
geom_point(PCA$x[,1:2],mapping = aes(x = PCA$x[,1],y = PCA$x[,2]),col=colorRampPalette(colors = c("red","green"))(length(as.numeric(t)))[rank(as.numeric(t))])#+
#geom_point(Projected_PCA_values[,1:2],mapping = aes(x = Projected_PCA_values[,1],y = Projected_PCA_values[,2]),col=colorRampPalette(colors = c("red","green"))(length(as.numeric(1:n_interpolations/n_interpolations)))[rank(as.numeric(1:n_interpolations/n_interpolations))])+
#geom_point(Projected_PCA_values3[,1:2],mapping = aes(x = Projected_PCA_values3[,1],y = Projected_PCA_values3[,2]),col=colorRampPalette(colors = c("black","black"))(length(as.numeric(t)))[rank(as.numeric(t))],size=2)

PCA = prcomp(gam_fitted_interpolated)
autoplot(PCA,col=colorRampPalette(colors = c("red","green"))(length(as.numeric(1:n_interpolations/n_interpolations)))[rank(as.numeric(1:n_interpolations/n_interpolations))])


#-----------------------------------------------------------------------------------------------------#
#							Make plot
#-----------------------------------------------------------------------------------------------------#
@RRR fix please!~!


Rvolcano(data = DE_res_all$PRO, title = "condition effects in PRO", s_pvalTH = s_pvalTH, s_fcTH = s_fcTH, showtopgenesnumber = 5)
Rvolcano(data = DE_res_all$DIFFY, title = "condition effects in DIFFY", s_pvalTH = s_pvalTH, s_fcTH = s_fcTH, showtopgenesnumber = 5)
Rvolcano(data = DE_res_all$DIFFM, title = "condition effects in DIFFM", s_pvalTH = s_pvalTH, s_fcTH = s_fcTH, showtopgenesnumber = 5)

Rplot(insert={Rvolcano(data = DE_res_all$PRO, title = "condition effects in PRO", s_pvalTH = s_pvalTH, s_fcTH = s_fcTH, showtopgenesnumber = 5)},title="DE_volcano_pro",resolution = 350, width = 480, height = 480)

Rplot(insert={Rvolcano(data = DE_res_all$DIFFY, title = "condition effects in DIFFY", s_pvalTH = s_pvalTH, s_fcTH = s_fcTH, showtopgenesnumber = 5)},title="DE_volcano_diffy",resolution = 350, width = 480, height = 480)

Rplot(insert={Rvolcano(data = DE_res_all$DIFFM, title = "condition effects in DIFFM", s_pvalTH = s_pvalTH, s_fcTH = s_fcTH, showtopgenesnumber = 5)},title="DE_volcano_diffm",resolution = 350, width = 480, height = 480)

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
#							output
#-----------------------------------------------------------------------------------------------------#
save(gam_fitted,file = paste0(s_OUT_dir,"GAM/gam_fitted.Rdata"))

