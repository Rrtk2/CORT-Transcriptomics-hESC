#-----------------------------------------------------------------------------------------------------#
#							Call settings
#-----------------------------------------------------------------------------------------------------#
source("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Settings.R")
#-----------------------------------------------------------------------------------------------------#
#							libs
#-----------------------------------------------------------------------------------------------------#
library("edgeR")
#-----------------------------------------------------------------------------------------------------#
#							Get data from LIBD Stem Cell Browser (http://stemcell.libd.org/scb/data_links.html)
#-----------------------------------------------------------------------------------------------------#
load("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Data_RAW/libd_stemcell_timecourse_rseGene_n157.rda")
# results rse_gene


#-----------------------------------------------------------------------------------------------------#
#							Input
#-----------------------------------------------------------------------------------------------------#
load(file = paste0(s_OUT_dir,"Prep/logcpm.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/Pheno.Rdata"))
load(file = paste0(s_OUT_dir,"Scorpius/Scorpius_obj.Rdata"))
load(file = paste0(s_OUT_dir,"LME/Genes_of_interest.Rdata"))


#-----------------------------------------------------------------------------------------------------#
#							Prepare LIBD Stem Cell data
#-----------------------------------------------------------------------------------------------------#
# get the count data and normalize to TMM log2(cpm+1)
rse_data = assay(rse_gene)
y = DGEList(counts=rse_data)

#keep= filterByExpr(y,min.count=5)

#y = y[keep, , keep.lib.sizes=FALSE]

y = calcNormFactors(y, method = "TMM")

rse_data_norm = edgeR::cpm(y, log=TRUE)

# get relevant days
#rse_data_norm_relevant = rse_data_norm[,rse_gene$DAY %in% c(15, 21,49 )]
#rse_gene_adjust = rse_gene[,rse_gene$DAY %in% c(15, 21,49 )]


# get relevant genes
rse_gene_relevant = rse_data_norm[gsub(rownames(rse_data_norm),pattern = "\\..$|\\...$",replacement = "") %in% paste0(c(Genes_of_interest$pro_diffy,Genes_of_interest$diffy_diffm)),]


# get the logcounts 
FoundGenes = unique(gsub(rownames(rse_gene_relevant),pattern = "\\..$|\\...$",replacement = ""))
logcpm_relevant  = logcpm[rownames(logcpm)%in%FoundGenes,]

# logcpm only DMSO
logcpm_relevant = logcpm_relevant[,grep(colnames(logcpm_relevant),pattern = "^D_")]

# Reorder rse and logcpm based on foundgenes
logcpm_relevant = logcpm_relevant[FoundGenes,]
rse_gene_relevant = rse_gene_relevant[match(FoundGenes,(gsub(rownames(rse_gene_relevant),pattern = "\\..$|\\...$",replacement = ""))),]

# Take averages of the timings
rse_df = data.frame(Day15 = apply(rse_gene_relevant[,rse_gene$DAY=="15"],1,mean),
					Day21 = apply(rse_gene_relevant[,rse_gene$DAY=="21"],1,mean),
					Day49 = apply(rse_gene_relevant[,rse_gene$DAY=="49"],1,mean))

logcounts_df = data.frame(Pro = apply(logcpm_relevant[,grep(colnames(logcpm_relevant),pattern = "pro")],1,mean),
					Diffy = apply(logcpm_relevant[,grep(colnames(logcpm_relevant),pattern = "diffy")],1,mean),
					Diffm = apply(logcpm_relevant[,grep(colnames(logcpm_relevant),pattern = "diffm")],1,mean))



corframes  = function(X,Y,method = "p"){

	out = data.frame(Feature = rownames(X), cor = NA)
	for(i in 1:dim(X)[1]){
		res = cor(as.numeric(X[i,]),as.numeric(Y[i,]),method = method)
		out[i,2] = c(res)
	}
	
	return(out)
	
}

res_cor = corframes(rse_df,logcounts_df)
hist(res_cor$cor)


which(res_cor$Feature=="ENSG00000146006.7")

#-----------------------------------------------------------------------------------------------------#
#							Extract the relevant Genes
#-----------------------------------------------------------------------------------------------------#
#LRRTM2, KCND3, KCNIP4 GRIA3 TSPAN5 
Getgenes = c("ENSG00000146006.7","ENSG00000171385.9","ENSG00000185774.14","ENSG00000125675.17","ENSG00000168785.7")

res_cor[match(Getgenes,res_cor$Feature),]
 
#               Feature       cor
#16  ENSG00000146006.7 0.9717350
#4   ENSG00000171385.9 0.9690624
#11 ENSG00000185774.14 0.8780239
#71 ENSG00000125675.17 0.9532571
#14  ENSG00000168785.7 0.9650923


table(round(abs(res_cor$cor),1))
#  0 0.2 0.4 0.5 0.6 0.7 0.8 0.9   1 
#  1   1   1   3   1   3  10  21  30 

30+21+10 / 30+21+10+3+1+3+1+1+1
# 92.33333



