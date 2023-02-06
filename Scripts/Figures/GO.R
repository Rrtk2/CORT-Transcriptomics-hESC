
#-----------------------------------------------------------------------------------------------------#
# 							GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#
# File description:
#	Name
#		Coloc.R
#
#	Purpose 
# 		Dedicated script to generate all functional annotations (GO and KEGG)
#
# Author comment:
#	Rick A. Reijnders 
#	ra.reijnders@maastrichtuniversity.nl
#
#	Personal comment:
#	!  This script now gerates all images; maybe more appropriate in "Functional annotation" block; then
#   only select good ones here! 
#
#-----------------------------------------------------------------------------------------------------#
#							Main settings or load
#-----------------------------------------------------------------------------------------------------#
# Auto generate high res images. PDF and TIFF at location given below
# default is 7 inch; this is now 480 (px?) in this function (so both pdf and tidd scale nice)
Rplot = function(insert=NA,title="Temp_title",resolution = 350, width = 480, height = 480){
	pdf(paste0("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Data_QC/Figures/",title,".pdf"),width = width/(480/7),height = height/(480/7)) 
		print({insert})
	dev.off()


	tiff(paste0("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Data_QC/Figures/",title,".tiff"),width = width*(resolution/72), height = height*(resolution/72),res = resolution)
		print({insert})
	dev.off()
}


#-----------------------------------------------------------------------------------------------------#
#							DAta
#-----------------------------------------------------------------------------------------------------#
load("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Data_QC/DE/DE_res_all.Rdata")

data  = DE_res_all$pro_diffy_cortVSpro_diffy_dmso

# needs SvR DE data
#load("~/GitHub/PRISMO_Transcriptomics/Scripts/POST only processing/data_DE.RData")

# Load packages
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library("org.Hs.eg.db")
library("enrichplot")

Plot_GSE_GO = function(GSE_obj = NA, title = "temp_title"){
	# title = go_ont
	if(class(GSE_obj)[1]!="gseaResult"){warning("Plotting failed, no gseaResult object.\n");break()}
	Rplot(cnetplot(GSE_obj),paste0("GSE_",title,"_cnetplot"))
	
	Rplot(dotplot(GSE_obj),paste0("GSE_",title,"_dotplot"))
	
	Rplot(plotGOgraph(GSE_obj),paste0("GSE_",title,"_plotGOgraph"))
	
	Rplot(upsetplot(GSE_obj),paste0("GSE_",title,"_upsetplot"))
	
	Rplot(emapplot(pairwise_termsim(GSE_obj)),paste0("GSE_",title,"_pairwise_emapplot"))
	
	Rplot(gseaplot(GSE_obj, title = GSE_obj$Description[1], geneSetID = 1),paste0("GSE_",title,"_gseaplot1"))
	
	#Rplot(gseaplot(GSE_obj, title = GSE_obj$Description[2], geneSetID = 2),paste0("GSE_",title,"_gseaplot2"))
	
	#Rplot(gseaplot(GSE_obj, title = GSE_obj$Description[3], geneSetID = 3),paste0("GSE_",title,"_gseaplot3"))
}

Plot_GSE_KEGG = function(GSE_obj = NA, title = "temp_title"){
	# title = go_ont
	if(class(GSE_obj)[1]!="gseaResult"){warning("Plotting failed, no gseaResult object.\n");break()}
	Rplot(cnetplot(GSE_obj),paste0("GSE_",title,"_cnetplot"))
	
	Rplot(dotplot(GSE_obj),paste0("GSE_",title,"_dotplot"))
	
	Rplot(upsetplot(GSE_obj),paste0("GSE_",title,"_upsetplot"))
	
	Rplot(emapplot(pairwise_termsim(GSE_obj)),paste0("GSE_",title,"_pairwise_emapplot"))
	
	Rplot(gseaplot(GSE_obj, title = GSE_obj$Description[1], geneSetID = 1),paste0("GSE_",title,"_gseaplot1"))
	
	#Rplot(gseaplot(GSE_obj, title = GSE_obj$Description[2], geneSetID = 2),paste0("GSE_",title,"_gseaplot2"))
	
	#Rplot(gseaplot(GSE_obj, title = GSE_obj$Description[3], geneSetID = 3),paste0("GSE_",title,"_gseaplot3"))
}


Plot_Enrich_GO = function(GSE_obj = NA, title = "temp_title"){
	# title = go_ont
	if(class(GSE_obj)[1]!="enrichResult"){stop("Plotting failed, no enrichResult object.\n")}
	Rplot(cnetplot(GSE_obj),paste0("ENRICH_",title,"_cnetplot"))
	
	Rplot(barplot(GSE_obj),paste0("ENRICH_",title,"_barplot"))
	
	Rplot(dotplot(GSE_obj),paste0("ENRICH_",title,"_dotplot"))
	
	Rplot(plotGOgraph(GSE_obj),paste0("ENRICH_",title,"_plotGOgraph"))
	
	Rplot(upsetplot(GSE_obj),paste0("ENRICH_",title,"_upsetplot"))
	
	Rplot(emapplot(pairwise_termsim(GSE_obj)),paste0("ENRICH_",title,"_pairwise_emapplot"))

}
	

Plot_Enrich_KEGG = function(GSE_obj = NA, title = "temp_title"){
	# title = go_ont
	if(class(GSE_obj)[1]!="enrichResult"){stop("Plotting failed, no enrichResult object.\n")}
	Rplot(cnetplot(GSE_obj),paste0("ENRICH_",title,"_cnetplot"))
	
	Rplot(barplot(GSE_obj),paste0("ENRICH_",title,"_barplot"))
	
	Rplot(dotplot(GSE_obj),paste0("ENRICH_",title,"_dotplot"))
	
	Rplot(upsetplot(GSE_obj),paste0("ENRICH_",title,"_upsetplot"))
	
	Rplot(emapplot(pairwise_termsim(GSE_obj)),paste0("ENRICH_",title,"_pairwise_emapplot"))

}	

kegg_prep = function(Obj){
	# this function converts SYMBOL to ENTREZID
	
	# Convert gene IDs for gseKEGG function
	# We will lose some genes here because not all IDs will be converted
	ids<-bitr(rownames(Obj), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
	 # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
	dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

	# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
	df2 = Obj[rownames(Obj) %in% dedup_ids$ENSEMBL,]

	# Create a new column in df2 with the corresponding ENTREZ IDs
	df2$Y = dedup_ids$ENTREZID

	# Create a vector of the gene unuiverse
	kegg_gene_list <- df2$logFC

	# Name vector with ENTREZ ids
	names(kegg_gene_list) <- df2$Y

	# omit any NA values 
	kegg_gene_list<-na.omit(kegg_gene_list)

	# sort the list in decreasing order (required for clusterProfiler)
	kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
	
	return(kegg_gene_list)
}
		
		
	

# Main function that will do it all
FuncAnno = function(DE_obj = DE_Results, .s_pvalTH = s_pvalTH, .s_fcTH = s_fcTH, .s_ROOT_dir = s_ROOT_dir,.s_out_folder = s_out_folder){
	# define top genes (pval & FC1.5)
	top_genes = DE_obj[DE_obj$P.Value<.s_pvalTH & abs(DE_obj$logFC) >= log2(.s_fcTH),]

	geneList_all = DE_obj
	geneList_top = top_genes

	for(go_ont in c("BP", "MF", "CC")){
		##################################
		#     GSE GO
		##################################

		# overrep (compares full data)
		egogenelist = (geneList_all$logFC)
		names(egogenelist) = rownames(geneList_all)
		egogenelist = sort(egogenelist,decreasing = T)   
		ego3 <- gseGO(geneList     = egogenelist, # needs: order ranked geneList
					  OrgDb        = org.Hs.eg.db,
					  ont          = go_ont,
					  minGSSize    = 10,
					  maxGSSize    = 100,
					 # scoreType = "pos",
					  #pvalueCutoff = 0.001,
					  #readable      = TRUE,
					  pAdjustMethod = "fdr",
					  pvalueCutoff  = s_pvalTH,
					  keyType = "ENSEMBL")
		head(ego3)
		
		if(dim(ego3)[1]!=0){
			ego4 <- simplify(ego3, cutoff=0.7, by="p.adjust", select_fun=min)

			Plot_GSE_GO(GSE_obj = ego4, title = paste0("GO_",go_ont))
		}else{message("No GO gse")}
		
		##################################
		#     ENRICH GO
		##################################

		# enrichment (compares 36 to all)
		ego <- enrichGO(gene          = rownames(geneList_top),
						universe      = rownames(geneList_all),
						OrgDb         = org.Hs.eg.db,
						ont           = go_ont,
						minGSSize = 10,
						maxGSSize = 100,
						#pvalueCutoff  = .s_pvalTH,
						pAdjustMethod = "fdr",
						pvalueCutoff  = s_pvalTH,
						#readable      = TRUE,
						keyType = "ENSEMBL")
		
		if(dim(ego)[1]!=0){		
			ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
			
			Plot_Enrich_GO(GSE_obj = ego2, title = paste0("GO_",go_ont))
		}else{message("No GO enrich")}
	}
	
	##################################
	#     GSE KEGG
	##################################
	
	#R.utils::setOption("clusterProfiler.download.method","auto") # if species is not found (should work with "hsa")
	kk2 <- gseKEGG(geneList     = kegg_prep(geneList_all),
		   organism     = "hsa",
		   nPerm        = 10000,
		   minGSSize    = 3,
		   maxGSSize    = 800,
		   pvalueCutoff = 0.05,
		   pAdjustMethod = "none",
		   keyType       = "ncbi-geneid")
		   
	#kk3<- simplify(kk2, cutoff=0.7, by="p.adjust", select_fun=min)
	if(dim(kk2)[1]!=0){
		Plot_GSE_KEGG(GSE_obj = kk2, title = paste0("KEGG"))
	}else{message("No kegg gse")}
	
	##################################
	#     ENRICH KEGG
	##################################
	
	R.utils::setOption("clusterProfiler.download.method","auto") # if species is not found (should work with "hsa")
	ekegg <- enrichKEGG(gene          = names(kegg_prep(geneList_top)),
						universe      = names(kegg_prep(geneList_all)),
						organism         = "hsa",
						pAdjustMethod = "none",
						minGSSize = 10,
						maxGSSize = 100,
						pvalueCutoff  = .s_pvalTH,
						keyType = "ncbi-geneid")
		   
	if(dim(ekegg)[1]!=0){
		Plot_Enrich_KEGG(GSE_obj = ekegg, title = paste0("KEGG"))
	}else{message("No kegg enrich")}
}

# Run the massive function. Takes a minute
# make IDs to SYMBOL; as these are requered
#ids<-bitr(rownames(data), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb=org.Hs.eg.db)
#data2 = data
#
#Newnames = ids[match(rownames(data),ids$ENSEMBL),"SYMBOL"]
#
#
#rownames(data2)[!is.na(Newnames)][!duplicated(Newnames[!is.na(Newnames)])] = Newnames[!is.na(Newnames)][!duplicated(Newnames[!is.na(Newnames)])]

FuncAnno(DE_obj = data)


# Full Cleanup
#rm(list=ls())

