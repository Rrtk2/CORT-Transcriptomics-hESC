#-----------------------------------------------------------------------------------------------------#
# 							GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#
# File description:
#	Name
#		Functions.R
#
#	Purpose 
#		This code was made as sub-master level regulation script.
#
# Author comment:
#	Rick A. Reijnders 
#	ra.reijnders@maastrichtuniversity.nl
#
#	Personal comment:
#
#
#-----------------------------------------------------------------------------------------------------#
#							Main settings or load
#-----------------------------------------------------------------------------------------------------#




#-----------------------------------------------------------------------------------------------------#
#							Main Functions
#-----------------------------------------------------------------------------------------------------#
normalizeMedianValues <- function(x) 
  # Jacked from limmma
#	Normalize columns of a matrix to have the same median value
#	Gordon Smyth
#	24 Jan 2011.  Last modified 24 Jan 2011.
{
	narrays <- NCOL(x)
	if(narrays==1) return(x)
	cmed <- log(apply(x, 2, median, na.rm=TRUE))
	cmed <- exp(cmed - mean(cmed))
	t(t(x)/cmed)
}

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
   
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
   
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

Extract_DEGs = function(DE_Obj=NA, .s_pvalTH = s_pvalTH, .s_fcTH = s_fcTH){
	rownames(DE_Obj[DE_Obj$P.Value<.s_pvalTH & abs(DE_Obj$logFC) > log2(s_fcTH),])
}


#samples as cols, features as rows
TMMnormalize = function(obs, plotting=TRUE,p = 0.75){
  # Jacked from edgeR, modified slightly with NO impact on the internal mechanism (TMM).
	# Just made a standalone function as im often searching for the normalized table. Takes in RAW counts; outputs corrected(TMM normalized) counts.
	
	#Sources:
	# https://rdrr.io/bioc/edgeR/src/R/cpm.R
	# https://rdrr.io/bioc/edgeR/src/R/calcNormFactors.R#sym-.calcFactorTMM

	lib.size = apply(obs,2,sum)#rep(1,dim(obs)[2])
	if(plotting) plot(lib.size,main="Lib sizes")


	# calculate the 75th quantile of each sample
		f <- rep_len(1,ncol(obs))
		for (j in seq_len(ncol(obs))) f[j] <- quantile(obs[,j], probs=p)
		if(plotting) plot(f,main="Q75 of all samples in RAW")
		if(min(f)==0) warning("One or more quantiles are zero")
		f75 = f / lib.size

	# if median quantile of samples is extremely small, do
		# find max summed sqrt column -> (sum(sqrt(col)) per col, use this as ref sample
	# else
		# Find the sample which the 75th quantile is equal to mean.	Assuming this sample is representable (as it is 75th quantile)
	if(median(f75)<1e-20){
		refColumn <- which.max(colSums(sqrt(x)))
					} else {
		refColumn <- which.min(abs(f75-mean(f75)))}
		if(plotting) print(refColumn)
		ref = obs[,refColumn]
		
	# TMM
	res_f = c()
	obs2 = obs
	for(i in 1:(dim(obs)[2])){

							
		libsize.temp_obs=NULL
		libsize.ref=NULL
		logratioTrim=.3
		sumTrim=0.05
		doWeighting=TRUE
		Acutoff=-1e10

		temp_obs <- as.numeric(obs[,i])
		ref <- as.numeric(ref)

		if( is.null(libsize.temp_obs) ) nO <- sum(temp_obs) else nO <- libsize.temp_obs
		if( is.null(libsize.ref) ) nR <- sum(ref) else nR <- libsize.ref

		logR <- log2((temp_obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
		absE <- (log2(temp_obs/nO) + log2(ref/nR))/2  # absolute expression
		v <- (nO-temp_obs)/nO/temp_obs + (nR-ref)/nR/ref   # estimated asymptotic variance

	#	remove infinite values, cutoff based on A
		fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)

		logR <- logR[fin]
		absE <- absE[fin]
		v <- v[fin]

		if(max(abs(logR)) < 1e-6){
			res_f[i] = 1
			obs2[,i] = (obs[,i]/(sum(obs[,i])*res_f[i])) * sum(ref)
			next
		}

	#	taken from the original mean() function
		n <- length(logR)
		loL <- floor(n * logratioTrim) + 1
		hiL <- n + 1 - loL
		loS <- floor(n * sumTrim) + 1
		hiS <- n + 1 - loS

	#	keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
	#	a fix from leonardo ivan almonacid cardenas, since rank() can return
	#	non-integer values when there are a lot of ties
		keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)

		if(doWeighting){
			f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
		}else{
			f <- mean(logR[keep], na.rm=TRUE)
		}
	#	Results will be missing if the two libraries share no features with positive counts
	#	In this case, return unity
		if(is.na(f)) f <- 0
		
		res_f[i] = 2^f
		
		#obs2[,i] = (obs[,i]/(sum(obs[,i])*res_f[i])) * 1e6 # adjusted for factor
		obs2[,i] = (obs[,i]/(sum(obs[,i])*res_f[i])) * sum(ref)	
	}
	# Show plots
	if(plotting) boxplot(log2(obs),main="log2(x+1) before TMM normalization")
	if(plotting) boxplot(log2(obs2),main="log2(x+1) after TMM normalization")
	
	# Return object
	return(obs2)
}

#-----------------------------------------------------------------------------------------------------#
#							PLOTTING 
#-----------------------------------------------------------------------------------------------------#

# Auto generate high res images. PDF and TIFF at location given below
# default is 7 inch; this is now 480 (px?) in this function (so both pdf and tidd scale nice)
Rplot = function(insert=NA,title="Temp_title",resolution = 350, width = 480, height = 480){
	pdf(paste0(s_ROOT,s_out_folder,".FnT/",title,".pdf"),width = width/(480/7),height = height/(480/7)) 
		print({insert})
	dev.off()


	tiff(paste0(s_ROOT,s_out_folder,".FnT/",title,".tiff"),width = width*(resolution/72), height = height*(resolution/72),res = resolution)
		print({insert})
	dev.off()
}

# plot2 is used
# Auto generate high res images. PDF and TIFF at location given below
# default is 7 inch; this is now 480 (px?) in this function (so both pdf and tidd scale nice)
Rplot2 = function(insert=NA,title="Temp_title",resolution = 350, width = 480, height = 480,s_figure_folder = "C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Data_QC/test/.FnT/Boxplots/"){
	pdf(paste0(s_figure_folder,title,".pdf"),width = width/(480/7),height = height/(480/7)) 
		print({insert})
	dev.off()


	tiff(paste0(s_figure_folder,title,".tiff"),width = width*(resolution/72), height = height*(resolution/72),res = resolution)
		print({insert})
	dev.off()
}


Rvolcano = function(data = NA, title = "", s_pvalTH = 0.05,s_fcTH = 1.5, showtopgenesnumber = 5){

	# Load packages
	require(tidyverse)
	require(ggrepel)
	require(clusterProfiler)
	require("org.Hs.eg.db")


	##################################
	#     Data formatting
	##################################
	col_nonsig = "grey70"
	col_sign= "grey30"
	col_up= "red2"
	col_down= "blue2"

	# define cols (sign)
	Cols = rep("black",length(data$P.Value))
	Cols[data$P.Value>s_pvalTH] = col_nonsig#"Non. Sign."#"grey20" #
	Cols[data$P.Value<s_pvalTH] = col_sign#"Sign."#"grey80" #
	Cols[data$P.Value<s_pvalTH & data$logFC>=log2(s_fcTH)] = col_up#"Sign. Up"#"red2" #"Sign. Up"
	Cols[data$P.Value<s_pvalTH & data$logFC<= -log2(s_fcTH)] = col_down#"Sign. Down"#"blue2" #"Sign. Down"

	# define sizes (sign)
	sizes = rep(0.1,length(data$P.Value))
	sizes[data$P.Value>s_pvalTH] = 1.5
	sizes[data$P.Value<s_pvalTH] = 1
	sizes[data$P.Value<s_pvalTH & data$logFC>=log2(s_fcTH)] = 2
	sizes[data$P.Value<s_pvalTH & data$logFC<= -log2(s_fcTH)] = 2


	# define top genes (pval & FCs_fcTH)
	top_genes = data[data$P.Value<s_pvalTH & abs(data$logFC) >= log2(s_fcTH),]

	  
	##################################
	#     Plotting volcano
	##################################
	# basic volcano
	p3 = ggplot(data, aes(x= logFC, y= -log(P.Value,10))) + #color = Cols #,
	  geom_point(size = sizes,colour = Cols) #+
	  #scale_color_manual(values = c("Sign." = "grey50", 
	   #                             "Non. Sign." = "grey70",
		#                            "Sign. Up" = "red2",
		 #                           "Sign. Down" = "blue2")) 
	#+  guides(size=s_fcTH)
	#scale_fill_manual(values=c("red", "blue", "green","blue"))
	p3
	
	# text at informative ones (36)
	p4 =  p3 +
	  geom_hline(yintercept  = -log10(s_pvalTH), linetype="dashed")+ #sign cutoff
	  geom_vline(xintercept  =  log2(s_fcTH), linetype="dashed",col=alpha("red2",0.3))+ # + logFC cutoff
	  geom_vline(xintercept  = - log2(s_fcTH), linetype="dashed",col=alpha("blue2",0.3))+
	  geom_label_repel(data = top_genes,
					   mapping = aes(logFC, -log(P.Value,10), label = c(rownames(top_genes)[1:showtopgenesnumber],rep(NA,length(rownames(top_genes))-showtopgenesnumber))),
					   size = 3,
					   max.overlaps = 100,
					   box.padding = 0.4,
					   label.padding = 0.3,
					   force=2,
					   label.size=0,
					   segment.size = 0.1)
	  p4
	p5 = p4 +
		xlab(expression("log"[2]*"FC")) + 
		ylab(expression("-log"[10]*"P.Value")) +
		ggtitle(paste0("Volcanoplot ",title)) # - logFC cutoff
	p5 

}
#-----------------------------------------------------------------------------------------------------#
#							GO functions
#-----------------------------------------------------------------------------------------------------#
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




#-----------------------------------------------------------------------------------------------------#
#							Cleanup
#-----------------------------------------------------------------------------------------------------#
Rclean = function(){
	rm(list=ls(envir = .GlobalEnv)[grep(ls(envir = .GlobalEnv),pattern = "^temp_")],envir = .GlobalEnv)
}

#-----------------------------------------------------------------------------------------------------#
#							output
#-----------------------------------------------------------------------------------------------------#

