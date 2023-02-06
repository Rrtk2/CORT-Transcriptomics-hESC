#-----------------------------------------------------------------------------------------------------#
#							Call settings
#-----------------------------------------------------------------------------------------------------#
source("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Settings.R")

#-----------------------------------------------------------------------------------------------------#
#							Packages
#-----------------------------------------------------------------------------------------------------#

# Load packages
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library("org.Hs.eg.db")
library("enrichplot")

#-----------------------------------------------------------------------------------------------------#
#							Input
#-----------------------------------------------------------------------------------------------------#
load(file = paste0(s_OUT_dir,"Prep/logcpm.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/Pheno.Rdata"))
load(file = paste0(s_OUT_dir,"Scorpius/Scorpius_obj.Rdata"))

load(file = paste0(s_OUT_dir,"LME/lme_table.Rdata"))
load(file = paste0(s_OUT_dir,"LME/Genes_of_interest.Rdata"))

#-----------------------------------------------------------------------------------------------------#
#							enrich seetings
#-----------------------------------------------------------------------------------------------------#
temp_s_pval = 0.3
temp_s_minGS_size = 3

#-----------------------------------------------------------------------------------------------------#
#							Enrichment pro_diffy
#-----------------------------------------------------------------------------------------------------#
# enrichment (compares selected to all)
ego <- enrichGO(gene          = lme_table[lme_table$pro_diffy.p.value<0.05,"feature"],
			universe      = rownames(logcpm),
			OrgDb         = org.Hs.eg.db,
			ont           = "all",
			minGSSize = temp_s_minGS_size,
			maxGSSize = 100,
			#pvalueCutoff  = .s_pvalTH,
			pvalueCutoff  = temp_s_pval,
			#readable      = TRUE,
			keyType = "ENSEMBL")

						
temp_ncat = min(dim(ego@result)[1],20)		
temp_npanel = length(unique(ego@result$ONTOLOGY))

# titles and such
width = 290 + 240*temp_npanel 
height = 100+ 37 * temp_ncat
resolution = 350


tiff(paste0(s_OUT_dir,"GO/Enrichment_LME_pro_diffy.tiff"),width = width*(resolution/72), height = height*(resolution/72),res = 350) 
	print({dotplot(ego, showCategory = temp_ncat) + facet_grid(~ONTOLOGY)})
dev.off()

write.table(format(ego@result, digits=3),paste0(s_OUT_dir,"GO/Enrichment_LME_pro_diffy.tsv"),col.names = NA,row.names = TRUE,sep = "\t",quote = FALSE)


#--
#-----------------------------------------------------------------------------------------------------#
#							Enrichment diffy_diffm
#-----------------------------------------------------------------------------------------------------#
# enrichment (compares selected to all)
ego <- enrichGO(gene          = lme_table[lme_table$diffy_diffm.p.value<0.05,"feature"],
			universe      = rownames(logcpm),
			OrgDb         = org.Hs.eg.db,
			ont           = "all",
			minGSSize = temp_s_minGS_size,
			maxGSSize = 100,
			#pvalueCutoff  = .s_pvalTH,
			pvalueCutoff  = temp_s_pval,
			#readable      = TRUE,
			keyType = "ENSEMBL")

						
temp_ncat = min(dim(ego@result)[1],20)		
temp_npanel = length(unique(ego@result$ONTOLOGY))

# titles and such
width = 290 + 240*temp_npanel 
height = 100+ 37 * temp_ncat
resolution = 350


tiff(paste0(s_OUT_dir,"GO/Enrichment_LME_diffy_diffm.tiff"),width = width*(resolution/72), height = height*(resolution/72),res = 350) 
	print({dotplot(ego, showCategory = temp_ncat) + facet_grid(~ONTOLOGY)})
dev.off()

write.table(format(ego@result, digits=3),paste0(s_OUT_dir,"GO/Enrichment_LME_diffy_diffm.tsv"),col.names = NA,row.names = TRUE,sep = "\t",quote = FALSE)



#-----------------------------------------------------------------------------------------------------#
#							Enrichment TRAJECTORY genes
#-----------------------------------------------------------------------------------------------------#
# enrichment (compares selected to all)
ego2 <- enrichGO(gene          = lme_table$feature,
			universe      = rownames(logcpm),
			OrgDb         = org.Hs.eg.db,
			ont           = "all",
			minGSSize = temp_s_minGS_size,
			maxGSSize = 100,
			#pvalueCutoff  = .s_pvalTH
			pvalueCutoff  = temp_s_pval,
			#readable      = TRUE,
			keyType = "ENSEMBL")

						
temp_ncat = min(dim(ego2@result)[1],20)		
temp_npanel = length(unique(ego2@result$ONTOLOGY))

# titles and such
width = 290 + 240*temp_npanel 
height = 100+ 37 * temp_ncat
resolution = 350


tiff(paste0(s_OUT_dir,"GO/Enrichment_trajectory.tiff"),width = width*(resolution/72), height = height*(resolution/72),res = 350) 
	print({dotplot(ego2, showCategory = temp_ncat) + facet_grid(~ONTOLOGY)})
dev.off()

write.table(format(ego2@result, digits=3),paste0(s_OUT_dir,"GO/Enrichment_trajectory.tsv"),col.names = NA,row.names = TRUE,sep = "\t",quote = FALSE)



