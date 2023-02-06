#-----------------------------------------------------------------------------------------------------#
#							Call settings
#-----------------------------------------------------------------------------------------------------#
source("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Settings.R")


#-----------------------------------------------------------------------------------------------------#
#							Libraries
#-----------------------------------------------------------------------------------------------------#


#-----------------------------------------------------------------------------------------------------#
#							Input
#-----------------------------------------------------------------------------------------------------#
load(file = paste0(s_OUT_dir,"Scorpius/Scorpius_obj.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/logcpm.Rdata"))
load(file = paste0(s_OUT_dir,"Prep/Pheno.Rdata"))


#-----------------------------------------------------------------------------------------------------#
#							GENIE3
#-----------------------------------------------------------------------------------------------------#
library(GENIE3)

# get the DMSO data only for genes of interes
DMSO_exp = logcpm[Scorpius_obj$modules$feature,Pheno$Corisol==FALSE]

weightMat_old <- GENIE3(DMSO_exp,nCores=8,nTrees = 1000)

# check if edge between one pair of nodes is higher in one, or other direction

weightMat_new = weightMat_old
genelist_row = rownames(weightMat_new)
genelist_col = colnames(weightMat_new)

for(i in genelist_row ){

	for(o in genelist_col){
	
	
		if(weightMat_new[i,o] > weightMat_new[o,i]){
			weightMat_new[o,i] = 0
			}else{
			weightMat_new[i,o] = 0
		}

	}
}



# top 5% of links selected
linkList_old <- getLinkList(weightMat_old,threshold = quantile(weightMat_old,0.99))
# top 5% of links selected
linkList_new <- getLinkList(weightMat_new,threshold = quantile(weightMat_old,0.99))



head(linkList_new)
plot(linkList_new,pch=19,col=alpha("black",0.1))
#  maybe tailor the TARGETS as all genes? see what comes outof it>
# Top 5% to cytoscape

#-----------------------------------------------------------------------------------------------------#
#							try with igraph
#-----------------------------------------------------------------------------------------------------#
library(igraph)
library(RCy3)
library(Rgraphviz)

edge_listsi <- linkList_new

Gsi <- graph.data.frame(edge_listsi,directed = T)
Asi <- get.adjacency(Gsi,sparse = TRUE,attr = "weight",type = "both")

g_arasi <- graph.adjacency(Asi,mode = "directed",weighted = T)
g.cyto <- igraph.to.graphNEL(g_arasi)


createNetworkFromGraph("net", graph=g.cyto)
loadTableData(
  data.frame(Scorpius_obj$modules),
  data.key.column = "feature",
  table = "node",
  table.key.column = "name")
  
analyzeNetwork(directed = TRUE)

# steps in cytoscape: 
# remove outdegree filter (1,1) AND indegree (0,0)
# remove loose ones
# color by module

#-----------------------------------------------------------------------------------------------------#
#							test with cort
#-----------------------------------------------------------------------------------------------------#

# get the DMSO data only for genes of interes
CORT_exp = logcpm[Scorpius_obj$modules$feature,Pheno$Corisol==TRUE]

weightMat_old_CORT <- GENIE3(CORT_exp,nCores=8,nTrees = 1000)

# check if edge between one pair of nodes is higher in one, or other direction

weightMat_new_CORT = weightMat_old_CORT
genelist_row = rownames(weightMat_new_CORT)
genelist_col = colnames(weightMat_new_CORT)

for(i in genelist_row ){

	for(o in genelist_col){
	
	
		if(weightMat_new_CORT[i,o] > weightMat_new_CORT[o,i]){
			weightMat_new_CORT[o,i] = 0
			}else{
			weightMat_new_CORT[i,o] = 0
		}

	}
}



# top 5% of links selected
linkList_old_CORT <- getLinkList(weightMat_old_CORT,threshold = quantile(weightMat_old_CORT,0.99))
# top 5% of links selected
linkList_new_CORT <- getLinkList(weightMat_new_CORT,threshold = quantile(weightMat_old_CORT,0.99))



head(linkList_new_CORT)
plot(linkList_new_CORT,pch=19,col=alpha("black",0.1))
#  maybe tailor the TARGETS as all genes? see what comes outof it>
# Top 5% to cytoscape

#compare
a = weightMat_new_CORT[rownames(weightMat_new),colnames(weightMat_new)] - weightMat_new

# top 5% of links selected
linkList_test <- getLinkList(a,threshold = quantile(weightMat_old_CORT,0.99))



#-----------------------------------------------------------------------------------------------------#
#							output
#-----------------------------------------------------------------------------------------------------#
save.image(file = paste0(s_OUT_dir,"GRN/GRN_img.Rdata"))
save(GRN,file = paste0(s_OUT_dir,"GRN/GRN.Rdata"))

#-----------------------------------------------------------------------------------------------------#
#							out
#-----------------------------------------------------------------------------------------------------#

# get upper 10% -> 1500 genes?
#linkList <- getLinkList(weightMat,threshold = quantile(weightMat,0.9))
#dim(linkList)