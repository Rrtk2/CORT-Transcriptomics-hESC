
#-----------------------------------------------------------------------------------------------------#
# 							GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#
# File description:
#	Name
#		Volcanoplot.R
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
#
#-----------------------------------------------------------------------------------------------------#
#							Input
#-----------------------------------------------------------------------------------------------------#
load("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Data_QC/DE/DE_res_all.Rdata")


# needs SvR DE data
#load("~/GitHub/PRISMO_Transcriptomics/Scripts/POST only processing/data_DE.RData")

# Load packages
library(tidyverse)
library(ggrepel)
library(clusterProfiler)
library("org.Hs.eg.db")

s_pvalTH = 0.05
s_fcTH = 1.5
# data
data  = DE_res_all$pro_diffy_cortVSpro_diffy_dmso
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
showtopgenesnumber = 5
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
    ggtitle("Volcanoplot, pro_diffy_cortVSpro_diffy_dmso") # - logFC cutoff
p5 


# titles and such


# save
#Rplot(p5,"DE_Volcano",width = 480*2)

	
	
# Full Cleanup
#rm(list=ls())

