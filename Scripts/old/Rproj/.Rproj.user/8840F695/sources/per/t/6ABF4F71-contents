# Heatmap of DE genes based on condition
library(edgeR)
library(ggplot2)
library(tidyverse)
library(ggsignif)
library(ggrepel)
library(clusterProfiler)
library("org.Hs.eg.db")

# set seed
set.seed(42)

# load data
exp = read.table("C:/DATA STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/gene_count.txt",sep = "\t",header = TRUE)
pheno_gene = exp[,c("gene_id","gene_name","gene_chr","gene_start","gene_end","gene_strand","gene_length","gene_biotype","gene_description","tf_family")]

rownames(pheno_gene) = exp[,"gene_id"]
rownames(exp) = exp[,"gene_id"]

exp = exp[,!colnames(exp)%in%c("gene_id","gene_name","gene_chr","gene_start","gene_end","gene_strand","gene_length","gene_biotype","gene_description","tf_family")]
pheno = data.frame(Corisol=grepl(colnames(exp),pattern = "C"),stage=factor(gsub(colnames(exp),pattern = "C_|D_|1|2|3",replacement = ""),levels = c("pro","diffy","diffm")))


y = DGEList(counts=exp)

keep= filterByExpr(y,min.prop = 0.25, min.count=10)

y = y[keep, , keep.lib.sizes=FALSE]

y = calcNormFactors(y, method = "TMM")

logcpm = cpm(y, log=TRUE)

PCA = prcomp(t(logcpm))

library(ggplot2)
library(ggfortify)
plot(PCA)
autoplot(PCA,col=as.numeric(pheno$stage),shape=ifelse(pheno$Corisol,16,17))
autoplot(PCA,x=c(3,4),col=as.numeric(pheno$stage),shape=ifelse(pheno$Corisol,16,17))

# make barplot of relative (C/D) of gene X per stage
gene_names_interest = c("ENSG00000113580", "ENSG00000151623", "ENSG00000096060", "ENSG00000206190", "ENSG00000105880", "ENSG00000119138", "ENSG00000109906", "ENSG00000118515", "ENSG00000112679")
for(s_gene in gene_names_interest){
  if(!s_gene%in%rownames(logcpm)){next}
  
  mean_pro = mean(logcpm[s_gene,pheno$Corisol==TRUE&pheno$stage=="pro"]/logcpm[s_gene,pheno$Corisol==FALSE&pheno$stage=="pro"]
  )
  sd_pro = sd(logcpm[s_gene,pheno$Corisol==TRUE&pheno$stage=="pro"]/logcpm[s_gene,pheno$Corisol==FALSE&pheno$stage=="pro"]
  )
  
  mean_diffy = mean(logcpm[s_gene,pheno$Corisol==TRUE&pheno$stage=="diffy"]/logcpm[s_gene,pheno$Corisol==FALSE&pheno$stage=="diffy"]
  )
  sd_diffy = sd(logcpm[s_gene,pheno$Corisol==TRUE&pheno$stage=="diffy"]/logcpm[s_gene,pheno$Corisol==FALSE&pheno$stage=="diffy"]
  )
  
  mean_diffm = mean(logcpm[s_gene,pheno$Corisol==TRUE&pheno$stage=="diffm"]/logcpm[s_gene,pheno$Corisol==FALSE&pheno$stage=="diffm"]
  )
  sd_diffm = sd(logcpm[s_gene,pheno$Corisol==TRUE&pheno$stage=="diffm"]/logcpm[s_gene,pheno$Corisol==FALSE&pheno$stage=="diffm"]
  )
  
  
  # create dummy data
  data <- data.frame(
    name=c("pro","diffy","diffm"),
    value=c(mean_pro,mean_diffy,mean_diffm),
    sd=c(sd_pro,sd_diffy,sd_diffm)
  )
  
  # rectangle
  g = ggplot(data) +
    geom_bar( aes(x=name, y=value), stat="identity", fill="skyblue", alpha=0.9) +
    geom_errorbar( aes(x=name, y=value, ymin=value-sd, ymax=value+sd), width=0.5, colour="black", alpha=0.8, size=1)+
    ggtitle(s_gene)
  
  print(g)
}

# Check genes in DMSO over stage
DMSO_genes= c("ENSG00000113580", "ENSG00000151623", "ENSG00000096060", "ENSG00000112679")
for(i in 1:length(DMSO_genes)){
  boxplot(logcpm[DMSO_genes[i],pheno$Corisol==FALSE]~pheno$stage[pheno$Corisol==FALSE],main=DMSO_genes[i])
}

logcpm_original = logcpm

#################### C_PRO vs D_PRO

# Load DE set
DE_pro = read.table("C:/DATA STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/1.deglist/C_PROvsD_PRO/C_PROvsD_PRO_deg.xls",sep = "\t",header=TRUE)

DE_pro = DE_pro[filterByExpr(DE_pro[,22:27],min.count = 5),]
plot(DE_pro$log2FoldChange,-log10(DE_pro$pvalue))
DE_genes_pro = DE_pro$gene_id[which(DE_pro$pvalue<0.001 & abs(DE_pro$log2FoldChange)>log2(2))]


DE_genes_pro = DE_genes_pro[DE_genes_pro%in%rownames(logcpm_original)]

logcpm_PRO = logcpm_original[DE_genes_pro,]



rownames(logcpm_PRO ) = DE_pro$gene_name[match(rownames(logcpm_PRO),DE_pro$gene_id)]


library("gplots")

# for pro
pdf(file = paste0("Heatmap_pro.pdf"))
heatmap.2(logcpm_PRO[,pheno$stage=="pro"],
          col = colorpanel(100,"red","black","green"), 
          symkey=F,symm=F,symbreaks=T, scale="none",
          trace = "none", density.info = "none")
dev.off()

pdf(file = paste0("Heatmap_pro_differentscaling.pdf"))
heatmap.2(logcpm_PRO[,pheno$stage=="pro"],
          col = colorpanel(100,"black","yellow"), 
          symkey=F,symm=F,symbreaks=FALSE, scale="none",
          trace = "none", density.info = "none")
dev.off()

#################### C_Diffy vs D_Diffy

# Load DE set
DE_diffy = read.table("C:/DATA STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/1.deglist/C_DiffyvsD_Diffy/C_DiffyvsD_Diffy_deg.xls",sep = "\t",header=TRUE)

DE_genes_diffy = DE_diffy$gene_id[which(DE_diffy$pvalue<0.001 & abs(DE_diffy$log2FoldChange)>log2(2))]


DE_genes_diffy = DE_genes_diffy[DE_genes_diffy%in%rownames(logcpm_original)]

logcpm_Diffy = logcpm_original[DE_genes_diffy,]



rownames(logcpm_Diffy ) = DE_diffy$gene_name[match(rownames(logcpm_Diffy),DE_diffy$gene_id)]


# for diffy
pdf(file = paste0("Heatmap_diffy.pdf"))
heatmap.2(logcpm_Diffy[,pheno$stage=="diffy"], col = colorpanel(100,"red","black","green"), 
          symkey=F,symm=F,symbreaks=T, scale="none",
          trace = "none", density.info = "none")
dev.off()

pdf(file = paste0("Heatmap_diffy_differentscaling.pdf"))
heatmap.2(logcpm_Diffy[,pheno$stage=="diffy"], col = colorpanel(100,"black","yellow"), 
          symkey=F,symm=F,symbreaks=FALSE, scale="none",
          trace = "none", density.info = "none")
dev.off()

#################### C_Diffm vs D_Diffm

# Load DE set
DE_diffm = read.table("C:/DATA STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/1.deglist/C_DiffmvsD_Diffm/C_DiffmvsD_Diffm_deg.xls",sep = "\t",header=TRUE)

DE_genes_diffm = DE_diffm$gene_id[which(DE_diffm$pvalue<0.001 & abs(DE_diffm$log2FoldChange)>log2(2))]


DE_genes_diffm = DE_genes_diffm[DE_genes_diffm%in%rownames(logcpm_original)]

logcpm_Diffm = logcpm_original[DE_genes_diffm,]



rownames(logcpm_Diffm ) = DE_diffm$gene_name[match(rownames(logcpm_Diffm),DE_diffm$gene_id)]


# for diffm
pdf(file = paste0("Heatmap_diffm.pdf"))
heatmap.2(logcpm_Diffm[,pheno$stage=="diffm"], col = colorpanel(100,"red","black","green"), 
          symkey=F,symm=F,symbreaks=T, scale="none",
          trace = "none", density.info = "none")
dev.off()

pdf(file = paste0("Heatmap_diffm_differentscaling.pdf"))
heatmap.2(logcpm_Diffm[,pheno$stage=="diffm"], col = colorpanel(100,"black","yellow"), 
          symkey=F,symm=F,symbreaks=F, scale="none",
          trace = "none", density.info = "none")
dev.off()


# overlap

# Trajectory using CoGAPS
# devtools::install_github("FertigLab/CoGAPS")
library("CoGAPS")

minrange= abs(range(logcpm)[1]) + 1
logcpm_edit = (logcpm+minrange)
#a= CoGAPS(logcpm+minrange, nIterations=1000)

# set the value for a specific parameter
D_GAPS = CoGAPS(logcpm_edit[,pheno$Corisol==FALSE], nIterations=10000)
C_GAPS = CoGAPS(logcpm_edit[,pheno$Corisol==TRUE], nIterations=10000)
plot(D_GAPS)
plot(C_GAPS)
# plot(C_GAPS)
save(D_GAPS,file="D_GAPS.Rdata")
save(C_GAPS,file="C_GAPS.Rdata")
# result <- CoGAPS(GIST.matrix, params, nIterations=5000, messages=FALSE)

pdf(file = paste0("CoGAPS_DMSO_Time_patterns.pdf"))
  plot(D_GAPS)
dev.off()
  
# pattern 2 is low,high,high
# get top 100 genes with pattern
patterncluster = data.frame(GeneID = rownames(logcpm_edit), PatternID = NA)
s_th_pattern = 0.05
for(i in 1:7){
  csum = cumsum(sort((D_GAPS@featureLoadings[,i]),decreasing = T))
  csum = (csum-min(csum))/(max(csum) - min(csum))
  
  # make contribution cumsum
  plot(csum,main=i)
  abline(h=s_th_pattern)
  
  
  Top_2_genes = names(which(csum<=s_th_pattern))#names(sort(-D_GAPS@featureLoadings[,i])[1:20])
  
  #make boxplot
  pdf(file = paste0("Boxplot_genes_pattern_",i,".pdf"))
    boxplot(logcpm_edit[Top_2_genes,],main=paste0("Pattern ",i,"\nTop ",s_th_pattern*100,"% contibuting genes"),las=2, col="white",
            border=ifelse(pheno$Corisol,yes = "red3",no = "blue3"))
  dev.off()
  
  # make cluster
  patterncluster[patterncluster$GeneID%in%Top_2_genes,"PatternID"] = i
}

for(i in 1:7){
  pdf(file = paste0("CoGAPS_DMSO_Time_pattern_",i,".pdf"))
    plot(D_GAPS@sampleFactors[,i],type = "l",lwd=2,ylim=c(0,1),
         main=paste0("Pattern ",i), ylab="Relative Amplitude", xlab="Samples")
    points(D_GAPS@sampleFactors[,i],lwd=2,ylim=c(0,1),pch=1)
  dev.off()
}


table(patterncluster$PatternID)

# check patten genes for DEGs of contrast.
pro_genes = DE_pro[DE_pro$padj<0.05,"gene_id"]
diffy_genes = DE_diffy[DE_diffy$padj<0.05,"gene_id"]
diffm_genes = DE_diffm[DE_diffm$padj<0.05,"gene_id"]

total_genes = c(pro_genes,diffy_genes,diffm_genes)

sign_patterncluster = patterncluster[patterncluster$GeneID %in% unique(total_genes),]
sign_patterncluster = sign_patterncluster[!is.na(sign_patterncluster$PatternID),]

#@RRR stats on ocurrens signifficance compared to random
# add significance to pattern
sign_patterncluster$contr_pattern = NA
sign_patterncluster$pval_pro = NA
sign_patterncluster$pval_diffy = NA
sign_patterncluster$pval_diffm = NA

for(i in 1:7){
    sign_patterncluster[sign_patterncluster$PatternID==i,"contr_pattern"]  = 
      as.numeric(D_GAPS@featureLoadings[match( sign_patterncluster[sign_patterncluster$PatternID==i,"GeneID"],rownames( D_GAPS@featureLoadings)),i]
  )
}

# check sign for stages
matchid = match( sign_patterncluster$GeneID,DE_pro$gene_id)
sign_patterncluster[!is.na(matchid),"pval_pro"]  =  DE_pro$padj[matchid][!is.na(matchid)]

matchid = match( sign_patterncluster$GeneID,DE_diffy$gene_id)
sign_patterncluster[!is.na(matchid),"pval_diffy"]  =  DE_diffy$padj[matchid][!is.na(matchid)]

matchid = match( sign_patterncluster$GeneID,DE_diffm$gene_id)
sign_patterncluster[!is.na(matchid),"pval_diffm"]  =  DE_diffm$padj[matchid][!is.na(matchid)]

# make function to plot
f_plottraj = function(s_pattern=1,s_genes=1){
  # Look at pattern 3, expected diff makes biological effect
  for( i in s_pattern){
    sign3_patterncluster_genes = sign_patterncluster$GeneID[sign_patterncluster$PatternID==i]
    
    for(o in s_genes){
      gene_index = rownames(logcpm)%in%sign3_patterncluster_genes[o]
      
      gene_name = rownames(logcpm)[gene_index]
      
      gene_name_symbol = c(DE_pro$gene_name,DE_diffy$gene_name,DE_diffm$gene_name)[match(gene_name,c(DE_pro$gene_id,DE_diffy$gene_id,DE_diffm$gene_id))]
      
      
      
      plotdata = data.frame(D_pro = as.numeric(logcpm[gene_index,pheno$Corisol==FALSE][1:3]),
                            D_diffy = as.numeric(logcpm[gene_index,pheno$Corisol==FALSE][4:6]),
                            D_diffm = as.numeric(logcpm[gene_index,pheno$Corisol==FALSE][7:9]),
                            C_pro = as.numeric(logcpm[gene_index,pheno$Corisol==TRUE][1:3]),
                            C_diffy = as.numeric(logcpm[gene_index,pheno$Corisol==TRUE][4:6]),
                            C_diffm = as.numeric(logcpm[gene_index,pheno$Corisol==TRUE][7:9]))
      
      plotdata2 = pivot_longer(plotdata,1:6)
      plotdata2$stage = factor(gsub(plotdata2$name,pattern = "._",replacement = ""),levels = c("pro","diffy","diffm"))
      plotdata2$cond = factor(gsub(plotdata2$name,pattern = "_.*",replacement = ""),levels = c("D","C"))
      # gg boxplot
      p <- ggplot(plotdata2, aes(x=stage, y=value, col = cond)) + 
        geom_boxplot() + 
        ggtitle(paste0("Pattern ", i, "\nGene ",gene_name," (",gene_name_symbol,")"))#+
        #geom_signif(comparisons = list(c("stage")), #c("D_pro", "C_pro"),c("D_diffy", "C_diffy"),c("D_diffm", "C_diffm")
                    #map_signif_level=TRUE)
        #+ geom_line(aes(x=stage,y=value,group=cond));p
      print(p)
      ggsave(plot = p, dpi = 300,filename = paste0("Trajectory_Pattern_", i, "_Gene_",gene_name_symbol,".jpeg"))
    }
  }
}


# Change order of object
sign_patterncluster_backup = sign_patterncluster
sign_patterncluster = sign_patterncluster[order(sign_patterncluster$contr_pattern,decreasing = T),]
# plot
f_plottraj(s_pattern=4,s_genes=1)
f_plottraj(s_pattern=4,s_genes=3)
f_plottraj(s_pattern=1:7,s_genes=1)

#################################################################
#             ONTOLOGIES                                        #
#################################################################


# next step ontology of top 5% for all patterns
# using clusterprofiler

Ont_results =list()
for( s_pattern in 1:7){
  pattern = s_pattern
  sample_gene = sign_patterncluster_backup$GeneID[sign_patterncluster_backup$PatternID==pattern]
  saple_gene_importance = sign_patterncluster_backup$contr_pattern[sign_patterncluster_backup$PatternID==pattern]
  universe = c(DE_pro$gene_name,DE_diffy$gene_name,DE_diffm$gene_name)
  sample_gene = universe[match(sample_gene,c(DE_pro$gene_id,DE_diffy$gene_id,DE_diffm$gene_id))]
  names(saple_gene_importance) = sample_gene
  
  geneList = sort(saple_gene_importance,decreasing = T)
  gene  = names(geneList)# names(geneList)[geneList>as.numeric(quantile(saple_gene_importance,0.5))]
    
  ggo3 <- groupGO(gene     = gene,
                  OrgDb    = org.Hs.eg.db,
                  ont      = "BP",
                  level    = 7,
                  keyType = "SYMBOL")
  head(ggo3[order(ggo3@result$Count,decreasing = T)])
  
  Ont_results[[s_pattern]] = ggo3[order(ggo3@result$Count,decreasing = T)]
}

lapply(Ont_results,function(x){x[1:10,2]})

# 5 is good level, might go deeper
# 7 is nice
Ont_1 = Ont_results[[1]]
Ont_2 = Ont_results[[2]]
Ont_3 = Ont_results[[3]]
Ont_4 = Ont_results[[4]]
Ont_5 = Ont_results[[5]]
Ont_6 = Ont_results[[6]]
Ont_7 = Ont_results[[7]]

write.table(Ont_1 ,file = "Ont_1_tab.txt",quote = FALSE,sep = "\t",row.names = FALSE)
write.table(Ont_2 ,file = "Ont_2_tab.txt",quote = FALSE,sep = "\t",row.names = FALSE)
write.table(Ont_3 ,file = "Ont_3_tab.txt",quote = FALSE,sep = "\t",row.names = FALSE)
write.table(Ont_4 ,file = "Ont_4_tab.txt",quote = FALSE,sep = "\t",row.names = FALSE)
write.table(Ont_5 ,file = "Ont_5_tab.txt",quote = FALSE,sep = "\t",row.names = FALSE)
write.table(Ont_6 ,file = "Ont_6_tab.txt",quote = FALSE,sep = "\t",row.names = FALSE)
write.table(Ont_7 ,file = "Ont_7_tab.txt",quote = FALSE,sep = "\t",row.names = FALSE)



#################################################################
#             volcanoplot                                        #
#################################################################
# redo the DE, low counts not filtered.. -> is OK
# do heatmap
# include pvalue col

# volcano of DE
de = DE_diffy
 
s_pval= 0.0005
s_lfc_th = 0.5

f_volcano = function(de,s_title="",s_pval = 0.05, s_lfc_th = 0.5){
  #de = de[filterByExpr(de[,2:7],min.count = 5),]
  de$gene_name[!(abs(de$log2FoldChange) > s_lfc_th & de$pvalue < s_pval)] = NA
  colors = ifelse(de$log2FoldChange>s_lfc_th,yes = "Sign. upregulated",
                  no = ifelse(de$log2FoldChange< -s_lfc_th,
                              yes = "Sign. downregulated",no = ""))
  
  colors[de$pvalue > s_pval] = ""
  
  ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), 
                      col=colors, label=gene_name)) +
  geom_point() + ggtitle(s_title)+
  #theme_minimal() +
  geom_text_repel(max.overlaps=10) +
  scale_color_manual(values=c("black", "red3", "green3"))+
  theme(legend.position="none")
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  #geom_hline(yintercept=-log10(0.05), col="red")
}

f_volcano(DE_pro,s_title="PRO")
f_volcano(DE_diffy,s_title="Diffy")
f_volcano(DE_diffm,s_title="Diffm")
# GO of DE

# PPI of DE


# save logCPM
write.table(logcpm,file = "C:/DATA STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Rproj/logcpm.txt",sep = "\t",quote = FALSE)
















