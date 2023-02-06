#-----------------------------------------------------------------------------------------------------#
#							Call settings
#-----------------------------------------------------------------------------------------------------#
source("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Scripts/Settings.R")

#-----------------------------------------------------------------------------------------------------#
#							PREP DATA
#-----------------------------------------------------------------------------------------------------#
# Get data and pheno
Data_raw = read.table("C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Data_RAW/gene_count.txt",sep = "\t",header = TRUE)

temp_dat = Data_raw[,!colnames(Data_raw)%in%c("gene_id","gene_name","gene_chr","gene_start","gene_end","gene_strand","gene_length","gene_biotype","gene_description","tf_family")]
temp_pheno = data.frame(Corisol=grepl(colnames(temp_dat),pattern = "C"),stage=factor(gsub(colnames(temp_dat),pattern = "C_|D_|1|2|3",replacement = ""),levels = c("pro","diffy","diffm")))

# Make final pheno and data
Data = temp_dat
Pheno = temp_pheno

# get names
rownames(Pheno) = colnames(temp_dat)
rownames(Data) = Data_raw[,"gene_id"]

#-----------------------------------------------------------------------------------------------------#
#							main
#-----------------------------------------------------------------------------------------------------#
y = DGEList(counts=Data)

keep= filterByExpr(y)

y = y[keep, , keep.lib.sizes=FALSE]

y = calcNormFactors(y, method = "TMM")

logcpm = edgeR::cpm(y, log=TRUE)


#-----------------------------------------------------------------------------------------------------#
#							Save 
#-----------------------------------------------------------------------------------------------------#
save(logcpm,file = paste0(s_OUT_dir,"Prep/logcpm.Rdata"))
save(Pheno,file = paste0(s_OUT_dir,"Prep/Pheno.Rdata"))
save(Data,file = paste0(s_OUT_dir,"Prep/Data.Rdata"))

