#-----------------------------------------------------------------------------------------------------#
#							Input
#-----------------------------------------------------------------------------------------------------#
temp_dat = Data # load(C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Data_QC/Salmon/Data.Rdata")
temp_pheno = Pheno # load(C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Data_QC/Salmon/Pheno.Rdata")

#-----------------------------------------------------------------------------------------------------#
#							main
#-----------------------------------------------------------------------------------------------------#
# code source: https://www.bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/conditionsVignette.html

# Reduce dim by pca
# first select genes of interest!
temp_pcadata = logcpm[DE_res2$Bestgenes,]

# get the correct PCs; select the ones that are relevant!
PC_stage = which.max(abs(cor(pca$x,as.numeric(temp_pheno$stage))));PC_stage;PC_stage=2#; PC_stage = 1
PC_Cort = which.max(abs(cor(pca$x,as.numeric(temp_pheno$Corisol))));PC_Cort #;PC_Cort = 2

pca = prcomp(t(temp_pcadata),scale = FALSE) #@RRR logcpm now
rd = pca$x[,c(PC_stage,PC_Cort)]
condition = factor(temp_pheno$Corisol)
cl = as.numeric(temp_pheno$stage)
# data
#data('slingshotExample')
#rd <- slingshotExample$rd # so thi is reduced dim; PC 1 and 2
#cl <- slingshotExample$cl # i think is cell stage
#condition <- factor(rep(c('A','B'), length.out = nrow(rd))) # condition is treatment
#condition[110:140] <- 'A'

# PCA like plot
	plot(rd, asp = 1, pch = 16, col = brewer.pal(3,'Set1')[condition], las=1)
	legend('topleft','(x,y)',legend = c('A','B'), title = 'Condition', pch=16, col = brewer.pal(3,'Set1')[1:2])

# Trajectory inference
	pto <- slingshot(rd, cl)
	
#Trajectory plot
	plot(rd, asp = 1, pch = 16, col = brewer.pal(3,'Set1')[condition], las=1)
	lines(SlingshotDataSet(pto), lwd=3)
	legend('topleft','(x,y)',legend = c('TRUE','FALSE'), title = 'Condition', pch=16, col = brewer.pal(3,'Set1')[1:2])
	
#lineage plot
	n <- nrow(rd); L <- ncol(slingPseudotime(pto))
	noise <- runif(n, -.1,.1)
	plot(as.numeric(slingPseudotime(pto)), rep(1:L, each=n)+noise,pch=16, col = brewer.pal(9,'Set1')[condition], axes=FALSE, xlab='Pseudotime', ylab='Lineage', las=1)
	axis(1); axis(2, at=1:L, las=1)


# lineage distir diff boxplot
	boxplot(slingPseudotime(pto)[,1] ~ condition, col = brewer.pal(3,'Set1')[1:2], main = 'Lineage 1', xlab='Condition', ylab='Pseudotime', las=1, pch = 16)


# Permutation test
	t1 <- slingPseudotime(pto, na=FALSE)[,1]
	w1 <- slingCurveWeights(pto)[,1]
	d1 <- weighted.mean(t1[condition=='TRUE'], w1[condition=='TRUE']) - 
		weighted.mean(t1[condition=='FALSE'], w1[condition=='FALSE'])
	dist1 <- replicate(1e4, {
		condition.i <- sample(condition)
		weighted.mean(t1[condition.i=='TRUE'], w1[condition.i=='TRUE']) - 
			weighted.mean(t1[condition.i=='FALSE'], w1[condition.i=='FALSE'])
	})


	paste0('Lineage 1 p-value: ', mean(abs(dist1) > abs(d1)))
	# if sign means that is difference between conditions; it is not


# Kolmogorov-Smirnov Test
	ks.test(slingPseudotime(pto)[condition=='TRUE',1], slingPseudotime(pto)[condition=='FALSE',1])
	# if sign means that is difference between conditions; it is not

if(FALSE){

	# More downstream analysis
	# source: https://www.bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html

	# fit negative binomial GAM
	# @RRR here things go wrong
	sce <- slingshot(rd, clusterLabels = 'GMM', reducedDim = 'PCA')

	sce <- fitGAM(sce)

	# test for dynamic expression
	ATres <- associationTest(pto)

	topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
	pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
	heatdata <- assays(sce)$counts[topgenes, pst.ord]
	heatclus <- sce$GMM[pst.ord]

	heatmap(log1p(heatdata), Colv = NA,
			ColSideColors = brewer.pal(9,"Set1")[heatclus])
			
	#-----------------------------------------------------------------------------------------------------#
	#							output
	#-----------------------------------------------------------------------------------------------------#
}
save(logcpm,file = "C:/DATA_STORAGE/Projects/Katherine_Trajectory_analysis/Data_QC/Slingshot/Slingshot.Rdata")