temp_requiredpackages=c("limma","edgeR","ggplot2","tidyverse","slingshot","SCORPIUS","RColorBrewer","tradeSeq")
#"slingshot","RColorBrewer","tradeSeq","statOmics",

# Install packages if needed and load, or just load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", ask = F)

for (i in temp_requiredpackages) {
	if (!requireNamespace(i, quietly = TRUE))
		BiocManager::install(i, ask = F)  # dependencies = c("Depends", "Imports")
	#print(i)
}

for (i in temp_requiredpackages) {
	require(as.character(i), character.only = TRUE)
	#print(i)
}
