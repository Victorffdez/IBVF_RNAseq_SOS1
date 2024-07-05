#################################################################
##           Differential Expression Analysis: limma           ##
#################################################################
library(limma)
library(dplyr)

# We load the normalized gene expression matrix
normalized.gene.expression = read.table(file="NORMALIZED_gene_expression_matrix.txt", header=T, row.names = 1)

# Specify the experimental design
experimental.design = model.matrix(~ -1+factor(c(1:16)))
samples = colnames(normalized.gene.expression)
colnames(experimental.design) = samples


# Fit the estimation of the expression levels of each gene to a linear model (means)
linear.fit = lmFit(normalized.gene.expression, experimental.design)


# Create as many comparisons as possible
combs = c()
for (i in 1:(ncol(normalized.gene.expression) - 1)) {
  for (j in (i + 1):ncol(normalized.gene.expression)) {
    combs = c(combs, paste0(colnames(normalized.gene.expression)[i], ".vs.", colnames(normalized.gene.expression)[j]))
    combs = c(combs, paste0(colnames(normalized.gene.expression)[j], ".vs.", colnames(normalized.gene.expression)[i]))
  }
}

# Create the contrast matrix
contrast.matrix = matrix(0, ncol = length(combs), nrow = ncol(experimental.design))
colnames(contrast.matrix) = combs
rownames(contrast.matrix) = colnames(experimental.design)

# Assign 1 and -1 according to the comparison
for (i in 1:length(combs)) {
  col1 = strsplit(combs[i], ".vs.")[[1]][1]
  col2 = strsplit(combs[i], ".vs.")[[1]][2]
  
  row1 = which(rownames(contrast.matrix) == col1)
  row2 = which(rownames(contrast.matrix) == col2)
  
  contrast.matrix[row1, i] = 1
  contrast.matrix[row2, i] = -1
}


# We calculate the log(Fold-Change)
contrast.linear.fit = contrasts.fit(linear.fit, contrast.matrix)
Resultado_logFC = as.data.frame(contrast.linear.fit[["coefficients"]])
genes = rownames(Resultado_logFC)
write.table(x = Resultado_logFC, file = "LogFC_Matrix.txt", quote = F, sep = "\t", row.names = T, col.names = T)


# We obtain the up and down regulated genes. Comparison of interest: H25_RNa.vs.WT_RNa
H25_RNa_vs_WT_RNa = as.data.frame(Resultado_logFC$`H25_RNa.vs.WT_RNa`)
colnames(H25_RNa_vs_WT_RNa) = "LogFC"
rownames(H25_RNa_vs_WT_RNa) = genes


# Activated genes
activated_genes = rownames(H25_RNa_vs_WT_RNa[H25_RNa_vs_WT_RNa$LogFC > 2, , drop = FALSE])
length(activated_genes)

activated_df = as.data.frame(H25_RNa_vs_WT_RNa[activated_genes, ])
colnames(activated_df) = "LogFC"
rownames(activated_df) = activated_genes


#Order the dataframe from largest to smallest# 
new_row_names = rownames(activated_df)[order(activated_df$LogFC, decreasing = TRUE)] 
activated_df = as.data.frame(activated_df[order(activated_df$LogFC,decreasing=TRUE),])
rownames(activated_df) = new_row_names
colnames(activated_df) = "LogFC"
activated_genes = row.names(activated_df)


# Repressed genes
repressed_genes = rownames(H25_RNa_vs_WT_RNa[H25_RNa_vs_WT_RNa$LogFC < -2, , drop = FALSE])
length(repressed_genes)
repressed_df = as.data.frame(H25_RNa_vs_WT_RNa[repressed_genes, ])
colnames(repressed_df) = "LogFC"
rownames(repressed_df) = repressed_genes


# Order the dataframe from largest to smallest
new_row_names = rownames(repressed_df)[order(repressed_df$LogFC, decreasing = TRUE)] # Importante no cambiar el orden de esta fila
repressed_df = as.data.frame(repressed_df[order(repressed_df$LogFC,decreasing=TRUE),])
rownames(repressed_df) = new_row_names
colnames(repressed_df) = "LogFC"
repressed_genes = row.names(repressed_df)


# Representing these results
all_regulated_genes = rbind(activated_df, repressed_df)

colors = ifelse(all_regulated_genes$LogFC > 0, "indianred", "lightskyblue")
y_range = range(all_regulated_genes$LogFC)

png("./Resultados de las comparaciones/H25_RNa vs WT_RNa/barplot.png", width = 40, height = 20, units = "cm", res = 1000)
barplot(all_regulated_genes$LogFC, 
        names.arg = rownames(all_regulated_genes),
        ylab = "log2(Fold-Change)",
        ylim = y_range,  
        col = colors,
        main = "H25_RNa_vs_WT_RNa",
        border = "gray1",   
        las = 3,        
        cex.names = 0.57,
        font = 2,
        space = 0,
        yaxt = "n")
axis(side = 2) 
dev.off()



# Save the genes and dataframes 
write.table(x = activated_genes, file = "./5. Annotation/input/activated_genes.txt", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(x = repressed_genes, file = "./5. Annotation/input/repressed_genes.txt", quote = F, row.names = F, col.names = F, sep = "\t")

write.table(x = activated_genes, file = "./Resultados de las comparaciones/H25_RNa vs WT_RNa/activated_genes.txt", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(x = repressed_genes, file = "./Resultados de las comparaciones/H25_RNa vs WT_RNa/repressed_genes.txt", quote = F, row.names = F, col.names = F, sep = "\t")


write.table(x = activated_df, file = "./5. Annotation/input/activated_df.txt", quote = F, row.names = T, col.names = T, sep = "\t")
write.table(x = repressed_df, file = "./5. Annotation/input/repressed_df.txt", quote = F, row.names = T, col.names = T, sep = "\t")

write.table(x = activated_df, file = "./Resultados de las comparaciones/H25_RNa vs WT_RNa/activated_df.txt", quote = F, row.names = T, col.names = T, sep = "\t")
write.table(x = repressed_df, file = "./Resultados de las comparaciones/H25_RNa vs WT_RNa/repressed_df.txt", quote = F, row.names = T, col.names = T, sep = "\t")
