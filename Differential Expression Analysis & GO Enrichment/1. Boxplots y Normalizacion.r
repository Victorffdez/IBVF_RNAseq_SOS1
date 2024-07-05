#################################################################
##          Boxplot representation of gene expression          ##
#################################################################

# Load the raw data matrix
gene.expression = read.table(file = "RAW_gene_expression_matrix.txt", header = T, sep = "\t", dec = ".", row.names = 1)
gene.expression = gene.expression[apply(gene.expression[,-1], 1, function(x) !all(x==0)),]


# Add 1 to all expression levels in order to normalize. The problem is caused by x < 1 --> log2(x) < 0
gene.expression.1 <- gene.expression + 1

# Save the raw gene expression data +1
write.table(x = gene.expression.1,file = "RAW_gene_expression_matrix_+1.txt",
            quote = F,row.names = T, col.names = T,
            sep = "\t")

# Representation of the global distribution of gene expression
boxplot(gene.expression, outline = FALSE, col = rep(c("seagreen", "orange3","tomato", "steelblue3"), each = 4),
        ylab = "Gene Expression (FPKM)",
        xlab = "Samples",
        cex.lab = 1.25,     # Aumenta el tamaño de la etiqueta de los ejes
        cex.axis = 0.8)    # Aumenta el tamaño de los labels del eje X)



# Let's see if the logarithm in base 2 improves the distribution
log2normalized.gene.expression <- read.table(file="./Normalyzer Results/log2-normalized.txt", header=T, row.names = 1)
boxplot(log2normalized.gene.expression, outline = FALSE, col = rep(c("seagreen", "orange3","tomato", "steelblue3"), each = 4),
        ylab = "Normalized Gene Expression (log2(FPKM+1))",
        xlab = "Samples",
        cex.lab = 1.25,     # Aumenta el tamaño de la etiqueta de los ejes
        cex.axis = 0.8)    # Aumenta el tamaño de los labels del eje X)



##################################################################
##                      Data normalization                      ##
##################################################################

# We are going to use the NormalyzerDE package to normalize
BiocManager::install("NormalyzerDE", force = T)
library(NormalyzerDE)


# To apply standardization we have to indicate which samples are replicates of which condition. 
# It is essential to redefine the experimental design. It has two columns: sample and group. 

design <- data.frame(sample=colnames(gene.expression),
                     group=colnames(gene.expression))

write.table(x = design,file = "./RESULTADOS/normalyzer_design.tsv",quote = F,row.names = F, sep = "\t")

# We use the normalyzer function to normalize
normalyzer(jobName = "NORMALIZACION", designPath = "./DESIGN.txt",
           dataPath = "./RAW_gene_expression_matrix_+1.txt",outputDir = "./",
           requireReplicates = FALSE)


# After observing the results of normalyzer, let's represent some of the results for comparison
log2normalized.gene.expression <- read.table(file="./Normalyzer Results/log2-normalized.txt", header=T, row.names = 1)
boxplot(log2normalized.gene.expression, outline = FALSE, col = rep(c("seagreen", "orange3","tomato", "steelblue3"), each = 4),
        ylab = "Normalized Gene Expression (log2(FPKM+1))",
        xlab = "Samples",
        cex.lab = 1.25,     # Aumenta el tamaño de la etiqueta de los ejes
        cex.axis = 0.8)    # Aumenta el tamaño de los labels del eje X)



vsnnormalized.gene.expression <- read.table(file="./Normalyzer Results/VSN-normalized.txt", header=T, row.names = 1)
boxplot(vsnnormalized.gene.expression, outline = FALSE, col = rep(c("seagreen", "orange3","tomato", "steelblue3"), each = 4),
        ylab = "Normalized Gene Expression (log(FPKM+1))",
        xlab = "Samples",
        cex.lab = 1.25,     # Aumenta el tamaño de la etiqueta de los ejes
        cex.axis = 0.8)    # Aumenta el tamaño de los labels del eje X)



cuantilenormalized.gene.expression <- read.table(file="./Normalyzer Results/Quantile-normalized.txt", header=T, row.names = 1)
boxplot(cuantilenormalized.gene.expression, outline = FALSE, col = rep(c("seagreen", "orange3","tomato", "steelblue3"), each = 4),
        ylab = "Normalized Gene Expression (log(FPKM+1))",
        xlab = "Samples",
        cex.lab = 1.25,     # Aumenta el tamaño de la etiqueta de los ejes
        cex.axis = 0.8)    # Aumenta el tamaño de los labels del eje X)


mediannormalized.gene.expression <- read.table(file="./Normalyzer Results/median-normalized.txt", header=T, row.names = 1)
boxplot(mediannormalized.gene.expression, outline = FALSE, col = rep(c("seagreen", "orange3","tomato", "steelblue3"), each = 4),
        ylab = "Normalized Gene Expression (log(FPKM+1))",
        xlab = "Samples",
        cex.lab = 1.25,     # Aumenta el tamaño de la etiqueta de los ejes
        cex.axis = 0.8)    # Aumenta el tamaño de los labels del eje X)
