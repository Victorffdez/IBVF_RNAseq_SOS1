#################################################################
##         Exploratory data analysis: PCA & Clustering         ##
#################################################################
library(FactoMineR)
library(factoextra)
library(dendextend)
library(ggplot2)

# We load the normalized gene expression matrix
normalized.gene.expression = read.table(file="./NORMALIZED_gene_expression_matrix.txt", header=T, row.names = 1)

# FactoMineR expects the data to be transposed (genes in columns and conditions in rows)
pca.gene.expression = data.frame(colnames(normalized.gene.expression), t(normalized.gene.expression))
colnames(pca.gene.expression)[1] = "Sample"

# Apply PCA function
res.pca = PCA(pca.gene.expression, graph = F,scale.unit = T,quali.sup = 1)
res.hcpc = HCPC(res.pca, graph=FALSE, nb.clust = 5)  

# We perform clustering by testing various values of k
fviz_dend(res.hcpc,k=5,
          cex = 0.75,                    # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
          labels_track_height = 1400     # Augment the room for labels
)


# We repeat the hierarchical clustering using the R Base hclust function

dist_matrix = dist(t(normalized.gene.expression)) 
hclust_result = hclust(dist_matrix)  

# Plot the resulting dendrogram
dend = as.dendrogram(hclust_result)
colored_dend = color_branches(dend, k = 5)  # Colorear en 5 grupos

# Customize the display
fviz_dend(colored_dend, k = 5,
          main = "Dendrograma de Clustering Jerárquico",
          type = "rectangle",
          k_colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),  # Colores para los grupos
          k_labels = TRUE,  # Mostrar etiquetas de grupos
          k_label_colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),  # Colores de etiquetas
          cex = 0.7,  # Tamaño de las etiquetas
          rect = TRUE, rect_fill = TRUE, rect_border = "white", # Agregar rectángulos
          linewidth = 4
)


# Principal Component Analysis (PCA)
samples = colnames(normalized.gene.expression)

# A easy form of representation
fviz_pca_ind(res.pca, col.ind = samples, 
             pointsize=2, pointshape=21,fill.ind = samples,
             repel = TRUE, 
             addEllipses = F,ellipse.type = "confidence",
             legend.title="Conditions",
             title="Principal Component Analysis",
             show_legend=TRUE,show_guide=TRUE)


# Extract the results of the PCA
pcaData = data.frame(res.pca$ind$coord)
pcaData$Condition = pca.gene.expression$Sample

# Calculate the variance explained by each principal component
percentVar = round(res.pca$eig[1:2, 2], 2)



# Define the colors and shapes for each condition (we have 16 conditions)
shapes = rep(16, 16)
colors = rep(c("#00AFBB", "#FC4E07", "#A6CEE3", "#1F78B4", 
                "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", 
                "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", 
                "#FFFF99", "#B15928", "#FF33A5", "#33FFA5"), each = 1)


# Create the PCA chart
ggplot(pcaData, aes(Dim.1, Dim.2, color=Condition, shape=Condition, group =Condition)) +
  ggtitle("PCA PLOT") +
  geom_point(size=10) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) + 
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  theme_minimal() +  
  theme(legend.position = "right",
        axis.text = element_text(size = 11, face = "bold", color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11, face = "bold", color = "black"),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"))
