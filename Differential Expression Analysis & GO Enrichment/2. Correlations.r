##################################################################
##                     Correlation analysis                     ##
##################################################################
library(corrplot)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)


# We load the normalized gene expression matrix
normalized.gene.expression = read.table(file="NORMALIZED_gene_expression_matrix.txt", header=T, row.names = 1)


# We use the cor function to create the correlations matrix
correlation = cor(normalized.gene.expression, method = "spearman")


# We will use two different types of representation styles: corrplot and pheatmap
col = colorRampPalette(c("#4477AA","#77AADD","#FFFFFF", "#EE9988","#BB4444"))

png("Heatmap_Correlation1.png", width = 30, height = 18, units = "cm", res = 500)
corrplot(correlation, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # hide correlation coefficient on the principal diagonal
         diag=T,
         is.corr = F,
         col.lim = c(min(correlation), max(correlation))
)
dev.off()


png("Heatmap_Correlation2.png", width = 30, height = 18, units = "cm", res = 500)
pheatmap(correlation,
         fontsize = 12,
         cluster_cols = T,
         cluster_rows = T, # Recommended: run the script twice: clustering and not clustering by rows.
         legend = T,
         show_rownames = T,
         col = colorRampPalette(c("blue", "white", "red"))(50))

dev.off()


# We can also create a heatmap with the annotation of our samples

Type = c(rep("WT",4),rep("H23",4),rep("H25",4),rep("SOS1",4))
Organ = rep(c("Leaf", "Root"), each = 1, 8)
Condition = rep(c("Control", "Stress"), each = 2, 4)


# Combine condition into a data frame
condition_data = data.frame(Type,Organ,Condition)
rownames(condition_data) = colnames(normalized.gene.expression)


# Create the heatmap
png("Heatmap_Correlation2.png", width = 30, height = 18, units = "cm", res = 500)

pheatmap(
  mat = as.matrix(normalized.gene.expression),
  scale = "column",
  cluster_cols = T,
  cluster_rows = T, 
  legend = T,
  col = colorRampPalette(c("blue", "white", "red"))(50),
  annotation = condition_data, 
  annotation_legend = T,
  fontsize = 12,
  border_color = F,
  labels_row = ""
)

dev.off()




# Correlation between each pair of samples: scatterplots

create_scatterplots = function(data, save_dir) {
  # Gets all column combinations
  column_combinations = combn(names(normalized.gene.expression), 2, simplify = FALSE)
  
  # Iterate over each combination and create the scatterplot.
  for (cols in column_combinations) {
    x_col = cols[1]
    y_col = cols[2]
    
    # Fit a linear model
    lm_model = lm(data[[y_col]] ~ data[[x_col]])
    r_squared = summary(lm_model)$r.squared
    
    # Create the scatterplot
    scatterplot = ggplot(data, aes_string(x = x_col, y = y_col)) +
      geom_point(color = "gray40") +  # Añade puntos para cada observación
      geom_smooth(method = "lm", se = TRUE, color = "red") +  # Añade la línea de regresión
      labs(x = x_col, y = y_col,
           title = paste("Correlation:", round(cor(data[[x_col]], data[[y_col]],method = "spearman"), 3), "| R-squared:", round(r_squared, 3))) +
      theme(plot.title = element_text(size = 12, hjust = 1),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    # Save the scatterplot as a PNG file
    ggsave(filename = paste0(save_dir, x_col, "_vs_", y_col, ".png"),
           plot = scatterplot, width = 8, height = 6)
  }
}

# Call the function with your dataframe and output directory
create_scatterplots(normalized.gene.expression, "./Plots/Scatterplots/")
