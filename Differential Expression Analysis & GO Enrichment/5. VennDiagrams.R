#################################################################
##                        Venn Diagrams                        ##
#################################################################
library("ggVennDiagram")
library("ggplot2")
library("gplots")

dataframes = c("H25_HNa","SOS1_HNa", "WT_HNa", "H23_HNa")

# We load activated and repressed genes

H25_HNa = read.table(file = "../H25_HNa vs H25_HC/repressed_genes.txt", header = F, sep = "\t")
H25_HNa = H25_HNa$V1

H23_HNa = read.table(file = "../H23_HNa vs H23_HC/repressed_genes.txt", header = F, sep = "\t")
H23_HNa = H23_HNa$V1

SOS1_HNa = read.table(file = "../SOS1_HNa vs SOS1_HC/repressed_genes.txt",header = F, sep = "\t")
SOS1_HNa = SOS1_HNa$V1


WT_HNa = read.table(file = "../WT_HNa vs WT_HC/repressed_genes.txt",header = F, sep = "\t")
WT_HNa = WT_HNa$V1



# Create a list object with gene sets
create_gene_sets = function(dataframes) {
  gene_sets = list()
  
  for (df_name in dataframes) {
    df = get(df_name)
    gene_sets[[df_name]] = df
  }
  
  return(gene_sets)
}

gene_sets = create_gene_sets(dataframes)


# Plot
figura_venn = ggVennDiagram(gene_sets, label_size = 6, set_size = 6) + 
  scale_fill_gradient(low="white",high = "blue") +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 16))  


ggsave(figura_venn,
       filename = "./HOJA_DOWN_ggVennDiagram.png",
       height = 25,width = 25,units = "cm") 


figura_upset = ggVennDiagram(gene_sets, force_upset = TRUE)
ggsave(figura_upset,
       filename = "./HOJA_DOWN_upset.png",
       height = 20,width = 30,units = "cm") 


# Access diagram elements
ItemsList = venn(gene_sets, show.plot = FALSE)
ItemsList = attributes(ItemsList)$intersections

# Ensure that all elements are of equal length.
max_length = max(lengths(ItemsList))
ItemsList = lapply(ItemsList, function(x) {
  if (length(x) < max_length) {
    c(x, rep(NA, max_length - length(x)))
  } else {
    x
  }
})

ItemsList_df = data.frame(ItemsList)

write.table(ItemsList_df, file = "HOJA_DOWN_genes.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

