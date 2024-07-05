#################################################################
##          Load and prepare environment for Enricher          ##
#################################################################
library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(AnnotationHub) 
library(data.table)

# Read genes list
library(readxl)
genes = read_xlsx(path = "../../Resultados de las comparaciones/5. COMPARACIONES (H23-SOS1-WT vs sus controles)/RAIZ_DOWN_genes.xlsx")
genes = genes$H23_RNa
genes = genes[!is.na(genes)]


# Read genes list
#genes = read.table("../input/activated_genes.txt", header=F)
#genes = genes[[1]]

# Read genes R
#genes = read.table("../../Resultados de las comparaciones/5. COMPARACIONES (H23-SOS1-WT vs sus controles)/HOJA_UP_genes.xlsx", header=T)
#genes = genes$H23_RNa.H25_RNa.WT_RNa
#genes = genes[!is.na(genes)]

##################################################################
##                   RUNNING ENRICHR FUNCTION                   ##
##################################################################
# Run the universal enrichment function from clusterProfiler using SELF-PROVIDED gene lists and gene annotation files. 


# TERM2GENE
term2gene_BP = read.table("../cache/goid2gene_BP.txt", sep="\t", header=T, quote="", stringsAsFactors = F)
term2gene_CC = read.table("../cache/goid2gene_CC.txt", sep="\t", header=T, quote="", stringsAsFactors = F)
term2gene_MF = read.table("../cache/goid2gene_MF.txt", sep="\t", header=T, quote="", stringsAsFactors = F)

# TERM2NAME
term2name_BP = read.table("../cache/godb_BP.txt", sep="\t", header=T, quote="", stringsAsFactors = F)
term2name_CC = read.table("../cache/godb_CC.txt", sep="\t", header=T, quote="", stringsAsFactors = F)
term2name_MF = read.table("../cache/godb_MF.txt", sep="\t", header=T, quote="", stringsAsFactors = F)

# Enricher function
enricher_results_BP = enricher(gene = genes, # a vector of gene id
               TERM2GENE = term2gene_BP, # user input annotation of TERM TO GENE mapping 
               TERM2NAME = term2name_BP, # user input of TERM TO NAME mapping
               pvalueCutoff = 0.05, # p-value cutoff (default)
               pAdjustMethod = "BH", # multiple testing correction method to calculate adjusted p-value (default)

              )

enricher_results_CC = enricher(gene = genes, # a vector of gene id
                            TERM2GENE = term2gene_CC, # user input annotation of TERM TO GENE mapping 
                            TERM2NAME = term2name_CC, # user input of TERM TO NAME mapping
                            pvalueCutoff = 0.05, # p-value cutoff (default)
                            pAdjustMethod = "BH", # multiple testing correction method to calculate adjusted p-value (default)
)

enricher_results_MF = enricher(gene = genes, # a vector of gene id
                            TERM2GENE = term2gene_MF, # user input annotation of TERM TO GENE mapping 
                            TERM2NAME = term2name_MF, # user input of TERM TO NAME mapping
                            pvalueCutoff = 0.05, # p-value cutoff (default)
                            pAdjustMethod = "BH", # multiple testing correction method to calculate adjusted p-value (default)
)

enricher_df_BP = as.data.frame(enricher_results_BP)
enricher_df_CC = as.data.frame(enricher_results_CC)
enricher_df_MF = as.data.frame(enricher_results_MF)

write.table(enricher_df_BP, "../../Resultados de las comparaciones/5. COMPARACIONES (H23-SOS1-WT vs sus controles)/ENRICHMENT/Hoja_UP_H23_H25_Comunes/enricher_BP_df.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(enricher_df_CC, "../../Resultados de las comparaciones/5. COMPARACIONES (H23-SOS1-WT vs sus controles)/ENRICHMENT/Hoja_UP_H23_H25_Comunes/enricher_CC_df.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(enricher_df_MF, "../../Resultados de las comparaciones/5. COMPARACIONES (H23-SOS1-WT vs sus controles)/ENRICHMENT/Hoja_UP_H23_H25_Comunes/enricher_MF_df.txt", sep="\t", row.names=FALSE, quote=FALSE)



#################################################################
##                  RUNNING ENRICHGO FUNCTION                  ##
#################################################################


# Find rice orgDb
hub = AnnotationHub()
query(hub, c("Oryza sativa","orgdb"))
rice = hub[["AH114587"]]

# ID conversion table
IDtable = read.csv("../input/riceIDtable.csv")

# RAPDB IDs to Entrez ID
genes_annotated = IDtable[match(genes, IDtable$rapdb), "entrezgene"]
genes_annotated = as.character(genes_annotated[!is.na(genes_annotated)])


# EnrichGO function
enrichGO_results_BP = enrichGO(gene = genes_annotated, # a vector of gene id
                             OrgDb = rice, # OrgDb object
                             ont = "BP", # One of "MF", "BP", and "CC" subontologies
                             pvalueCutoff = 0.05, # p-value cutoff (default)
                             pAdjustMethod = "BH", # multiple testing correction method to calculate adjusted p-value (default)
)

enrichGO_results_CC = enrichGO(gene = genes_annotated, # a vector of gene id
                            OrgDb = rice, # OrgDb object
                            ont = "CC", # One of "MF", "BP", and "CC" subontologies
                            pvalueCutoff = 0.05, # p-value cutoff (default)
                            pAdjustMethod = "BH", # multiple testing correction method to calculate adjusted p-value (default)
)

enrichGO_results_MF = enrichGO(gene = genes_annotated, # a vector of gene id
                            OrgDb = rice, # OrgDb object
                            ont = "MF", # One of "MF", "BP", and "CC" subontologies
                            pvalueCutoff = 0.05, # p-value cutoff (default)
                            pAdjustMethod = "BH", # multiple testing correction method to calculate adjusted p-value (default)
)

enrichGO_df_BP = as.data.frame(enrichGO_results_BP)
enrichGO_df_CC = as.data.frame(enrichGO_results_CC)
enrichGO_df_MF = as.data.frame(enrichGO_results_MF)

write.table(enrichGO_df_BP, "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/DOWN/EnrichGO/enrichGO_BP_df.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(enrichGO_df_CC, "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/DOWN/EnrichGO/enrichGO_CC_df.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(enrichGO_df_MF, "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/UP/EnrichGO/enrichGO_MF_df.txt", sep="\t", row.names=FALSE, quote=FALSE)



##################################################################
##                      RUNNING ENRICHKEGG                      ##
##################################################################

# ID conversion table
IDtable_KEGG = as.data.frame(fread("../cache/IRGSP-1.0_representative_annotation_2021-11-11.tsv", quote=""))

# RAPDB IDs to transcript ID (KEGG ID)
genes_KEGG = IDtable_KEGG[match(genes, IDtable_KEGG$Locus_ID),"Transcript_ID"]

# EnrichKEGG function
enrichKEGG_results = enrichKEGG(gene = genes_KEGG, # a vector of gene id
                   organism = "dosa", # kegg organism
                   keyType = "kegg", # keytype of input gene
                   pvalueCutoff = 0.05, # p-value cutoff (default)
                   pAdjustMethod = "BH", # multiple testing correction method to calculate adjusted p-value (default)
)

enrichKEGG_df = as.data.frame(enrichKEGG_results)
write.table(enrichKEGG_df, "../../Resultados de las comparaciones/5. COMPARACIONES (H23-SOS1-WT vs sus controles)/ENRICHMENT/Hoja_UP_H23_H25_Comunes/enrichKEGG_df.txt", sep="\t", row.names=FALSE, quote=FALSE)
