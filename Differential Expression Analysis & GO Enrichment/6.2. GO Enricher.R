#################################################################
##                   Enricher representation                   ##
#################################################################


#################################################################
##                           Dotplot                           ##
#################################################################
#BP
figure1 = dotplot(enricher_results_BP,
                  font.size = 12,
                  showCategory = 20)
ggsave(figure1,
       filename = "../../Resultados de las comparaciones/5. COMPARACIONES (H23-SOS1-WT vs sus controles)/ENRICHMENT/RAIZ_DOWN_H23_ESP_dotplot_BP.png",
       height = 20,width = 25,units = "cm") 

figure2 = barplot(enricher_results_BP,
                  font.size = 12,
                  showCategory = 20)
ggsave(figure2,
       filename = "../../Resultados de las comparaciones/5. COMPARACIONES (H23-SOS1-WT vs sus controles)/ENRICHMENT/RAIZ_DOWN_H23_ESP_barplot_BP.png",
       height = 22,width = 25,units = "cm") 

#CC
figure2 = dotplot(enricher_results_CC,
                  title = "Top 10 most statistically significant enriched GO terms (CC)",
                  font.size = 10,
                  showCategory = 10)
ggsave(figure2,
       filename = "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/UP/Enricher/dotplot_CC.png",
       height = 14,width = 24,units = "cm") 

#MF
figure3 = dotplot(enricher_results_MF,
                  title = "Top 20 most statistically significant enriched GO terms (MF)",
                  font.size = 10,
                  showCategory = 20)
ggsave(figure3,
       filename = "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/UP/Enricher/dotplot_MF.png",
       height = 14,width = 20,units = "cm") 




#################################################################
##                           Cnetplot                          ##
#################################################################

#BP
figure1 = cnetplot(enricher_results_BP,
                  font.size = 10,
                  showCategory = 20,
                  node_label = "all",
                  cex_label_category = 0.6,
                  cex_label_gene = 0.25)
ggsave(figure1,
       filename = "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/UP/Enricher/cnetplot_BP.png",
       height = 16,width = 24,units = "cm") 


figure2 = cnetplot(enricher_results_CC,
                  font.size = 10,
                  showCategory = 10,
                  node_label = "all",
                  cex_label_category = 0.5,
                  cex_label_gene = 0.35)
ggsave(figure2,
       filename = "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/UP/Enricher/cnetplot_CC.png",
       height = 16,width = 24,units = "cm") 


figure3 = cnetplot(enricher_results_MF,
                  font.size = 10,
                  showCategory = 20,
                  node_label = "all",
                  cex_label_category = 0.6,
                  cex_label_gene = 0.25)
ggsave(figure3,
       filename = "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/UP/Enricher/cnetplot_MF.png",
       height = 16,width = 24,units = "cm") 


#################################################################
##                           Heatplot                          ##
#################################################################

df = read.table(file = "../input/activated_df.txt", header = T, sep = "\t", row.names = 1)
fold_changes = df$LogFC
 
names(fold_changes) = row.names(df)

#BP
figure1 = heatplot(enricher_results_BP, showCategory = 20) #foldChange = fold_changes)

ggsave(figure1,
       filename = "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/UP/Enricher/heatplot_BP.png",
       height = 14,width = 50,units = "cm") 
 
# #CC
figure2 = heatplot(enricher_results_CC, showCategory = 10, foldChange = fold_changes)

ggsave(figure2,
      filename = "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/UP/Enricher/heatplot_CC.png",
      height = 12,width = 40,units = "cm") 

#MF
figure3 = heatplot(enricher_results_MF, showCategory = 20, foldChange = fold_changes)

ggsave(figure3,
      filename = "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/UP/Enricher/heatplot_MF.png",
      height = 12,width = 50,units = "cm") 

#################################################################
##                           UPsetplot                         ##
#################################################################

#BP
figure1 = upsetplot(enricher_results_BP,
                  font.size = 10,
                  showCategory = 20)
ggsave(figure1,
       filename = "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/UP/Enricher/upsetplot_BP.png",
       height = 16,width = 20,units = "cm") 

#CC
figure2 = upsetplot(enricher_results_CC,
                  font.size = 10,
                  showCategory = 10)
ggsave(figure2,
       filename = "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/UP/Enricher/upsetplot_CC.png",
       height = 12,width = 20,units = "cm") 

#MF
figure3 = upsetplot(enricher_results_MF,
                  font.size = 10,
                  showCategory = 20)
ggsave(figure3,
       filename = "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/UP/Enricher/upsetplot_MF.png",
       height = 16,width = 40,units = "cm") 


#################################################################
##                           Treeplot                          ##
#################################################################
object_emapplot_BP = pairwise_termsim(enricher_results_BP)
object_emapplot_CC = pairwise_termsim(enricher_results_CC)
object_emapplot_MF = pairwise_termsim(enricher_results_MF)


figure1 = treeplot(object_emapplot_BP)
ggsave(figure1,
       filename = "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/UP/Enricher/treemap_BP.png",
       height = 16,width = 40,units = "cm") 

figure2 = treeplot(object_emapplot_CC)
ggsave(figure2,
       filename = "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/UP/Enricher/treemap_CC.png",
       height = 16,width = 40,units = "cm") 

figure3 = treeplot(object_emapplot_MF)
ggsave(figure3,
       filename = "../../Resultados de las comparaciones/3. ENRIQUECIMIENTOS DE INTERES/Genes específicos SOS1 Raíz/UP/ENricher/treemap_MF.png",
       height = 16,width = 40,units = "cm") 
