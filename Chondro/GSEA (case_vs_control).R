custom_gene_chondro_df <- read.table(file = "~/Desktop/GitHub/Chondro/case_vs_control/case_vs_control.tsv", 
                                     sep = "\t",
                                     header = TRUE)
custom_gene_chondro_df_s <- setNames(custom_gene_chondro_df$cluster, custom_gene_chondro_df$gene_name)
#создаем ранк для второй группы (case_vs_control) 
vres_chondro_2 <- as.data.frame(vres_chondro_2)
vres_chondro_2$gene_id <- rownames(vres_chondro_2)
data_set <- read.table(file = "~/Downloads/ensbl2geneid.tsv", header = TRUE, sep = "\t")
vres_chondro_2 <- merge(vres_chondro_2, data_set, by = 'gene_id')

vres_chondro_2_LFC <- vres_chondro_2$log2FoldChange
vres_chondro_2pvalue <- -log10(vres_chondro_2$pvalue)
vres_chondro_2_LFC <- sign(vres_chondro_2_LFC)
rank_chondro_2 <- vres_chondro_2_LFC * vres_chondro_2pvalue

rank_chondro_2 <- sort(rank_chondro_2, decreasing = TRUE)
rank_chondro_2 <- data.frame(gene_name = vres_chondro_2$gene_name, rank_chondro_2)
rank_chondro_2_df <- rank_chondro_2
rank_chondro_2$gene_name <- NULL

rank_chondro_2 <- setNames(rank_chondro_2$rank_chondro_2, as.character(rank_chondro_2_df$gene_name))
rank_chondro_2 <- rank_chondro_2[!duplicated(names(rank_chondro_2))]

library(clusterProfiler)
library(biomaRt)
library(GseaVis)
gsea_1_chondro <- GSEA(geneList = rank_chondro_2, TERM2GENE = custom_gene_chondro_df[c(2,1)], 
                       eps = 0, pAdjustMethod = "fdr",
                       pvalueCutoff = 1)
gsea_1_chondro.df <- as.data.frame(gsea_1_chondro)
gseaNb(gsea_1_chondro, geneSetID = gsea_1_chondro.df$ID[1])
gseaNb(gsea_1_chondro, geneSetID = gsea_1_chondro.df$ID[2])
gseaNb(gsea_1_chondro, geneSetID = gsea_1_chondro.df$ID[3])

#применение сырого датасета для GSEA
vres_chondro_1_df <- vres_chondro_1

vres_chondro_1_df <- as.data.frame(vres_chondro_1_df)
vres_chondro_1_df$gene_id <- rownames(vres_chondro_1_df)
vres_chondro_1_df <- merge(vres_chondro_1_df, data_set, by = 'gene_id')

vres_chondro_1_df_LFC <- vres_chondro_1_df$log2FoldChange
vres_chondro_1_df_pvalue <- - log10(vres_chondro_1_df$pvalue)
vres_chondro_1_df_LFC <- sign(vres_chondro_1_df_LFC)
rank <-  vres_chondro_1_df_pvalue * vres_chondro_1_df_LFC
rank <- sort(rank, decreasing = TRUE)

rank <- data.frame(gene_name = vres_chondro_1_df$gene_name, rank)
rank_df <- rank
rank$gene_name <- NULL

rank <- setNames(rank$rank, as.character(rank_df$gene_name))
rank <- rank[!duplicated(names(rank))]

#применение втроенной кластеризации
library(enrichplot) # на 22 версии BiocManager (аналогично и дргуие пакеты)
library(clusterProfiler)
library(org.Hs.eg.db) # база данных человеческого организма
library(DOSE)
#GO
rank_for_GO <- rank_df$gene_name # нужны именованные гены в формате character (НЕ factor)
#биологические процессы
ggo_bp <- groupGO(gene = rank_for_GO,
               OrgDb = org.Hs.eg.db,
               ont = 'BP',
               level = 1,
               readable = FALSE)

ego_bp <- enrichGO(
  gene          = rank_for_GO,              # значимые гены (кастомный набор через ранки/сырой датасет)
  universe      = vres_chondro_1_df$gene_name,   # фоновый набор (все гены эксперимента)
  OrgDb         = org.Hs.eg.db,      # база данных организма
  ont           = "BP",              # онтология: "BP", "MF", "CC"
  pAdjustMethod = "BH",              # метод коррекции (Benjamini-Hochberg)
  pvalueCutoff  = 0.05,              # порог p-value
  qvalueCutoff  = 0.25,              # порог q-value
  keyType = 'SYMBOL',
  readable      = FALSE
)

#молекулярные функции
ggo_mf <- groupGO(gene = rank_for_GO,
                  OrgDb = org.Hs.eg.db,
                  ont = 'MF',
                  level = 1,
                  readable = FALSE)

ego_mf <- enrichGO(
  gene          = rank_for_GO,              # значимые гены (кастомный набор через ранки/сырой датасет)
  universe      = vres_chondro_1_df$gene_name,   # фоновый набор (все гены эксперимента)
  OrgDb         = org.Hs.eg.db,      # база данных организма
  ont           = "MF",              # онтология: "BP", "MF", "CC"
  pAdjustMethod = "BH",              # метод коррекции (Benjamini-Hochberg)
  pvalueCutoff  = 0.05,              # порог p-value
  qvalueCutoff  = 0.25,              # порог q-value
  keyType = 'SYMBOL',
  readable      = FALSE
)

#клеточные компоненты
ggo_cc <- groupGO(gene = rank_for_GO,
                  OrgDb = org.Hs.eg.db,
                  ont = 'CC',
                  level = 1,
                  readable = FALSE)

ego_cc <- enrichGO(
  gene          = rank_for_GO,              # значимые гены (кастомный набор через ранки/сырой датасет)
  universe      = vres_chondro_1_df$gene_name,   # фоновый набор (все гены эксперимента)
  OrgDb         = org.Hs.eg.db,      # база данных организма
  ont           = "CC",              # онтология: "BP", "MF", "CC"
  pAdjustMethod = "BH",              # метод коррекции (Benjamini-Hochberg)
  pvalueCutoff  = 0.05,              # порог p-value
  qvalueCutoff  = 0.25,              # порог q-value
  keyType = 'SYMBOL',
  readable      = FALSE
)

#применение GO_GSEA
ego_precise_MF <- gseGO(
  geneList     = rank,
  OrgDb        = org.Hs.eg.db,
  ont          = "MF",
  minGSSize    = 15,      # меньше минимальный размер
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH", 
  keyType = 'SYMBOL',
  verbose      = FALSE
)
ego_precise_MF_df <- as.data.frame(ego_precise_MF)
goplot(ego_precise_MF) # создание кластера

ego_precise_BP <- gseGO(
  geneList     = rank,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  minGSSize    = 15,      
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH", 
  keyType = 'SYMBOL',
  verbose      = FALSE
)
ego_precise_BP_df <- as.data.frame(ego_precise_BP)
goplot(ego_precise_BP)

ego_precise_СС <- gseGO(
  geneList     = rank,
  OrgDb        = org.Hs.eg.db,
  ont          = "CC",
  minGSSize    = 15,      
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  keyType = 'SYMBOL', 
  verbose      = FALSE
)
ego_precise_CC_df <- as.data.frame(ego_precise_СС)
goplot(ego_precise_СС)

#KEGG/подготавливаем именованный ранк
human <- search_kegg_organism('hsa', by='kegg_code')
print(human)
vres_chondro_1_df_LFC <- sign(vres_chondro_1_df$log2FoldChange)
vres_chondro_1_df_pvalue <- -log10(vres_chondro_1_df$pvalue)
kegg_rank <-  vres_chondro_1_df_pvalue * vres_chondro_1_df_LFC
kegg_rank <- sort(kegg_rank, decreasing = TRUE)
kegg_rank <- data.frame(gene_name = vres_chondro_1_df$gene_name, kegg_rank)
rank_for_kegg_names_values_df <- bitr(kegg_rank$gene_name, 
                       fromType = "SYMBOL",
                       toType   = "ENTREZID", 
                       OrgDb    = org.Hs.eg.db)
colnames(rank_for_kegg_names_values_df)[1] <- 'gene_name'
kegg_rank_prefinal <- merge(rank_for_kegg_names_values_df, kegg_rank, by = 'gene_name')
kegg_rank_prefinal <- kegg_rank_prefinal[order(kegg_rank_prefinal$kegg_rank, decreasing = TRUE),] 
kegg_rank_final <- setNames(kegg_rank_prefinal$kegg_rank, as.character(kegg_rank_prefinal$ENTREZID))
kegg_rank_final <- kegg_rank_final[!duplicated(names(kegg_rank_prefinal))]
#KEGG/GSEA
kegg_analysis_GSEA <- gseKEGG(geneList = kegg_rank_final,
               organism     = 'hsa',
               minGSSize    = 10,
               maxGSSize = 800,
               pvalueCutoff = 0.1,
               pAdjustMethod = 'BH',
               verbose = FALSE)
kegg_analysis_GSEA_df <- as.data.frame(kegg_analysis_GSEA)
#no term enriched under specific pvalueCutoff...(не дал результатов при pvalue = 0.05, но дал при pvalue = 1)
browseKEGG(kegg_analysis_GSEA_df, 'hsa03420')
#WikiPathways
wp_gsea <- gseWP(
  geneList = kegg_rank_final,
  organism = "Homo sapiens",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05,
  verbose = FALSE)
wp_gsea_df <- as.data.frame(wp_gsea)

#Reactome pathway gene set enrichment analysis
library(ReactomePA)
reactome_gsea <- gsePathway(
  geneList = kegg_rank_final,
  organism = "human",
  minGSSize = 10,
  maxGSSize = 800,
  pvalueCutoff = 0.1,
  pAdjustMethod = "BH",
  verbose = FALSE
)
reactome_gsea <- as.data.frame(reactome_gsea)
#no term enriched under specific pvalueCutoff...(не дал результатов вообще)
