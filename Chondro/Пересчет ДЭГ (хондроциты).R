# импорт датасета
library(ggplot2)
library(vegan)
library(DESeq2)
library(openxlsx)
library(volcanoPlot)
library(clusterProfiler)
library(BioNERO)
library(biomaRt)
library(GseaVis)
library(enrichplot) # на 22 версии BiocManager (аналогично и дргуие пакеты)
library(org.Hs.eg.db) # база данных человеческого организма
library(DOSE)
meta_data_chondro <- read.xlsx("~/Desktop/GitHub/Chondro/chondro_metadata/chondro_metadata_1.xlsx", 
                               sheet= T, sep="\t")
meta_data_chondro_a <- meta_data_chondro
colnames(meta_data_chondro_a)

# Перевести значения в factor()
meta_data_chondro_a <- subset(meta_data_chondro_a, meta_data_chondro_a$group == 'control' | 
                                meta_data_chondro_a$group == 'case')
meta_data_chondro_a[] <- lapply(meta_data_chondro_a, factor)
summary(meta_data_chondro_a)
meta_data_chondro_a$group <- relevel(meta_data_chondro_a$group, ref = "control")
# фильтруем метаданные для оценки case or control
# DEseq2
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = meta_data_chondro_a,
                                       directory = "~/Desktop/chondro_htseq",
                                       design = ~ run + type + group)
rownames(ddsHTSeq) <- sapply(strsplit(rownames(ddsHTSeq), "\\."), function(x) x[1])

# фильтрация данных
ddsHTSeq <- DESeq(ddsHTSeq)
colData(ddsHTSeq)
min_count <- 10
min_samples <- ncol(ddsHTSeq) / 3 * 2

norm_counts <- counts(ddsHTSeq, normalized = TRUE)
keep <- rowSums(norm_counts >= min_count) >= min_samples
ddsHTSeq <- ddsHTSeq[keep, ]
summary(keep)

ddsHTSeq <- DESeq(ddsHTSeq, test = "Wald")
resultsNames(ddsHTSeq)

# Применение модуля vst для нормализации данных
vsd_1 <- vst(object = ddsHTSeq, blind = F)
vst_counts_1 <- assay(vsd_1)

sample_info_1 <- as.data.frame(colData(ddsHTSeq))
dist_matrix_1 <- vegdist(t(vst_counts_1), method = "euclidean")
permanova_result_1 <- adonis2(dist_matrix_1 ~ run + type + group,
                              sample_info_1, 
                              permutations = 9999, by = "margin")
print(permanova_result_1)
print(colData(ddsHTSeq))

pca_data_1 <- plotPCA(vsd_1, intgroup = c("group"), returnData = TRUE)

# Проверяем структуру
print(head(pca_data_1))
print(class(pca_data_1))

# Строим кастомизированный график
percentVar <- round(100 * attr(pca_data_1, "percentVar"))

ggplot(pca_data_1, aes(x = PC1, y = PC2, color = group, shape = cell, )) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  scale_color_manual(values = c("control" = "green", "case" = "red"))+
  theme_minimal() +
  scale_shape_manual(values = 1:length(unique(pca_data_1$group)))

norm_counts <- counts(ddsHTSeq, normalized = TRUE)
summary(norm_counts)
# Case vs Control
res_chondro_GG <- results(ddsHTSeq, contrast = c("group", "case", "control")) 
res_chondro_GG <- as.data.frame(res_chondro_GG)
res_chondro_GG <- res_chondro_GG[order(res_chondro_GG$log2FoldChange, decreasing = TRUE),]
# For GSEA

vres_chondro_GG <- res_chondro_GG
vres_chondro_GG$gene_id <- rownames(vres_chondro_GG) # (add gene_name column)

# Filtered data
res_chondro_GG <- subset(res_chondro_GG, abs(res_chondro_GG$log2FoldChange) > 1 & res_chondro_GG$padj < 0.05)
# Graphics
res_chondro_GG$regulation <- ifelse(
  is.na(res_chondro_GG$padj) | res_chondro_GG$padj >= 0.05 | abs(res_chondro_GG$log2FoldChange) <= 1,
  "Not significant",
  ifelse(res_chondro_GG$log2FoldChange > 0, "Upregulated", "Downregulated")
)

volcano_plot_1 <- ggplot(res_chondro_GG, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = regulation), size = 2) + 
  scale_color_manual(
    name = "Regulation",
    values = c("Upregulated" = "red", "Downregulated" = "green", "Not significant" = "gray")
  ) +
  theme_minimal() +
  labs(
    title = "case vs control",
    x = "Log2 Fold Change/Expression", 
    y = "-log10 Adjusted p-value"
  ) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),                   
    axis.title = element_text(size = 14),                    
    plot.title = element_text(size = 14)        
  )
print(volcano_plot_1)
# Dataset
write.table(cbind(gene_id = rownames(res_chondro_GG), res_chondro_GG),
            file = "~/Desktop/GitHub/Chondro/case_vs_control.tsv", sep = "\t", quote = F)
# Points (>= | <=)
res_chondro_GG_UR <- subset(res_chondro_GG, res_chondro_GG$log2FoldChange >= 0)
tail(res_chondro_GG_UR)

res_chondro_GG_DR <- subset(res_chondro_GG, res_chondro_GG$log2FoldChange <= 0)
head(res_chondro_GG_DR)
# STRING (custome_genes)
# UR_1
custome_UR_1 <- read.xlsx(xlsxFile = "~/Desktop/Хондроциты(проверка результатов)/UR_chondro/Mixed, incl. Amplification of signal from the kinetochores, and Condensin complex_UR.xlsx",
                          colNames = T)
colnames(custome_UR_1)[1] <- 'node1'
custome_UR_1_vector <- c(unique(custome_UR_1$node1), unique(custome_UR_1$node2))
custome_UR_1_vector <- as.vector(custome_UR_1_vector)
custome_UR_1_vector_df <- as.data.frame(custome_UR_1_vector)
custome_UR_1_vector_df$cluster <- 'cl1'
colnames(custome_UR_1_vector_df)[1] <- 'gene_name'
# DR_1
custome_DR_1 <- read.xlsx(xlsxFile = "~/Desktop/Хондроциты(проверка результатов)/DR_chondro/Extracellular matrix organization_DR_cluster_1.xlsx",
                          colNames = T)
colnames(custome_DR_1)[1] <- 'node1'
custome_DR_1_vector <- c(unique(custome_DR_1$node1), unique(custome_DR_1$node2))
custome_DR_1_vector <- as.vector(custome_DR_1_vector)
custome_DR_1_vector_df <- as.data.frame(custome_DR_1_vector)
custome_DR_1_vector_df$cluster <- 'cl2'
colnames(custome_DR_1_vector_df)[1] <- 'gene_name'
# DR_2
custome_DR_2 <- read.xlsx(xlsxFile = "~/Desktop/Хондроциты(проверка результатов)/DR_chondro/Mixed, incl. Interferon alpha:beta signaling, and Guanylate-binding protein, C-terminal_DR_cluster_2.tsv.xlsx",
                          colNames = T)
colnames(custome_DR_2)[1] <- 'node1'
custome_DR_2_vector <- c(unique(custome_DR_2$node1), unique(custome_DR_2$node2))
custome_DR_2_vector <- as.vector(custome_DR_2_vector)
custome_DR_2_vector_df <- as.data.frame(custome_DR_2_vector)
custome_DR_2_vector_df$cluster <- 'cl3'
colnames(custome_DR_2_vector_df)[1] <- 'gene_name'
# Genset
gene_set <- rbind(custome_UR_1_vector_df, custome_DR_1_vector_df, custome_DR_2_vector_df)
write.table(file = "~/Desktop/Хондроциты(проверка результатов)/case_vs_control_gene_set.tsv",
            quote = F, sep = '\t', x = gene_set)
gene_set_list <- setNames(gene_set$cluster, as.character(gene_set$gene_name))
# Data from DESeq2
vres_chondro_GG <- as.data.frame(vres_chondro_GG)
metadata_for_chondro <- read.csv(file = "~/Desktop/Хондроциты(проверка результатов)/ensbl2geneid_first.tsv",
                                 sep = '\t')
head(metadata_for_chondro)
colnames(vres_chondro_GG)
vres_chondro_GG <- merge(vres_chondro_GG, metadata_for_chondro, by = 'gene_id')
vres_chondro_GG_LFC <- sign(vres_chondro_GG$log2FoldChange)
vres_chondro_GG_pvalue <- -log10(vres_chondro_GG$pvalue)
chondro_rank <- vres_chondro_GG_LFC * vres_chondro_GG_pvalue
chondro_rank <- sort(chondro_rank, decreasing = T)
rank_df <- data.frame(gene_name = vres_chondro_GG$gene_name, chondro_rank)
rank_frame <- rank_df
rank_df$gene_name <- NULL
# Rank for GSEA
rank_for_case_vs_control <- setNames(rank_df$chondro_rank, as.character(rank_frame$gene_name))
rank_for_case_vs_control <- rank_for_case_vs_control[!duplicated(names(rank_for_case_vs_control))]
# Custome GSEA
gsea_chondro_case_vs_control <- GSEA(geneList = rank_for_case_vs_control, TERM2GENE = gene_set[c(2,1)], 
                       eps = 0, pAdjustMethod = "fdr",
                       pvalueCutoff = 0.5) # none significant
gsea_chondro_case_vs_control_df <- as.data.frame(gsea_chondro_case_vs_control)
gseaNb(gsea_chondro_case_vs_control, geneSetID = gsea_chondro_case_vs_control_df$ID[1])

# GSEA without custome
head(rank_for_case_vs_control)
length(rank_for_case_vs_control)
rank_for_GO <- rank_frame$gene_name # нужны именованные гены в формате character (НЕ factor)
# биологические процессы
ggo_bp <- groupGO(gene = rank_for_GO,
                  OrgDb = org.Hs.eg.db,
                  ont = 'BP',
                  level = 1,
                  readable = FALSE)

ego_bp <- enrichGO(
  gene          = rank_for_GO,              # значимые гены (кастомный набор через ранки/сырой датасет)
  universe      = vres_chondro_GG$gene_name,   # фоновый набор (все гены эксперимента)
  OrgDb         = org.Hs.eg.db,      # база данных организма
  ont           = "BP",              # онтология: "BP", "MF", "CC"
  pAdjustMethod = "BH",              # метод коррекции (Benjamini-Hochberg)
  pvalueCutoff  = 0.05,              # порог p-value
  qvalueCutoff  = 0.25,              # порог q-value
  keyType = 'SYMBOL',
  readable      = FALSE
)

# молекулярные функции
ggo_mf <- groupGO(gene = rank_for_GO,
                  OrgDb = org.Hs.eg.db,
                  ont = 'MF',
                  level = 1,
                  readable = FALSE)

ego_mf <- enrichGO(
  gene          = rank_for_GO,              
  universe      = vres_chondro_GG$gene_name,
  OrgDb         = org.Hs.eg.db,
  ont           = "MF",        
  pAdjustMethod = "BH",        
  pvalueCutoff  = 0.05,        
  qvalueCutoff  = 0.25,        
  keyType = 'SYMBOL',
  readable      = FALSE
)

# клеточные компоненты
ggo_cc <- groupGO(gene = rank_for_GO,
                  OrgDb = org.Hs.eg.db,
                  ont = 'CC',
                  level = 1,
                  readable = FALSE)

ego_cc <- enrichGO(
  gene          = rank_for_GO,              
  universe      = vres_chondro_GG$gene_name,
  OrgDb         = org.Hs.eg.db,      
  ont           = "CC",              
  pAdjustMethod = "BH",              
  pvalueCutoff  = 0.05,              
  qvalueCutoff  = 0.25,              
  keyType = 'SYMBOL',
  readable      = FALSE
)

# применение GO_GSEA
ego_precise_MF <- gseGO(
  geneList     = rank_for_case_vs_control,
  OrgDb        = org.Hs.eg.db,
  ont          = "MF",
  minGSSize    = 15,      # минимальный размер
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH", 
  keyType = 'SYMBOL',
  verbose      = FALSE
)
ego_precise_MF_df <- as.data.frame(ego_precise_MF)
dotplot(ego_precise_MF) # создание кластера

ego_precise_BP <- gseGO(
  geneList     = rank_for_case_vs_control,
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
dotplot(ego_precise_BP)

ego_precise_СС <- gseGO(
  geneList     = rank_for_case_vs_control,
  OrgDb        = org.Hs.eg.db,
  ont          = "CC",
  minGSSize    = 15,      
  maxGSSize    = 500,
  pvalueCutoff = 0.25,
  pAdjustMethod = "BH",
  keyType = 'SYMBOL', 
  verbose      = FALSE
)
ego_precise_CC_df <- as.data.frame(ego_precise_СС)
dotplot(ego_precise_СС)

# KEGG (подготавливаем именованный ранк)
human <- search_kegg_organism('hsa', by='kegg_code')
print(human)
rank_for_KEGG <- bitr(rank_frame$gene_name,
                      fromType = "SYMBOL",
                      toType   = "ENTREZID", 
                      OrgDb    = org.Hs.eg.db)
colnames(rank_for_KEGG)[1] <- 'gene_name'
kegg_rank_prefinal <- merge(rank_for_KEGG, rank_frame, by = 'gene_name')
kegg_rank_prefinal <- kegg_rank_prefinal[order(kegg_rank_prefinal$chondro_rank, decreasing = T),]
kegg_rank_final <- setNames(kegg_rank_prefinal$chondro_rank, as.character(kegg_rank_prefinal$ENTREZID))
kegg_rank_final <- kegg_rank_final[!duplicated(names(kegg_rank_final))]

# KEGG with GSEA
kegg_analysis_GSEA <- gseKEGG(geneList = kegg_rank_final,
                              organism     = 'hsa',
                              minGSSize    = 10,
                              maxGSSize = 800,
                              pvalueCutoff = 0.25,
                              pAdjustMethod = 'BH',
                              verbose = FALSE)
kegg_analysis_GSEA_df <- as.data.frame(kegg_analysis_GSEA)
browseKEGG(kegg_analysis_GSEA_df, 'hsa03420') # можно менять индексы из фрейма данных

# WikiPathways
options(timeout = 300)
wp_gsea <- gseWP(
  geneList = kegg_rank_final,
  organism = "Homo sapiens",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.25,
  verbose = FALSE)
wp_gsea_df <- as.data.frame(wp_gsea)

# Reactome pathway gene set enrichment analysis
library(ReactomePA)
reactome_gsea <- gsePathway(
  geneList = kegg_rank_final,
  organism = "human",
  minGSSize = 10,
  maxGSSize = 800,
  pvalueCutoff = 0.5,
  pAdjustMethod = "BH",
  verbose = FALSE
)
reactome_gsea <- as.data.frame(reactome_gsea)
