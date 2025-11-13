# импорт датасета
library(ggplot2)
library(vegan)
library(DESeq2)
library(openxlsx)
meta_data_chondro <- read.xlsx("~/Desktop/GitHub/Chondro/chondro_metadata/chondro_metadata_1.xlsx", sheet= T, sep="\t")

meta_data_chondro_a <- meta_data_chondro
colnames(meta_data_chondro_a)

# Перевести значения в factor()
meta_data_chondro_a$sampleid <- factor(meta_data_chondro_a$sampleid)
meta_data_chondro_a$File <- factor(meta_data_chondro_a$File)
meta_data_chondro_a$cell <- factor(meta_data_chondro_a$cell)
meta_data_chondro_a$type <- factor(meta_data_chondro_a$type)
meta_data_chondro_a$paired <- factor(meta_data_chondro_a$paired)
meta_data_chondro_a$group <- factor(meta_data_chondro_a$group)
meta_data_chondro_a$run <- factor(meta_data_chondro_a$run)

# фильтруем метаданные для оценки case or control
meta_data_chondro_a <- subset(meta_data_chondro_a, meta_data_chondro_a$group == 'control' | 
                                meta_data_chondro_a$group == 'case')

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

ggplot(pca_data_1, aes(x = PC1, y = PC2, color = cell, shape = type, )) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  scale_shape_manual(values = 1:length(unique(pca_data_1$type)))

norm_counts <- counts(ddsHTSeq, normalized = TRUE)
summary(norm_counts)

# Получение результатов и сравнение выбранных групп
res_chondro_1 <- results(ddsHTSeq, contrast = c("type", "2d", "3d")) 
res_chondro_1 <- as.data.frame(res_chondro_1)

res_chondro_2 <- results(ddsHTSeq, contrast = c("group", "case", "control")) 
res_chondro_2 <- as.data.frame(res_chondro_2)

# Для GSEA
vres_chondro_1 <- res_chondro_1
vres_chondro_2 <- res_chondro_2

# Сортировка значений по уменьшению параметра LFC

res_chondro_1 <- res_chondro_1[order(res_chondro_1$log2FoldChange, decreasing = T),]
res_chondro_2 <- res_chondro_2[order(res_chondro_2$log2FoldChange, decreasing = T),]

# Сортировка значений по заданным условиям
res_chondro_1 <- subset(res_chondro_1, (res_chondro_1$log2FoldChange > 1 | res_chondro_1$log2FoldChange < -1) & res_chondro_1$padj < 0.05)
res_chondro_2 <- subset(res_chondro_2, (res_chondro_2$log2FoldChange > 1 | res_chondro_2$log2FoldChange < -1) & res_chondro_2$padj < 0.05)


# Графическая часть модуль 1
res_chondro_1$regulation <- ifelse(
  is.na(res_chondro_1$padj) | res_chondro_1$padj >= 0.05 | abs(res_chondro_1$log2FoldChange) <= 1, # Выбор условий
  "Not significant",
  ifelse(res_chondro_1$log2FoldChange > 0, "Upregulated", "Downregulated")
)

library(volcanoPlot)
library(ggplot2)

volcano_plot_1 <- ggplot(res_chondro_1, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = regulation), size = 2) + 
  scale_color_manual(
    name = "Regulation",  # Название легенды
    values = c("Upregulated" = "red", "Downregulated" = "green", "Not significant" = "gray")
  ) +
  theme_minimal() +
  labs(
    title = "2d vs 3d",
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

# Графическая часть модуль 2
res_chondro_2$regulation <- ifelse(
  is.na(res_chondro_2$padj) | res_chondro_2$padj >= 0.05 | abs(res_chondro_2$log2FoldChange) <= 1, # Выбор условий
  "Not significant",
  ifelse(res_chondro_2$log2FoldChange > 0, "Upregulated", "Downregulated")
)

volcano_plot_2 <- ggplot(res_chondro_2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = regulation), size = 2) + 
  scale_color_manual(
    name = "Regulation",  # Название легенды
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

print(volcano_plot_2)

# сохранение данных по группам 2d_vs_3d and case_vs_control
write.table(cbind(gene_id = rownames(res_chondro_1), res_chondro_1), 
            file = "~/Desktop/GitHub/Chondro/res_chondro_1_2d_vs_3d.tsv", sep = "\t", quote = F)
write.table(cbind(gene_id = rownames(res_chondro_2), res_chondro_2),
            file = "~/Desktop/GitHub/Chondro/res_chondro_2_case_vs_control_3d.tsv", sep = "\t", quote = F)

# датасет 1
# > 0
data_vres_1_1 = subset(vres_chondro_1, vres_chondro_1$log2FoldChange > 0)
head(data_vres_1_1)
# < 0
data_vres_1_2 = subset(vres_chondro_1, vres_chondro_1$log2FoldChange < 0)
head(data_vres_1_2)

# датасет 2
# > 0
data_vres_2_1 = subset(vres_chondro_2, vres_chondro_2$log2FoldChange > 0)
head(data_vres_2_1)
# < 0
data_vres_2_2 = subset(vres_chondro_2, vres_chondro_2$log2FoldChange < 0)
head(data_vres_2_2)
#резервная копия
vres_chondro_1_df <- vres_chondro_1


