library(clusterProfiler)
library(openxlsx)
library(DESeq2)
library(vegan)
library(ggplot2)
library(BioNERO)
library(biomaRt)
library(GseaVis)
# meta data import
meta_data_case_vs_control <- read.xlsx(xlsxFile = "/home/veretinegor/R/metachondro.xlsx", sheet = T, sep = '\t')
head(meta_data_case_vs_control)
# DESeq2
meta_data_analisys <- meta_data_case_vs_control
colnames(meta_data_analisys)
meta_data_analisys <- subset(meta_data_analisys, meta_data_analisys$group == 'case' | meta_data_analisys$group == 'control')
meta_data_analisys[] <- lapply(meta_data_analisys, factor)
meta_data_analisys$group <- relevel(meta_data_analisys$group, ref = 'control')

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = meta_data_analisys,
                                       directory = "/home/veretinegor/R/chondro_htseq",
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
# COexpression
set.seed(123)

norm_counts.t <- as.data.frame(t(norm_counts))

final_exp <- exp_preprocess(
  t(norm_counts.t), min_exp = 10, variance_filter = TRUE, n = 3000, 
)

mem.maxVSize(vsize = 16000)

sft.p <- SFT_fit(final_exp, net_type = "signed hybrid", cor_method = "pearson")
power.p <- sft.p$power

net.p <- exp2gcn(
  final_exp, net_type = "signed hybrid", SFTpower = power.p, 
  cor_method = "pearson"
)

plot_cor_module.p <- plot_eigengene_network(net.p)
plot_cor_module.p
plot_ngenes.p <- plot_ngenes_per_module(net.p)
plot_ngenes.p
meta_data_case_vs_control$run <- NULL
summary(meta_data_case_vs_control)
colnames(meta_data_case_vs_control)
dfg <- meta_data_case_vs_control[c(1,6)]
rownames(dfg) <- dfg$sampleid
dfg <- dfg[-1]

MEtrait.p <- module_trait_cor(exp = final_exp, MEs = net.p$MEs, metadata = dfg)
module_trait_cor_plot.p <- plot_module_trait_cor(MEtrait.p)
module_trait_cor_plot.p

# Построение графиков коэкспрессии
colnames(final_exp)
colnames(dfg)
dfg <- subset(dfg, dfg$group == 'case' | dfg$group == 'control')
plot_expression_profile(
  exp = final_exp, 
  net = net.p, 
  plot_module = TRUE,
  metadata = dfg,
  modulename = "green"
)

plot_expression_profile(
  exp = final_exp, 
  net = net.p, 
  plot_module = TRUE,
  metadata = dfg,
  modulename = "greenyellow"
)
summary(net.p)
net_Mes_chondro <- data.frame(net.p$MEs)
net_GM_chondro <- data.frame(net.p$genes_and_modules)
write.table(net_Mes_chondro, file = "/home/veretinegor/R/net_Mes_chondro.tsv", sep ="\t", quote = F)
write.table(net_GM_chondro, file = "/home/veretinegor/R/net_GM_chondro.tsv", sep ="\t", quote = F)
data_genes <- read.csv("/home/veretinegor/R/ensbl2geneid.tsv", header = T, sep = "\t")
net_GM_chondro$gene_id <- net_GM_chondro$Genes
net_GM_chondro$Genes <- NULL
data_chondro <- merge(net_GM_chondro, data_genes, by = 'gene_id')
summary(data_chondro)
data_chondro$description <- NULL
# Создаем группу целевых модулей, для более глубокого анализа необходимо разделить модули!
data_group <- subset(data_chondro, data_chondro$Modules == 'greenyellow' | data_chondro$Modules == 'green')
data_green_chondro <- subset(data_group, data_group$Modules == 'green')
write.table(x = data_green_chondro, file = "/home/veretinegor/R/data_green_chondro.tsv", sep = "\t", quote = F)
data_greenyellow_chondro <- subset(data_group, data_group$Modules == 'greenyellow')
write.table(x = data_greenyellow_chondro, file = "/home/veretinegor/R/data_yellowgreen_chondro.tsv", sep = "\t", quote = F)
# GSEA
head(vres_chondro_GG)
df_analysis <- setNames(object = net_GM_chondro$Modules, nm = net_GM_chondro$gene_id)
colnames(data_genes)

res_GSEA <- merge(vres_chondro_GG, data_genes, by = 'gene_id')
res_GSEA_LFC <- sign(res_GSEA$log2FoldChange)
res_GSEA_pvalue <- -log10(res_GSEA$pvalue)
rank_GSEA <- res_GSEA_LFC * res_GSEA_pvalue
rank_GSEA <- sort(rank_GSEA, decreasing = T)
rank_GSEA_df <- data.frame(gene_name = res_GSEA$gene_name, rank_GSEA)

rank_final_GSEA <- setNames(object = rank_GSEA_df$rank_GSEA, nm = as.factor(rank_GSEA_df$gene_name))
print(length(rank_final_GSEA))

rank_final_GSEA <- rank_final_GSEA[!duplicated(names(rank_final_GSEA))]
print(length(rank_final_GSEA))

gsea_green <- GSEA(geneList = rank_final_GSEA, TERM2GENE = data_green_chondro[c(2,3)], 
                   eps = 0, pAdjustMethod = "fdr")
gsea_df_green <- data.frame(gsea_green)
write.table(x = gsea_df_green,file = "/home/veretinegor/R/green_module_coex.tsv")
gseaNb(gsea_green, geneSetID = gsea_df_green$ID[1])

gsea_greenyellow <- GSEA(geneList = rank_final_GSEA, TERM2GENE = data_greenyellow_chondro[c(2,3)], 
                         eps = 0, pAdjustMethod = "fdr", pvalueCutoff = 0.5)
gsea_df_greenyellow <- data.frame(gsea_greenyellow)
write.table(x = gsea_df_greenyellow,file = "/home/veretinegor/R/greenyellow_module_coex.tsv")
gseaNb(gsea_greenyellow, geneSetID = gsea_greenyellow$ID[1])




