custom_gene_chondro_df <- read.table(file = "~/Desktop/GitHub/Chondro/case_vs_control/case_vs_control.tsv", 
                                     sep = "\t",
                                     header = TRUE)
custom_gene_chondro_df_s <- setNames(custom_gene_chondro_df$cluster, custom_gene_chondro_df$gene_name)
#создаем ранк
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
