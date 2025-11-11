library(openxlsx)
#case_vs_control (UR)
df_chondro_case_vs_control <- subset(res_chondro_2, res_chondro_2$log2FoldChange > 0)
tail(df_chondro_case_vs_control)
df_chondro_case_vs_control_UR_1 <- read.xlsx(sep = '\t', 
                                              xlsxFile = "~/Desktop/GitHub/Chondro/case_vs_control/UR/(+)Mixed, incl. Amplification of signal from the kinetochores, and Condensin complex.xlsx",
                                             colNames = TRUE)
colnames(df_chondro_case_vs_control_UR_1)[1] <- 'node1'
gene_list_cluster_1_UR <- c(unique(df_chondro_case_vs_control_UR_1$node1), unique(df_chondro_case_vs_control_UR_1$node2))
gene_list_cluster_1_UR <- as.vector(gene_list_cluster_1_UR)
gene_list_cluster_1_UR <- as.data.frame(gene_list_cluster_1_UR)
gene_list_cluster_1_UR$cluster <- 1 : nrow(gene_list_cluster_1_UR)
gene_list_cluster_1_UR$cluster <- 'cl1'
names(gene_list_cluster_1_UR)[names(gene_list_cluster_1_UR) == 'gene_list_cluster_1_UR'] <- 'gene_name'
#ссылки на кластеры
#https://string-db.org/cgi/network?taskId=bltSV1cYQxMC&sessionId=bAs0oUG8zjlZ
#https://string-db.org/cgi/network?taskId=bMlKWMW118L2&sessionId=bAs0oUG8zjlZ

#case_vs_control (DR)
df_chondro_case_vs_control <- subset(res_chondro_2, res_chondro_2$log2FoldChange < 0)
head(df_chondro_case_vs_control)
#https://string-db.org/cgi/network?taskId=bXuWK5RhDUnU&sessionId=bAs0oUG8zjlZ
#https://string-db.org/cgi/network?taskId=bAqHpAJbFutW&sessionId=bAs0oUG8zjlZ
df_chondro_case_vs_control_DR_1 <- read.xlsx(sep = '\t', xlsxFile = "~/Desktop/GitHub/Chondro/case_vs_control/DR/(-)ECM-receptor interaction+Extracellular matrix organization.xlsx", 
                                             colNames = TRUE)
colnames(df_chondro_case_vs_control_DR_1)[1] <- 'node1'
gene_list_cluster_1_DR <- c(unique(df_chondro_case_vs_control_DR_1$node1), unique(df_chondro_case_vs_control_DR_1$node2))
gene_list_cluster_1_DR <- as.vector(gene_list_cluster_1_DR)
gene_list_cluster_1_DR <- as.data.frame(gene_list_cluster_1_DR)
gene_list_cluster_1_DR$cluster <- 1 : nrow(gene_list_cluster_1_DR)
gene_list_cluster_1_DR$cluster <- 'cl2'
names(gene_list_cluster_1_DR)[names(gene_list_cluster_1_DR) == 'gene_list_cluster_1_DR'] <- 'gene_name'

# ссылка на другие кластеры в области DR
#file:///Users/egorka.1veretinicloud.com/Downloads/string_vector_graphic.svg 
#https://string-db.org/cgi/network?taskId=bXuWK5RhDUnU&sessionId=bAs0oUG8zjlZ
#https://string-db.org/cgi/network?taskId=bAqHpAJbFutW&sessionId=bAs0oUG8zjlZ

df_chondro_case_vs_control_DR_2 <- read.xlsx(sep = '\t', xlsxFile = "~/Desktop/GitHub/Chondro/case_vs_control/DR/(-)Interferon alphabeta signaling+Defense response to virus+Antiviral defense.xlsx", 
                                             colNames = TRUE)
colnames(df_chondro_case_vs_control_DR_2)[1] <- 'node1'
gene_list_cluster_2_DR <- c(unique(df_chondro_case_vs_control_DR_2$node1), unique(df_chondro_case_vs_control_DR_2$node2))
gene_list_cluster_2_DR <- as.vector(gene_list_cluster_2_DR)
gene_list_cluster_2_DR <- as.data.frame(gene_list_cluster_2_DR)
gene_list_cluster_2_DR$cluster <- 1 : nrow(gene_list_cluster_2_DR)
gene_list_cluster_2_DR$cluster <- 'cl3'
names(gene_list_cluster_2_DR)[names(gene_list_cluster_2_DR) == 'gene_list_cluster_2_DR'] <- 'gene_name'

#выполнение рассчетов
custom_gen_set_chondro_1 <- rbind(gene_list_cluster_1_UR, gene_list_cluster_1_DR, gene_list_cluster_2_DR)
write.table(custom_gen_set_chondro_1, "~/Desktop/GitHub/Chondro/case_vs_control.tsv", quote = FALSE, sep = "\t")
