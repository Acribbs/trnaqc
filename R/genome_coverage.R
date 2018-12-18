library("tidyverse")
library("ggpubr")
library("forcats")
setwd("/ifs/research-groups/botnar/proj025/analyses/test_trna_3")

file_1 <- read_tsv(file = "genome_statistics.dir/test.genomecov", col_names = FALSE)
file_2 <- read_tsv(file= "post_mapping_bams.dir/test_trna.idxstats", col_names =FALSE)
colnames(file_1) <- c("cluster", "position", "coverage")
colnames(file_2) <- c("cluster", "total_read_length", "total_coverage", "not_sure")

combined <- left_join(file_1,file_2, by = "cluster")
combined$percentage_coverage <- ifelse(combined$total_coverage > 0, (combined$coverage / combined$total_coverage )*100, 0)

sorted_combined <- arrange(combined, plyr::desc(total_coverage))
sorted_combined <- sorted_combined %>% dplyr::filter(cluster %in% head(unique(sorted_combined$cluster),50))
sorted_combined$cluster <- factor(sorted_combined$cluster, levels = unique(sorted_combined$cluster))

coverage_plot <- ggplot(data = sorted_combined, aes(x= position, y = fct_reorder(cluster, total_coverage, desc= TRUE))) + 
  geom_point(aes(colour = percentage_coverage)) +theme_bw() + labs(x= "Position in tRNA", y = "Cluster") +
  scale_colour_gradient(low= "white", high = "blue") + xlim(0, max(sorted_combined$total_read_length)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour= "black"))
print(coverage_plot)
# Need to incorporate DESeq2 and set for highest 50 differentially expressed
# Using total coverage and ranking by read depth of cluster overall