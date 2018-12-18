library("tidyverse")
library("ggpubr")
setwd("/ifs/research-groups/botnar/proj025/analyses/test_trna_3")

file_1 <- read_tsv(file = "genome_statistics.dir/test.genomecov", col_names = FALSE)
file_2 <- read_tsv(file= "post_mapping_bams.dir/test_trna.idxstats", col_names =FALSE)
colnames(file_1) <- c("cluster", "position", "coverage")
colnames(file_2) <- c("cluster", "total_read_length", "total_coverage", "not_sure")

combined <- left_join(file_1,file_2, by = "cluster")
combined$percentage_coverage <- ifelse(combined$total_coverage > 0, (combined$coverage / combined$total_coverage )*100, 0)


## 
# Just using 50 unique clusters, need to find top 50 differentially expressed
small_combined <- combined %>% dplyr::filter(cluster %in% head(unique(combined$cluster),50))

subset_combined <- dplyr::select(small_combined, -c(not_sure)) 
coverage_plot <- ggplot(data = subset_combined, aes(x= position, y = cluster)) + 
  geom_point(aes(colour = percentage_coverage)) +theme_bw() + labs(x= "Position in tRNA", y = "Cluster") + scale_colour_gradient(low= "white", high = "blue")
print(coverage_plot)

# Need to incorporate DESeq2 and set for highest 50 differentially expressed not just top 50 listed