library(dplyr)
library(ggplot2)

fastq_lines <- readLines("small_SRR1971253.fastq")
n_lines <- length(fastq_lines)
reads <- fastq_lines[seq(2, n_lines, 4)] 
qualities <- fastq_lines[seq(4, n_lines, 4)]  

stats <- data.frame(
  Total_Reads = length(reads),
  Avg_Length = mean(nchar(reads)),
  GC_Content = mean(sapply(reads, function(x) {
    gc <- sum(strsplit(x, "")[[1]] %in% c("G", "C"))
    gc / nchar(x)
  }))
)

phred_scores <- sapply(qualities, function(q) {
  scores <- as.integer(charToRaw(q)) - 33
  mean(scores, na.rm = TRUE) 
})
stats$Avg_Phred <- mean(phred_scores)

lengths <- nchar(reads)
len_df <- data.frame(Length = lengths)
p1 <- ggplot(len_df, aes(x = Length)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Distribution of Read Lengths", x = "Length", y = "Count") +
  theme_minimal()

qual_df <- data.frame(Mean_Phred = phred_scores)
p2 <- ggplot(qual_df, aes(x = Mean_Phred)) +
  geom_histogram(binwidth = 1, fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Mean Phred Scores", x = "Mean Phred Score", y = "Count") +
  theme_minimal()

ggsave("length_plot.png", p1, width = 8, height = 5)
ggsave("quality_plot.png", p2, width = 8, height = 5)
write.csv(stats, "read_stats_simple.csv", row.names = FALSE)
print(stats)