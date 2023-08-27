library(reshape2)
abd <- read.table("../data/Abundance_Stat.all.xls", sep = "\t", header = T, comment.char = "")
meta <- read.table("../data/meta.txt", sep = "\t", header = T, row.names = 1)

melt_abd <- melt(abd)
colnames(melt_abd)[1] <- "Kingdom"

melt_abd$group <- meta[match(melt_abd$variable, rownames(meta)), "Group"]
result <- aggregate(value ~ Kingdom + variable, melt_abd, sum)
dcast_result <- dcast(result, variable ~ Kingdom)
dcast_result$ratio <- dcast_result$Fungi / dcast_result$Bacteria

dcast_result$group <- meta[match(dcast_result$variable, rownames(meta)), "Group"]

fungi_bacteria_ratio <- dcast_result[, c("Bacteria", "Fungi", "ratio")]
rownames(fungi_bacteria_ratio) <- dcast_result$variable
fungi_bacteria_ratio <- t(fungi_bacteria_ratio)
write.table(fungi_bacteria_ratio, "../results/fungi_bacteria_ratio.txt", sep = "\t", quote = F, row.names = T, col.names = NA)

dcast_result2 <- dcast_result[, c("variable", "group", "ratio")]

x <- subset(dcast_result, group == "Clean_Prosthesis")$ratio
y <- subset(dcast_result, group == "Unclean_Prosthesis")$ratio

shapiro.test(x)
shapiro.test(y)

wilcoxon_result <- wilcox.test(x, y)
wilcoxon_result

x <- log10(x + 1e-60)
y <- log10(y + 1e-60)

max_length <- max(length(x), length(y))

# 补齐列表长度为最长长度
x <- c(x, rep(median(x), max_length - length(x)))
y <- c(y, rep(median(y), max_length - length(y)))

data <- data.frame(Clean_Prosthesis = x, Unclean_Prosthesis = y, stringsAsFactors = FALSE)
df <- melt(data)
# 绘制盒图
p <- ggplot(df, aes(x = variable, y = value, fill = variable)) + 
  geom_boxplot(fill = c("#3F5F9E", "#AC363E"), color = "black", outlier.shape = NA) +
  labs(x = "Group", y = "Ratio") +
  geom_text(x = 1.5, y = 1, label = paste("p =", round(wilcoxon_result$p.value, 3))) +
  theme_bw()

ggsave("../results/fungi_bacteria_ratio_boxplot.pdf", p, width = 4, height = 5)

