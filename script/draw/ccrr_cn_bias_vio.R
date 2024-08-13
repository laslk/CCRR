library(jsonlite)
library(ggplot2)
library(reshape2)
library(optparse)


option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "Input file path", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL, 
              help = "Output path", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

consensus_cn_statistic <- read.delim(opt$input, header=FALSE, quote="'")
consensus_cn_statistic <- na.omit(consensus_cn_statistic)


consensus_cn_statistic$V6_range <- cut(consensus_cn_statistic$V6, 
                                       breaks = c(1, 1e4, 1e5, 1e6, 1e7, 1e8, Inf), 
                                       labels = c("(1,1e4)", "(1e4,1e5)", "(1e5,1e6)", "(1e6,1e7)", "(1e7,1e8)", ">1e8"),
                                       include.lowest = TRUE)

range_counts <- table(consensus_cn_statistic$V6_range)
range_counts_df <- as.data.frame(range_counts)
colnames(range_counts_df) <- c("Range", "Count")

sl <- ggplot(range_counts_df, aes(x=Range, y=Count, fill=Range)) + 
  geom_bar(stat="identity", width=0.75) + 
  scale_x_discrete(limits=c("(1,1e4)", "(1e4,1e5)", "(1e5,1e6)", "(1e6,1e7)", "(1e7,1e8)", ">1e8")) +
  scale_fill_brewer(palette="Set3") + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(hjust=1, colour="black"),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=14),
    panel.border = element_blank(),
    axis.line = element_line(colour="black")
  ) + 
  ylab("Count") + 
  xlab("Range") + 
  ggtitle("Distribution of Segment Length") +
  geom_text(aes(label=Count), vjust=-0.3) +
  coord_cartesian(ylim=c(0, max(range_counts_df$Count) * 1.1))+scale_y_continuous(expand = c(0,0))

ggsave(file.path(opt$output, "segment_count.pdf"), plot = sl, width = 6, height = 8, units = "in", dpi = 600)
parse_json_to_df <- function(json_str) {
  df <- as.data.frame(fromJSON(json_str))
  return(df)
}

all_colnames <- unique(unlist(lapply(consensus_cn_statistic$V7, function(x) names(fromJSON(x)))))
all_colnames_v8 <- unique(unlist(lapply(consensus_cn_statistic$V8, function(x) names(fromJSON(x)))))


combined_long <- data.frame()

for (range in levels(consensus_cn_statistic$V6_range)) {
  try({
    subset_consensus <- consensus_cn_statistic[consensus_cn_statistic$V6_range == range, ]
    bias <- do.call(rbind, lapply(subset_consensus$V7, function(x) {
      df <- parse_json_to_df(x)
      for (col in all_colnames) {
        if (!col %in% colnames(df)) {
          df[[col]] <- NA
        }
      }
      df <- df[all_colnames]
      return(df)
    }))
    volatility <- do.call(rbind, lapply(subset_consensus$V8, function(x) {
      df <- parse_json_to_df(x)
      for (col in all_colnames_v8) {
        if (!col %in% colnames(df)) {
          df[[col]] <- NA
        }
      }
      df <- df[all_colnames_v8]
      return(df)
    }))
    
    bias_long <- reshape2::melt(bias)
    volatility_long <- reshape2::melt(volatility)
    
    bias_long$type <- "Bias"
    volatility_long$type <- "Volatility"
    
    bias_long$range <- range
    volatility_long$range <- range
    
    combined_long <- rbind(combined_long, bias_long, volatility_long)
  }, silent=TRUE)
}

bav <- ggplot(combined_long, aes(x=variable, y=value, fill=type)) + 
  geom_boxplot(width=0.4, position=position_dodge(0.8), alpha=0.7, outlier.shape=NA) + 
  geom_violin(width=1.6, position=position_dodge(0.8), alpha=0.5) + 
  scale_fill_manual(values = c("Bias" = "red", "Volatility" = "blue")) + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(hjust=1, colour="black", family="Times", size=12),
    axis.text.y = element_text(family="Times", face="plain", size=12),
    axis.title.y = element_text(family="Times", face="plain", size=14),
    panel.border = element_blank(),
    axis.line = element_line(colour="black"),
    legend.position = "right"
  ) + 
  ylab("Value") + 
  xlab("Tool") + 
  ggtitle("Bias and Volatility for Each Tool by Range") + 
  facet_grid(. ~ range, scales = "fixed")

bav <- bav + theme(
  plot.margin = unit(c(1, 1, 1, 1), "cm"), 
  plot.background = element_rect(fill = "white", color = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_text(size=16, face="bold", hjust=0.5) 
)

ggsave(file.path(opt$output, "Bias_and_volatility_for_CN_all_ranges.pdf"), plot = bav, width = 24, height = 8, units = "in", dpi = 600)
