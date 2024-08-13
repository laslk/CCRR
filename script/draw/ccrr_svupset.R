library(ggplot2)
library(ComplexUpset)
library(optparse)


option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "Input file path", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL, 
              help = "Output file path", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

generate_upset_plot <- function(input_file, output_file) {
  data <- read.csv(input_file, sep="\t")
  parse_from_column <- function(value) {
    if (startsWith(value, "[")) {
      value <- gsub("\\[|\\]|'", "", value)
      return(strsplit(value, ",\\s*")[[1]])
    } else {
      return(trimws(value))
    }
  }
  data$from <- lapply(data$from, parse_from_column)
  tools <- c('gridss', 'manta', 'svaba', 'lumpy', 'soreca', 'delly')
  upset_data <- data.frame(matrix(0, nrow=nrow(data), ncol=length(tools)))
  colnames(upset_data) <- tools
  for (i in 1:nrow(data)) {
    for (tool in tools) {
      if (tool %in% data$from[[i]]) {
        upset_data[i, tool] <- 1
      }
    }
  }
  upset_data$svtype <- data$svtype
  p <- upset(
    upset_data,
    tools,
    base_annotations = list(
      'Intersection Counts' = intersection_size(
        mapping = aes(fill = svtype)
      ) + scale_fill_manual(values = c(
        'DEL' = '#E41A1C', 'DUP' = '#377EB8',
        'h2hINV' = '#4DAF4A', 't2tINV' = '#FF7F00', 'TRA' = '#213695'
      ), guide = "none")
    ),
    set_sizes = (
      upset_set_size(
        geom = geom_bar(
          aes(fill = svtype),))
      + geom_text(aes(label = after_stat(count)), width = 0.8, hjust = "inward", stat = 'count') 
      + scale_fill_manual(values = c(
        'DEL' = '#E41A1C', 'DUP' = '#377EB8',
        'h2hINV' = '#4DAF4A', 't2tINV' = '#FF7F00', 'TRA' = '#213695')
      )
    ), guides = 'over',wrap=TRUE
  )  + ggtitle('Consensus SVs of All Tools')
  ggsave(output_file, plot = p, device = "pdf",width = 15, height = 10, dpi = 800)
}

generate_upset_plot(opt$input, opt$output)