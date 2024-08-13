args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_file> <output_file>", call. = FALSE)
}

input_file <- args[1]
output_file <- args[2]

sv <- read.table(input_file, header = TRUE, sep = '\t')
sv <- sv[, !(names(sv) %in% c("from", "id", "svtype"))]
sv$chrom1 <- gsub("chr", "", sv$chrom1)
sv$chrom2 <- gsub("chr", "", sv$chrom2)
sv$pos1_copy <- sv$pos1
sv$pos2_copy <- sv$pos2
sv$sequence <- seq_len(nrow(sv))
sv <- sv[, c('chrom1', 'pos1', 'pos1_copy', 'chrom2', 'pos2', 'pos2_copy', 'sequence', 'level', 'strand1', 'strand2')]


while(nrow(sv) > 1000) {
  if(length(unique(sv$level)) <= 2) break 
  sv <- sv[order(-sv$level),] 
  lowest_level <- tail(unique(sv$level), 1)
  lowest_level_rows <- which(sv$level == lowest_level)
  if(length(lowest_level_rows) <= 20) { 
    sv <- sv[-lowest_level_rows,]
  } else {
    rows_to_remove <- sample(lowest_level_rows, length(lowest_level_rows) %/% 2)
    sv <- sv[-rows_to_remove,]
  }
}

write.table(sv, output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
