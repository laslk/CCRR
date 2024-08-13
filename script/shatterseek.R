library(ShatterSeek)
library(gridExtra)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("path to sv file, path to cnv file, genome version (hg19 or hg38), work dir")
}

sv_file <- args[1]
cnv_file <- args[2]
genome_version <- args[3]
output_dir <- args[4]

if (!file.exists(sv_file)) {
  stop("SV file does not exist: ", sv_file)
}
if (!file.exists(cnv_file)) {
  stop("CNV file does not exist: ", cnv_file)
}

if (dir.exists(output_dir)) {
  unlink(output_dir, recursive = TRUE)
}

dir.create(output_dir, recursive = TRUE)
setwd(output_dir)
sv_data=read.table(sv_file,header=TRUE,sep='\t')
sv_data$chrom1=gsub("chr","",sv_data$chrom1)
sv_data$chrom2=gsub("chr","",sv_data$chrom2)
equal_pos_indices = sv_data$pos1 == sv_data$pos2
sv_data = sv_data[!equal_pos_indices, ]
sv_data <- subset(sv_data, chrom1 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X") & chrom2 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"))

SV <- SVs(chrom1=as.character(sv_data$chrom1), pos1=as.numeric(sv_data$pos1), chrom2=as.character(sv_data$chrom2), pos2=as.numeric(sv_data$pos2), SVtype=as.character(sv_data$svtype), strand1=as.character(sv_data$strand1), strand2=as.character(sv_data$strand2))
cn_data=read.table(cnv_file,header=FALSE,sep='\t')
colnames(cn_data) <- c("chrom", "start", "end", "cn","ccn")
cn_data$chrom=gsub("chr","",cn_data$chrom)
cn_data <- subset(cn_data, chrom %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"))

CN <- CNVsegs(chrom=as.character(cn_data$chrom), start=cn_data$start, end=cn_data$end, total_cn=cn_data$cn)
chromothripsis <- shatterseek(SV.sample=SV, seg.sample=CN,genome=genome_version)
write.table(chromothripsis@chromSummary, file.path(output_dir, "chromothripsis_summary.csv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
chr_list <- c(1:22, "X")
for (chr in chr_list) {
  try({
    plots_chr = plot_chromothripsis(ShatterSeek_output = chromothripsis, chr = chr)
    plot_chr = arrangeGrob(grobs = plots_chr, nrow = 4, ncol = 1, heights = c(0.2, 0.6, 0.6, 0.6))
    ggsave(filename = file.path(output_dir, paste0("chromothripsis_chr", chr, ".png")), plot = plot_chr, width = 10, height = 8)  
  }, silent = TRUE)  
}
unlink(file.path(output_dir,"Rplots.pdf"))
