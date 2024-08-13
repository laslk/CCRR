args <- commandArgs(trailingOnly = TRUE)
cnv_file <- args[1]
sv_file <- args[2]
genome_version <- args[3]
output_path <- args[4]
TOOL_DIR <- args[5]

library(GenomicRanges)
source(file.path(TOOL_DIR,"/SA/seismic_amplification_detection/seismic_amplification_detection.R"))
if (dir.exists(output_path)) {
  unlink(output_path, recursive = TRUE)
}
dir.create(output_path)
setwd(output_path)
if (genome_version == "hg19") {
  bands_file <- file.path(TOOL_DIR,"/SA/seismic_amplification_detection/chromosome_bands_hg19.txt")
} else if (genome_version == "hg38") {
  bands_file <- file.path(TOOL_DIR,"/SA/seismic_amplification_detection/chromosome_bands_hg38.txt")
} else {
  stop("Invalid genome version. Please specify either 'hg19' or 'hg38'.")
}


cnv_data <- read.table(cnv_file, header = FALSE, sep = "\t")
colnames(cnv_data) <- c("chrom", "start", "end", "cn")
cnv_gr <- GRanges(seqname = cnv_data$chrom, ranges = IRanges(start = cnv_data$start, end = cnv_data$end), cn = cnv_data$cn)


bp_data <- read.table(sv_file, header = TRUE, sep = "\t")
bp <- bp_data[, c("chrom1", "pos1", "chrom2", "pos2")]
colnames(bp) <- c("chr1", "bp1", "chr2", "bp2")


bands <- read.table(bands_file, header = TRUE, sep = "\t")
bands_gr <- GRanges(seqname = bands$chrom, ranges = IRanges(start = bands$chromStart, end = bands$chromEnd), name = bands$name, gieStain = bands$gieStain)


sa <- detect_seismic_amplification(cnv = cnv_gr, sv = bp, chrBands = bands_gr)

write.table(sa$svs, paste0(output_path, "/SA_svs.csv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(sa$amplicons), paste0(output_path, "/SA_amplicons.csv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
