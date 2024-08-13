Packages <- c("ShatterSeeky", "GenomeInfoDb", "plyr", "data.table", "GenomicRanges", "IRanges", "MASS", "ggplot2", "grid", "gridExtra", "dplyr", "ConsensusClusterPlus", "factoextra", "gplots", "ggpubr", "reshape2", "cowplot", "scales", "patchwork", "Cairo", "ggforce")
lapply(Packages, library, character.only = TRUE)
library(Starfish)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("6 arguments are required: path to SV file,path to CNV file,genome version (hg19 or hg38),sample id and gender(Male or Female),id.")
  }
cnv_file <- args[2]
sv_file <- args[1]
genome_version <- args[3]
sample_id =args[4]
gender=args[5]
output_dir <- args[6]


if (dir.exists(output_dir)) {
  unlink(output_dir, recursive = TRUE)
}
dir.create(output_dir, recursive = TRUE)
setwd(output_dir)
sv=read.table(sv_file,header=TRUE,sep='\t')
sv <- sv[, !(names(sv) %in% c("from", "id", "level"))]
sv$sample=sample_id
cnv=read.table(cnv_file,header=FALSE,sep="\t")
names(cnv)=c("chromosome","start","end","total_cn","cnn")
cnv <- cnv[, !(names(cnv) %in% c("cnn"))]
cnv$sample=sample_id
sample=data.frame(sample=c(sample_id),gender=c(gender))
tryCatch({
    starfish_all(sv, cnv, sample, prefix=sample_id, genome_v=genome_version, cnv_factor="auto", arm_del_rm=TRUE, plot=TRUE, cmethod="class")
}, error=function(e){
    message("No CGR region is identified, starfish is done")
})

