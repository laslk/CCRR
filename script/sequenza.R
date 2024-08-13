library(sequenza)
library(jsonlite)
args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])
id=args[2]
parallel=args[3]
chromosome_list <- paste0("chr", 1:22)
data <- sequenza.extract("small.seqz", verbose = FALSE,chromosome.list=chromosome_list,parallel=as.numeric(parallel))
CP <- sequenza.fit(data)
sequenza.results(sequenza.extract = data,
    cp.table = CP, sample.id = id,
    out.dir=".")
cint <- get.ci(CP)
cellularity <- cint$max.cellularity
ploidy <- cint$max.ploidy
jsonlite::write_json( list(cellularity = cellularity, ploidy = ploidy), "cellularity_ploidy.json")
