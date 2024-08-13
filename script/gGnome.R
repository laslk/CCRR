library(gGnome)
args <- commandArgs(trailingOnly = TRUE)
jabba_path <- args[1]
gg=gG(jabba=file.path(jabba_path,"jabba.rds"))
gg_events = events(gg, verbose = FALSE)
event_counts <- gg_events$meta$event[, table(type)]

setwd(file.path(jabba_path, "..", "gGnome"))
output_file <- "event_footprints.txt"
write.table(data.frame(Chromosome = character(), Start = integer(), End = integer(), Description = character()), file = output_file, row.names = FALSE, sep = "\t", quote = FALSE, col.names = TRUE)
plot_event <- function(event_type, i, gg_events) {
  event_folder <- paste0(event_type, "_events")
  dir.create(event_folder, showWarnings = FALSE)
  file_path <- paste0(event_folder, "/", event_type, i, ".pdf")
  pdf(file_path, width = 8, height = 6)
  if (event_type == 'inv') {
    selected_event <- eval(parse(text = paste0("gg_events$edges[simple=='INV", i, "']")))
    rang <- sum(width(ranges(selected_event$shadow)))*0.15+5e4
    plot(gg_events$gt, selected_event$shadow %>% streduce(rang))
    title(paste(event_type, " ", i))
    dev.off()
  } else if(event_type=='tra'){
    selected_event <- eval(parse(text = paste0("gg_events$edges[simple=='TRA", i, "']")))
    rang <- sum(width(ranges(selected_event$shadow)))*0.15+5e4
    plot(gg_events$gt, selected_event$shadow %>% streduce(rang))
    title(paste(event_type, " ", i))
    dev.off()
  } else if(event_type=='invdup'){
    selected_event <- eval(parse(text = paste0("gg_events$edges[simple=='INVDUP", i, "']")))
    rang <- sum(width(ranges(selected_event$shadow)))*0.15+5e4
    plot(gg_events$gt, selected_event$shadow %>% streduce(rang))
    title(paste(event_type, i))
    dev.off()
  } else {
    selected_event <- eval(parse(text = paste0("gg_events[",event_type, "==", i, "]")))
    rang <- width(ranges(selected_event$footprint))*0.15+5e4
    plot(gg_events$gt, selected_event$footprint +rang)
    title(paste(event_type,i))
    dev.off()
    if (!(event_type %in% c('del', 'dup'))) {
        to_write <- data.frame(
            Chromosome = as.character(seqnames(selected_event$footprint)),
            Start = start(ranges(selected_event$footprint)),
            End = end(ranges(selected_event$footprint)),
            Description = paste(event_type, i, sep="")
        )
    write.table(to_write, file = output_file, append = TRUE, row.names = FALSE, sep = "\t", quote = FALSE, col.names = FALSE)
    }
  }
}

for (event_type in names(event_counts)) {
    event_count=event_counts[event_type]
    for (i in 1:event_count) {
        try(plot_event(event_type, i, gg_events))
    }
}
