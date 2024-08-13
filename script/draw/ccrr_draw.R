library(circlize)
library(scales)
library(splines)
library(dplyr)
library(tidyr)
library(jsonlite)
library(ComplexHeatmap)
library(grid)
library(optparse)

SCRIPT_DIR <- Sys.getenv("SCRIPT_DIR")

option_list <- list(
  make_option(c("--region"), type="character", default=NULL, help="Path to the region file or region"),
  make_option(c("--shatterseek"), type="character", default=NULL),
  make_option(c("--ctlpscanner"), type="character", default=NULL),
  make_option(c("--sa"), type="character", default=NULL),
  make_option(c("--starfish"), type="character", default=NULL),
  make_option(c("--gGnome"), type="character", default=NULL),
  make_option(c("--aa"), type="character", default=NULL),
  make_option(c("--cn"), type="character", default=NULL),
  make_option(c("--out"), type="character", default=".", help="Output directory path"),
  make_option(c("--ref"), type="character", default=NULL),
  make_option(c("--sv"), type="character", default=NULL),
  make_option(c("--gapdegree"), type="integer", default=2, help="Gap degree [default=2]"),
  make_option(c("--startdegree"), type="integer", default=86, help="start degree [default=86]"),
  make_option(c("--cnsmoothspar"), type="double", default=0.3, help="CN smooth spar parameter [default=0.3]"),
  make_option(c("--singlechr_link_h"), type="double", default=-0.1, help="The height of lines representing structural variations on a single chromosome [default=-0.1]"),
  make_option(c("--trachr_link_h"), type="double", default=0.5, help="The height of lines representing TRA [default=0.5]"),
  make_option(c("--show_sector_region"), type="logical", default=FALSE, help="Display where the sector starts and ends [default=FALSE]"),
  make_option(c("--show_related_sv"), type="logical", default=TRUE, help="Display SV on complex event regions [default=FALSE]"),
  make_option(c("--show_unrelated_sv"), type="logical", default=FALSE, help="Display SV not on complex event regions [default=FALSE]"),
  make_option(c("--colorconfig "), type="character", default=file.path(SCRIPT_DIR,"draw/color.json"),help="Path to the color.json")
)
args <- parse_args(OptionParser(option_list=option_list))

##############################################################################
colorConfig <- fromJSON(args$colorconfig)
if (!is.null(args$ref)) {
  chromosome_data <- read.table(args$ref, quote="\"", comment.char="")
  colnames(chromosome_data) <- c('Genome','Length','Centromere')
} else {
   stop("ref error")
}

if (!is.null(args$cn)) {
  custom_cn <- read.delim(args$cn, header = FALSE)
  if (ncol(custom_cn) > 4) {
    custom_cn <- custom_cn[, -4]
  }
  colnames(custom_cn) <- c("chr", "start", "end", "copy_number")
} else {
   stop("cn error")
}

#############################################
if (!is.null(args$region)) {
  if (file.exists(args$region)) {
    region_data <- read.table(args$region, stringsAsFactors = FALSE)
    colnames(region_data) <- c('chr', 'start', 'end')
  } else {
    regions <- strsplit(args$region, "[\",\"]")[[1]]
    region_data <- do.call(rbind, lapply(regions, function(region) {
      parts <- strsplit(region, "[:-]")[[1]]
      if(length(parts) != 3) {
        stop("Region format error. Expected format 'chr:start-end'. Got: ", region)
      }
      tryCatch({
        chr <- parts[1]
        start <- as.integer(parts[2])
        end <- as.integer(parts[3])
        if (is.na(start) || is.na(end)) {
          stop("Start or end position is not a valid integer.")
        }
        return(data.frame(chr = chr, start = start, end = end, stringsAsFactors = FALSE))
      }, error = function(e) {
        stop("Error parsing region: ", region, ". Error message: ", e$message)
      })
    }))
    colnames(region_data) <- c('chr', 'start', 'end')
  }
} else {
  stop("Region argument is missing.")
}

if (!exists("region_data") || nrow(region_data) == 0) {
  stop("Region error: No valid regions specified.")
}

region_data$index <- row.names(region_data)
region_data$Centromere <- NA
region_data$Length <- NA
for (i in 1:nrow(region_data)) {
  row <- region_data[i, ]
  chrom_info <- chromosome_data[chromosome_data$Genome == row$chr, ]
  if (row$start >= row$end || row$start < 1 || row$end > chrom_info$Length) {
    stop("region error")
  }
  region_data$Length[i] <- row$end - row$start + 1
  if (!is.na(chrom_info$Centromere) && chrom_info$Centromere >= row$start && chrom_info$Centromere <= row$end) {
    region_data$Centromere[i] <- chrom_info$Centromere
  } else {
    region_data$Centromere[i] <- NA  
  }
}

for (chr in unique(region_data$chr)) {
  sub_data <- region_data[region_data$chr == chr, ]
  sub_data <- sub_data[order(sub_data$start), ]
  if (nrow(sub_data) > 1) {
    for (i in 1:(nrow(sub_data) - 1)) {
      if (sub_data$end[i] >= sub_data$start[i + 1]) {
        stop("region error")
      }
    }
  }
}

total_length <- sum(region_data$Length, na.rm = TRUE)
collected_regions <- data.frame(index = character(), start = numeric(), end = numeric(), stringsAsFactors = FALSE)
deduplicated_regions  <- data.frame(index = character(), start = numeric(), end = numeric(), stringsAsFactors = FALSE)



####################################################
try({
  png(filename = args$out, width = 8, height = 6, units = 'in', res = 1800)
  circos.clear()
  ##############################################################################
  #cn
  col_text <- colorConfig$main$text
  circos.par(canvas.xlim = c(-1.7, 1), canvas.ylim = c(-1, 1))
  circos.par("track.height" = 0.8, gap.degree = as.numeric(args$gapdegree), start.degree = as.numeric(args$startdegree), clock.wise = TRUE, cell.padding = c(0, 0, 0, 0))
  circos.initialize(factors = region_data$index, xlim = cbind(rep(0, nrow(region_data)), region_data$end-region_data$start))
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
      current_data <- region_data[region_data$index == CELL_META$sector.index, ]
      Genome = current_data$chr
      start = current_data$start
      xlim = CELL_META$xlim
      ylim = CELL_META$ylim
      if (args$show_sector_region){
        circos.text(mean(xlim), mean(ylim)+2.5, paste(Genome," : ",start,"-",current_data$end, sep = ""), cex = 0.5, col = col_text,facing="bending.outside",niceFacing = TRUE)
      } else {
        circos.text(mean(xlim), mean(ylim)+2.5,Genome, cex = 0.5, col = col_text,niceFacing = TRUE)
      }
      centromere_pos <- current_data$Centromere
      if (!is.na(centromere_pos)) {
          circos.rect(centromere_pos - 2*10^round(log10(total_length / 1000)) - start, ylim[1], centromere_pos + 2*10^round(log10(total_length / 1000)) - start, ylim[2], col = colorConfig$main$centromere, border = NA)
      }
  }, bg.col = colorConfig$main$background, bg.border = FALSE, track.height = 0.03)
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
      sector_data <- region_data[region_data$index == CELL_META$sector.index, ]
      magnitude <- 10^round(log10(total_length))
      step_choices <- c(1,2,3,4,5) * magnitude / 100
      steps_filtered <- step_choices[step_choices <= total_length/40]
      if (length(steps_filtered) > 0) {
          step <- max(steps_filtered)
      } else {
          step <- min(step_choices)
      }
      if ((sector_data$end - sector_data$start) >= step) {
          brk <- unique(c(sector_data$start, seq(ceiling(sector_data$start / step) * step, sector_data$end, by=step), sector_data$end))
          if (total_length>1000000){
            labels <- round(brk/1000000,2)
          }
          else {
            labels <- round(brk/1000,2)
          }

          if (sector_data$start %% step != 0) {
              labels[1] <- NA
          }
          if (sector_data$start == 1) {
              labels[1] <- 0
          }
          labels[length(labels)] <- NA
      } else {
          brk <- c(sector_data$start, sector_data$end)
          if (total_length>1000000){
            labels <- round(brk/1000000,2)
          }
          else {
            labels <- round(brk/1000,2)
          }
      }
      circos.axis(h="top",minor.ticks=0,major.at=brk-sector_data$start, labels=labels,labels.cex=0.3,col=col_text, labels.col=col_text, lwd=0.2, major.tick.length = 0.2)
  }, bg.border=FALSE)
  ###############################################################

  draw_cn <- data.frame(index=character(), start=integer(), end=integer(), copy_number=numeric(), stringsAsFactors=FALSE)
  for (i in 1:nrow(custom_cn)) {
    for (j in 1:nrow(region_data)) {
      if (custom_cn[i,]$chr == region_data[j,]$chr) {
        region_data[j,]$start
        if (custom_cn[i,]$start <= region_data[j,]$end && custom_cn[i,]$end >= region_data[j,]$start) {
          intersect_start <- max(custom_cn[i,]$start, region_data[j,]$start)
          intersect_end <- min(custom_cn[i,]$end, region_data[j,]$end)
          draw_cn <- rbind(draw_cn, data.frame(
            index = region_data[j,]$index,
            start = intersect_start-region_data[j,]$start,
            end = intersect_end-region_data[j,]$start,
            copy_number = custom_cn[i,]$copy_number
          ))
        }
      }
    }
  }
  min_cn <- min(draw_cn$copy_number)
  max_cn <- max(draw_cn$copy_number)

  if (min_cn >= 1.95) {
    min_cn <- 1.5
  }
  if (max_cn <= 2.05) {
    max_cn <- 2.5
  }

  draw_cn$adjusted_copy_number = log10((draw_cn$copy_number+1)/(max_cn*2))
  min_adjusted_cn<- min(draw_cn$adjusted_copy_number)
  max_adjusted_cn <- max(draw_cn$adjusted_copy_number)


  color_mapping_high <- col_numeric(palette = c(colorConfig$copynumber$high_cn_start, colorConfig$copynumber$high_cn_end), domain = c(log10(3/(max_cn*2)), max_adjusted_cn))
  color_mapping_low <- col_numeric(palette = c(colorConfig$copynumber$low_cn_start, colorConfig$copynumber$low_cn_end), domain = c(min_adjusted_cn, log10(3/(max_cn*2))))
  draw_cn$color <- mapply(function(adjusted_cn) {
    if (adjusted_cn > log10(3/(max_cn*2))) {
      color_mapping_high(adjusted_cn)
    } else if (abs(adjusted_cn - log10(3/(max_cn*2))) < 1e-5) {
      colorConfig$main$background
    } else {
      color_mapping_low(adjusted_cn)
    }
  }, adjusted_cn = draw_cn$adjusted_copy_number)

  circos.genomicTrack(draw_cn, panel.fun = function(region, value, ...) {
    for (i in 1:nrow(region)) {
      circos.genomicRect(region[i, ], ytop = log10((max_cn+1)/(max_cn*2))+0.1*abs(log10((max_cn+1)/(max_cn*2))), ybottom = log10((min_cn+1)/(max_cn*2))-0.1*abs(log10((min_cn+1)/(max_cn*2))),
                        col = value[i, "color"], border = NA)
    }
    adjusted_copy_number = log10((value[["copy_number"]]+1)/(max_cn*2))
    points = data.frame(x = numeric(0), y = numeric(0))
    
    for(i in 1:nrow(region)) {
      length = region$end[i] - region$start[i]
      if(length < 10^round(log10(total_length / 1000))) {
        mid_point = (region$start[i] + region$end[i]) / 2
        points = rbind(points, data.frame(x = mid_point, y = adjusted_copy_number[i]))
      } else {
        for(j in seq(from = region$start[i], to = region$end[i], by = 10^round(log10(total_length / 1000)))) {
          points = rbind(points, data.frame(x = j, y = adjusted_copy_number[i]))
        }
      }
    }
    if (nrow(region) > 0) {
      points = rbind(
        data.frame(x = region$start[1], y = adjusted_copy_number[1]),
        points,
        data.frame(x = region$end[nrow(region)], y = adjusted_copy_number[nrow(region)])
      )
    }
    points = unique(points)
    if (nrow(points) > 3) {
      spline_data = smooth.spline(points$x, points$y,cv=TRUE,spar=args$cnsmoothspar)
      for (i in 1:(length(spline_data$x) - 1)) {
        circos.lines(
          x = spline_data$x[c(i, i + 1)],
          y = spline_data$y[c(i, i + 1)],
          sector.index = CELL_META$sector.index,
          col = colorConfig$copynumber$cn_line,
          lwd = 0.3,
        )
      }
    } else {
      for (i in 1:(length(points$x) - 1)) {
        circos.lines(
          x = points$x[c(i, i + 1)],
          y = points$y[c(i, i + 1)],
          sector.index = CELL_META$sector.index,
          col = colorConfig$copynumber$cn_line,
          lwd = 0.3,
        )
      }
    }
    circos.lines(CELL_META$cell.xlim, c(log10(3/(max_cn*2)), log10(3/(max_cn*2))), lty = 3, lwd = 0.7, col = colorConfig$copynumber$cn_is_2_line)
  }, track.height = 0.1, ylim = c(log10((min_cn+1)/(max_cn*2))-0.1*abs(log10((min_cn+1)/(max_cn*2))), log10((max_cn+1)/(max_cn*2))+0.1*abs(log10((max_cn+1)/(max_cn*2)))),bg.border=NA)


  #############################################################################
  # shatterseek

  if (!is.null(args$shatterseek)) {
    shatterseek_summary <- read.delim(args$shatterseek)
    check_shatterseek <- function(row) {
      sum1_cols <- sapply(row[c("number_DEL", "number_DUP", "number_h2hINV", "number_t2tINV", "number_TRA", "clusterSize_including_TRA")],
                          function(x) ifelse(is.na(x), NA, as.numeric(x)))
      sum2_cols <- sapply(row[c("inter_number_DEL", "inter_number_h2hINV", "inter_number_t2tINV", "inter_number_DUP")],
                          function(x) ifelse(is.na(x), NA, as.numeric(x)))
      sum1 <- sum(sum1_cols, na.rm = TRUE)
      sum2 <- sum(sum2_cols, na.rm = TRUE)
      max_osc_CN_seg <- ifelse(is.na(row["max_number_oscillating_CN_segments_2_states"]), NA, as.numeric(row["max_number_oscillating_CN_segments_2_states"]))
      if (!is.na(sum1) && !is.na(max_osc_CN_seg)) {
        if (sum1 >= 6 && max_osc_CN_seg >= 7) {
          return("high")
        } else if (sum1 >= 3 && sum2 >= 4 && max_osc_CN_seg >= 7) {
          return("high")
        } else if (sum1 >= 6 && max_osc_CN_seg >= 4) {
          return("low")
        }
      }
      return(NA)
    }

    shatterseek_summary$classification <- apply(shatterseek_summary, 1, check_shatterseek)
    if (nrow(shatterseek_summary)>0){
      shatterseek_summary$chrom <- paste("chr", shatterseek_summary$chrom, sep = "")
    }
    shatterseek_color_mapping <- function(classification) {
      if (classification == "high") {
        return(colorConfig$shatterseek$high_confidence)
      } else if (classification == "low") {
        return(colorConfig$shatterseek$low_confidence)
      } else {
        return(colorConfig$main$background)
      }
    }

    circos.track(ylim=c(0,1),panel.fun = function(x, y) {
      chrom <- region_data[region_data$index==CELL_META$sector.index,]$chr
      region_start <- region_data[region_data$index==CELL_META$sector.index,]$start
      region_end <- region_data[region_data$index==CELL_META$sector.index,]$end
      chr_data <- shatterseek_summary[shatterseek_summary$chrom == chrom, ]

      if (!is.na(chr_data$start[1]) && !is.na(chr_data$end[1]) && !is.na(chr_data$classification[1]) &&
        region_start <= chr_data$end[1] && region_end >= chr_data$start[1]) {
      intersection_start <- max(region_start, chr_data$start[1])
      intersection_end <- min(region_end, chr_data$end[1])
      circos.rect(0, 0, intersection_start-region_start, 1, col = colorConfig$main$background, border = NA)
      circos.rect(intersection_start-region_start, 0, intersection_end-region_start, 1, 
                  col = shatterseek_color_mapping(chr_data$classification[1]), border = NA)
      circos.rect(intersection_end-region_start, 0, region_end-region_start, 1, col = colorConfig$main$background, border = NA)
      collected_regions <<- rbind(collected_regions, data.frame(index = CELL_META$sector.index, 
                                                                start = intersection_start, 
                                                                end = intersection_end))
    } else {
      circos.rect(0, 0, region_end-region_start, 1, col = colorConfig$main$background, border = NA)
    }
    }, bg.border = NA,track.height=0.04)
  }


  #############################################################################
  #ctlpscanner
  if (!is.null(args$ctlpscanner)) {
    CTLPRegion <- read.csv(args$ctlpscanner, sep="")
    CTLPRegion <- CTLPRegion[, -c(1:3)]
    if (nrow(CTLPRegion)>0){
      CTLPRegion$Chrom <- paste("chr", CTLPRegion$Chrom, sep = "")
    }
    CTLPRegion[,6]=log10(CTLPRegion[,6])
    if (nrow(CTLPRegion)>0){
      CTLPcolor <- col_numeric(palette = c(colorConfig$ctlpscanner$start, colorConfig$ctlpscanner$end), domain = c(5, ceiling(max(CTLPRegion[,6], na.rm = TRUE))))
    }
    CTLPRegion$color <- mapply(function(LR) {
      if (LR > 5) {
        CTLPcolor(LR)
      } else {
        colorConfig$main$background
      }
    }, LR = CTLPRegion[,6])
    for(chrom in c(paste0("chr", 1:22), "chrX", "chrY")) {
      if(!chrom %in% CTLPRegion$Chrom) {
        new_row <- setNames(as.list(rep(NA, ncol(CTLPRegion))), names(CTLPRegion))
        new_row$Chrom <-chrom
        CTLPRegion <- rbind(CTLPRegion, new_row)
      }
    }
    circos.track(ylim=c(0,1),panel.fun = function(x, y) {
      chrom <- region_data[region_data$index==CELL_META$sector.index,]$chr
      region_start <- region_data[region_data$index==CELL_META$sector.index,]$start
      region_end <- region_data[region_data$index==CELL_META$sector.index,]$end
      chr_data <- CTLPRegion[CTLPRegion$Chrom == chrom, ]
      if (!is.na(chr_data$Start[1]) && !is.na(chr_data$End[1]) && !is.na(chr_data$color[1]) &&
        region_start <= chr_data$End[1] && region_end >= chr_data$Start[1]) {
          intersection_start <- max(region_start, chr_data$Start[1])
          intersection_end <- min(region_end, chr_data$End[1])
          circos.rect(0, 0, intersection_start-region_start, 1, col = colorConfig$main$background, border = NA)
          circos.rect(intersection_start-region_start, 0, intersection_end-region_start, 1, col = chr_data$color[1], border = NA)
          circos.rect(intersection_end-region_start, 0,region_end-region_start
                      ,1, col = colorConfig$main$background, border = NA)
          if (chr_data$color[1] != colorConfig$main$background) {
            collected_regions <<- rbind(collected_regions, data.frame(index = CELL_META$sector.index, start = intersection_start, end = intersection_end))
          }  
        } else {
          circos.rect(0, 0, region_end-region_start, 1, col = colorConfig$main$background, border = NA)
        }
    }, bg.border = NA,track.height=0.04)
  }
  ##############################################################################
  #SA
  if(!is.null(args$sa)){
    SA_amplicons <- read.delim(args$sa)
    SA_amplicons <- SA_amplicons[, -c(4:18)]
    circos.track(ylim=c(0,1), panel.fun = function(x, y) {
      chrom <- region_data[region_data$index==CELL_META$sector.index,]$chr
      region_start <- region_data[region_data$index==CELL_META$sector.index,]$start
      region_end <- region_data[region_data$index==CELL_META$sector.index,]$end
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = colorConfig$main$background, border = NA)
      if (nrow(SA_amplicons) > 0) {
          for(i in 1:nrow(SA_amplicons)) {
            if (!is.na(SA_amplicons$seqnames[i]) && SA_amplicons$seqnames[i] == chrom && region_start <= SA_amplicons$end[i] && region_end >= SA_amplicons$start[i]) {
              intersection_start <- max(region_start, SA_amplicons$start[i])
              intersection_end <- min(region_end, SA_amplicons$end[i])
              circos.rect(intersection_start-region_start, 0, intersection_end-region_start, 1, col = colorConfig$SA$amplicon, border = NA)
              collected_regions <<- rbind(collected_regions, data.frame(index = CELL_META$sector.index, start = intersection_start, end = intersection_end))
            }
          }
        }
    }, bg.border = NA,track.height=0.04)
  }


  ##############################################################################
  #starfish
  if(!is.null(args$starfish)){
    tumor_connected_CGR_event <- read.csv(args$starfish)
    tumor_connected_CGR_event<-tumor_connected_CGR_event[, -c(4:7)]
    if(nrow(tumor_connected_CGR_event) > 0) {
      tumor_connected_CGR_event$chr <- paste("chr", tumor_connected_CGR_event$chr, sep = "")
    }
    circos.track(ylim=c(0,1), panel.fun = function(x, y) {
      chrom <- region_data[region_data$index==CELL_META$sector.index,]$chr
      region_start <- region_data[region_data$index==CELL_META$sector.index,]$start
      region_end <- region_data[region_data$index==CELL_META$sector.index,]$end
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = colorConfig$main$background, border = NA)
      if(nrow(tumor_connected_CGR_event) > 0) {
          for(i in 1:nrow(tumor_connected_CGR_event)) {
            if (!is.na(tumor_connected_CGR_event$chr[i]) && tumor_connected_CGR_event$chr[i] == chrom && region_start <= tumor_connected_CGR_event$end[i] && region_end >= tumor_connected_CGR_event$start[i]) {
              intersection_start <- max(region_start, tumor_connected_CGR_event$start[i])
              intersection_end <- min(region_end, tumor_connected_CGR_event$end[i])
              circos.rect(intersection_start-region_start, 0, intersection_end-region_start, 1, col = colorConfig$starfish$CGR, border = NA)
              collected_regions <<- rbind(collected_regions, data.frame(index = CELL_META$sector.index, start = intersection_start, end = intersection_end))
            }
          }
        }
    }, bg.border = NA,track.height=0.04)
  }

  ##############################################################################
  #gGnome
  if(!is.null(args$gGnome)){
    gGnome_event_footprints <- read.delim(args$gGnome)
    if(nrow(gGnome_event_footprints)>0){
      gGnome_event_footprints$Chromosome<- paste("chr", gGnome_event_footprints$Chromosome, sep = "")
    }
    unique_event <- unique(gGnome_event_footprints$Description)
    gGnome_colors <- colorConfig$gGnome$gGnome_colors
    while(length(gGnome_colors) < length(unique_event)) {
      gGnome_colors <- c(gGnome_colors, gGnome_colors)
    }
    gGnome_color_map <- setNames(gGnome_colors[1:length(unique_event)], unique_event)
    gGnome_event_footprints <- gGnome_event_footprints %>%
      mutate(color = gGnome_color_map[Description])
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
      Chromosome <- region_data[region_data$index==CELL_META$sector.index,]$chr
      region_start <- region_data[region_data$index==CELL_META$sector.index,]$start
      region_end <- region_data[region_data$index==CELL_META$sector.index,]$end
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = colorConfig$main$background, border = NA)

      if(nrow(gGnome_event_footprints)>0){
        for(i in 1:nrow(gGnome_event_footprints)) {
          if(!is.na(gGnome_event_footprints$Chromosome[i]) && gGnome_event_footprints$Chromosome[i] == Chromosome && region_start <= gGnome_event_footprints$End[i] && region_end >= gGnome_event_footprints$Start[i]) {
            intersection_start <- max(region_start, gGnome_event_footprints$Start[i])
            intersection_end <- min(region_end, gGnome_event_footprints$End[i])
            circos.rect(intersection_start-region_start-10^round(log10(total_length / 10000)), 0, intersection_end-region_start+10^round(log10(total_length / 10000)), 1, col = gGnome_event_footprints$color[i], border = NA)
            collected_regions <<- rbind(collected_regions, data.frame(index = CELL_META$sector.index, start = intersection_start, end = intersection_end))
          }
        }
      }
    }, bg.border = NA, track.height = 0.04)
  }
  # ##############################################################################
  #AA
  if(!is.null(args$aa)){
    AA_data <- fromJSON(args$aa)
    classification <- AA_data$Classification
    location <- AA_data$Location
    location_list <- lapply(location, function(x) {
      if (!is.na(x) && x != "[]" && x != "") {
        x <- gsub("'", "\"", x)
        return(fromJSON(x))
      } else {
        return(NA)
      }
    })
    AA_df <- data.frame(
      type = rep(classification, lengths(location_list)),
      location = unlist(location_list, use.names = FALSE)
    ) %>%
      mutate(location = ifelse(is.na(location), NA, location)) %>%
      separate(location, into = c("chr", "range"), sep = ":", fill = "right", extra = "drop") %>%
      separate(range, into = c("start", "end"), sep = "-", fill = "right", extra = "drop")
    AA_df <- AA_df %>%
      mutate(color = case_when(
        type == "ecDNA" ~ colorConfig$AA$ecDNA,
        type == "Linear" ~ colorConfig$AA$Linear,
        type == "BFB" ~ colorConfig$AA$BFB,
      )) %>%
      filter(complete.cases(.))

    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
      chrom <- region_data[region_data$index==CELL_META$sector.index,]$chr
      region_start <- region_data[region_data$index==CELL_META$sector.index,]$start
      region_end <- region_data[region_data$index==CELL_META$sector.index,]$end
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = colorConfig$main$background, border = NA)
      if(nrow(AA_df)>0){
        for(i in 1:nrow(AA_df)) {
          if(!is.na(AA_df$chr[i]) && AA_df$chr[i] == chrom&& region_start <= as.numeric(AA_df$end[i]) && region_end >= as.numeric(AA_df$start[i])) {
            intersection_start <- max(region_start, as.numeric(AA_df$start[i]))
            intersection_end <- min(region_end, as.numeric(AA_df$end[i]))
            circos.rect(intersection_start-region_start-10^round(log10(total_length / 10000)), 0, intersection_end-region_start+10^round(log10(total_length / 10000)), 1, col = AA_df$color[i], border = NA)
            collected_regions <<- rbind(collected_regions, data.frame(index = CELL_META$sector.index, start = intersection_start, end = intersection_end))
          }
        }
      }
    }, bg.border = NA, track.height = 0.04)
  }

  ##############################################################################
  if (nrow(collected_regions)>0){
    collected_regions <- collected_regions[order(collected_regions$index, collected_regions$start), ]
    deduplicated_regions <- data.frame(index = character(), start = numeric(), end = numeric())
    for(index in unique(collected_regions$index)) {
      index_regions <- collected_regions[collected_regions$index == index, ]
      i <- 1
      while(i <= nrow(index_regions)) {
        current_start <- index_regions$start[i]
        current_end <- index_regions$end[i]
        while(i < nrow(index_regions) && index_regions$start[i + 1] <= current_end + 1) {
          current_end <- max(current_end, index_regions$end[i + 1])
          i <- i + 1
        }
        deduplicated_regions <- rbind(deduplicated_regions, data.frame(index = index, start = current_start, end = current_end))
        i <- i + 1
      }
    }

    circos.track(ylim=c(0,1), panel.fun = function(x, y) {
      index <- CELL_META$sector.index
      circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = colorConfig$main$background, border = NA)
      region_start <- region_data[region_data$index==index,]$start
      region_end <- region_data[region_data$index==index,]$end
      if (nrow(deduplicated_regions)>0){
        for(i in 1:nrow(deduplicated_regions)) {
          if(deduplicated_regions$index[i] == index) {
            circos.rect(deduplicated_regions$start[i]-region_start, 0, deduplicated_regions$end[i]-region_start, 1, col = colorConfig$sv$deduplicated_region, border = NA)
          }
        }
      }
    }, bg.border = NA,track.height=0.01,track.margin = c(-0.01,0.12))
  }


  ##############################################################################
  #sv
  if(!is.null(args$sv)){
    custom_sv <- read.csv(args$sv, row.names=NULL, sep="\t")
    keep_row <- rep(TRUE, nrow(custom_sv))
    custom_sv$index1 <- NA
    custom_sv$index2 <- NA
    custom_sv$pos1_trans <- NA
    custom_sv$pos2_trans <- NA
    custom_sv$in_chrom <- FALSE
    custom_sv$in_region <- FALSE
    if (nrow(custom_sv)>0){
      for (i in 1:nrow(custom_sv)) {
        if (nrow(deduplicated_regions)>0){
          for (j in 1:nrow(deduplicated_regions)) {
            chrom1 <- region_data[region_data$index == deduplicated_regions$index[j],]$chr
            if (custom_sv$chrom1[i] == chrom1 &&
                custom_sv$pos1[i] >= deduplicated_regions$start[j] &&
                custom_sv$pos1[i] <= deduplicated_regions$end[j]) {
                  for (k in 1:nrow(region_data)) {
                    chrom2 <- region_data$chr[k]
                    if (custom_sv$chrom2[i] == chrom2 &&
                        custom_sv$pos2[i] >= region_data$start[k] &&
                        custom_sv$pos2[i] <= region_data$end[k]) {
                          custom_sv$index1[i] <- deduplicated_regions$index[j]
                          custom_sv$index2[i] <- region_data$index[k]
                          custom_sv$pos1_trans[i] <- custom_sv$pos1[i]-region_data[region_data$index == deduplicated_regions$index[j],]$start
                          custom_sv$pos2_trans[i] <- custom_sv$pos2[i]-region_data$start[k]
                          custom_sv$in_region[i] = TRUE
                          break
                    }
                  }
            }
          }
          if (!custom_sv$in_region[i]) {
            for (j in 1:nrow(deduplicated_regions)) {
              chrom2 <- region_data[region_data$index == deduplicated_regions$index[j],]$chr
              if (custom_sv$chrom2[i] == chrom2 &&
                  custom_sv$pos2[i] >= deduplicated_regions$start[j] &&
                  custom_sv$pos2[i] <= deduplicated_regions$end[j]) {
                    for (k in 1:nrow(region_data)) {
                      chrom1 <- region_data$chr[k]
                      if (custom_sv$chrom1[i] == chrom1 &&
                          custom_sv$pos1[i] >= region_data$start[k] &&
                          custom_sv$pos1[i] <= region_data$end[k]) {
                            custom_sv$index1[i] <- region_data$index[k]
                            custom_sv$index2[i] <- deduplicated_regions$index[j]
                            custom_sv$pos1_trans[i] <- custom_sv$pos1[i]-region_data$start[k]
                            custom_sv$pos2_trans[i] <- custom_sv$pos2[i]-region_data[region_data$index == deduplicated_regions$index[j],]$start
                            custom_sv$in_region[i] = TRUE
                            break
                      }
                    }
              }
            }
          }
        }
      }
      for (i in 1:nrow(custom_sv)) {
        if (nrow(region_data)>0){
          for (j in 1:nrow(region_data)) {
            chrom1 <- region_data$chr[j]
            if (custom_sv$chrom1[i] == chrom1  &&
                custom_sv$pos1[i] >= region_data$start[j] &&
                custom_sv$pos1[i] <= region_data$end[j]) {
                  for (k in 1:nrow(region_data)) {
                    chrom2 <- region_data$chr[k]
                    if (custom_sv$chrom2[i] == chrom2 &&
                        custom_sv$pos2[i] >= region_data$start[k] &&
                        custom_sv$pos2[i] <= region_data$end[k]) {
                          custom_sv$index1[i] <- region_data$index[j]
                          custom_sv$index2[i] <- region_data$index[k]
                          custom_sv$pos1_trans[i] <- custom_sv$pos1[i]-region_data$start[j]
                          custom_sv$pos2_trans[i] <- custom_sv$pos2[i]-region_data$start[k]
                          custom_sv$in_chrom[i] = TRUE
                          break
                    }
                  }
            }
          }
        }
      }
    }  


    inchrom_custom_sv <- custom_sv[custom_sv$in_chrom, ]    
    inregion_custom_sv <- inchrom_custom_sv[inchrom_custom_sv$in_region, ]
    notinregion_custom_sv <- inchrom_custom_sv[!inchrom_custom_sv$in_region, ]
    if (args$show_related_sv) {
      if (nrow(inregion_custom_sv)>0){
        for(i in 1:nrow(inregion_custom_sv)) {
          row <- inregion_custom_sv[i, ]
          if (row["svtype"]=='DEL'){
            circos.link(row["index1"], row["pos1_trans"], row["index2"], row["pos2_trans"],lwd=0.4,col = colorConfig$sv$line_del,h=args$singlechr_link_h)
          } else if (row["svtype"]=='DUP'){
            circos.link(row["index1"], row["pos1_trans"], row["index2"], row["pos2_trans"],lwd=0.4,col = colorConfig$sv$line_dup,h=args$singlechr_link_h)
          } else if (row["svtype"]=='h2hINV' | row["svtype"]=='t2tINV'){
            circos.link(row["index1"], row["pos1_trans"], row["index2"], row["pos2_trans"],lwd=0.4,col = colorConfig$sv$line_inv,h=args$singlechr_link_h)
          } else if (row["svtype"]=='TRA'){
            circos.link(row["index1"], row["pos1_trans"], row["index2"], row["pos2_trans"],lwd=0.6,col = colorConfig$sv$line_tra,h=args$trachr_link_h)
          }  
        }
      }
    }
    if (args$show_unrelated_sv) {
      if (nrow(notinregion_custom_sv)>0){
        for(i in 1:nrow(notinregion_custom_sv)) {
          row <- notinregion_custom_sv[i, ]
          if (row["svtype"]=='DEL'){
            circos.link(row["index1"], row["pos1_trans"], row["index2"], row["pos2_trans"],lwd=0.4,col = colorConfig$sv$nline_del,h=args$singlechr_link_h)
          } else if (row["svtype"]=='DUP'){
            circos.link(row["index1"], row["pos1_trans"], row["index2"], row["pos2_trans"],lwd=0.4,col = colorConfig$sv$nline_dup,h=args$singlechr_link_h)
          } else if (row["svtype"]=='h2hINV' | row["svtype"]=='t2tINV'){
            circos.link(row["index1"], row["pos1_trans"], row["index2"], row["pos2_trans"],lwd=0.4,col = colorConfig$sv$nline_inv,h=args$singlechr_link_h)
          } else if (row["svtype"]=='TRA'){
            circos.link(row["index1"], row["pos1_trans"], row["index2"], row["pos2_trans"],lwd=0.6,col = colorConfig$sv$nline_tra,h=args$trachr_link_h)
          }  
        }
      }
    }
  }

  ################################################

  y_start = 134

  ###################################################
  trackcount=1
  grid.text("Track1:CN", x = unit(6, "mm"), y = unit(y_start, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 10))
  lgd_cn_low = Legend(
    at = round(c(min_cn, 2),1), 
    col_fun = colorRamp2(c(min_cn, 2), c(colorConfig$copynumber$low_cn_start, colorConfig$copynumber$low_cn_end)),
    title_position = "topleft", labels_gp = gpar(fontsize = 4),
    grid_height = unit(3, "mm"),
    grid_width = unit(5, "mm"),
    direction = 'horizontal'
  )
  lgd_cn_high = Legend(
    at = round(c(2, max_cn),1),  
    col_fun = colorRamp2(c(2, max_cn), c(colorConfig$copynumber$high_cn_start, colorConfig$copynumber$high_cn_end)),
    title_position = "topleft", labels_gp = gpar(fontsize = 4),
    grid_height = unit(3, "mm"),
    grid_width = unit(5, "mm"),
    direction = 'horizontal'
  )

  lgd_cn_combined = packLegend(lgd_cn_low, lgd_cn_high, direction = "horizontal", column_gap = unit(1, "mm"))

  draw(lgd_cn_combined, x = unit(23, "mm"), y = unit(y_start-7, "mm"), just = c("left", "bottom"))
  grid.text("Gain or loss:", x = unit(8, "mm"), y = unit(y_start-4, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 7))

  lgd_cn_lines = Legend(at = c("cn"), type = "lines",
                        legend_gp = gpar(col = 'black', lwd = 2),labels_gp = gpar(fontsize = 8))
  lgd_cn_base_lines = Legend(at = c("cn=2"), type = "lines",
                            legend_gp = gpar(lty = 3, lwd = 1.5, col = colorConfig$copynumber$cn_is_2_line),labels_gp = gpar(fontsize = 8),)
  lgd_cn_line_combined = packLegend(lgd_cn_lines, lgd_cn_base_lines, direction = "horizontal", column_gap = unit(2, "mm"))
  draw(lgd_cn_line_combined, x = unit(23, "mm"), y = unit(y_start-12, "mm"), just = c("left", "bottom"))
  grid.text("CN line:", x = unit(8, "mm"), y = unit(y_start-11, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 7))
  y_start=y_start-12
  ################################################
  if (!is.null(args$shatterseek)){
    trackcount=trackcount+1
    grid.text(paste("Track",trackcount,": Shatterseek",sep = ""), x = unit(6, "mm"), y = unit(y_start-6, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 10))
    lgd_Shatterseek_high_points = Legend(at = c("High confidence"), type = "grid",
                                        legend_gp = gpar(fill = colorConfig$shatterseek$high_confidence, col = NA),labels_gp = gpar(fontsize = 8))
    lgd_Shatterseek_low_points = Legend(at = c( "Low confidence"), type = "grid",
                                        legend_gp = gpar(fill = colorConfig$shatterseek$low_confidence, col = NA),labels_gp = gpar(fontsize = 8))
    lgd_Shatterseek_points_combined = packLegend(lgd_Shatterseek_high_points, lgd_Shatterseek_low_points, direction = "vertical", row_gap = unit(1, "mm"))
    draw(lgd_Shatterseek_points_combined, x = unit(8, "mm"), y = unit(y_start-17, "mm"), just = c("left", "bottom"))
    y_start=y_start-17
  }

  ################################################
  if (!is.null(args$ctlpscanner)){
    trackcount=trackcount+1
    grid.text(paste("Track",trackcount,": CTLPScanner",sep = ""), x = unit(6, "mm"), y = unit(y_start-6, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 10))
    lrmax=ceiling(max(CTLPRegion[,6], na.rm = TRUE))
    if(is.infinite(lrmax)) { 
        lrmax <- 10 
    }
    lgd_ctlpscanner = Legend(
      at = c(5,(5+(5+lrmax)/2)/2,(5+lrmax)/2,(lrmax+(5+lrmax)/2)/2 ,lrmax), 
      col_fun = colorRamp2(c(5, lrmax), c(colorConfig$ctlpscanner$start, colorConfig$ctlpscanner$end)), 
      title_position = "topleft",labels_gp = gpar(fontsize = 4),
      grid_height = unit(3, "mm"),
      grid_width = unit(8, "mm"),
      direction = 'horizontal'
    )
    draw(lgd_ctlpscanner, x = unit(21, "mm"), y = unit(y_start-13, "mm"), just = c("left", "bottom"))
    grid.text("lg(LR)(>5):", x = unit(8, "mm"), y = unit(y_start-11, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 7))
    y_start = y_start-13
  }
  ################################################
  if (!is.null(args$sa)){
    trackcount=trackcount+1
    grid.text(paste("Track",trackcount,": SeismicAmplification",sep = ""), x = unit(6, "mm"), y = unit(y_start-6, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 10))
    lgd_sa_points = Legend(at = c("SA region"), type = "grid",
                        legend_gp = gpar(fill = colorConfig$SA$amplicon, col = NA),labels_gp = gpar(fontsize = 8))
    draw(lgd_sa_points, x = unit(8, "mm"), y = unit(y_start-12, "mm"), just = c("left", "bottom"))
    y_start = y_start-12
  }

  ################################################
  if (!is.null(args$starfish)){
    trackcount=trackcount+1
    grid.text(paste("Track",trackcount,": Starfish",sep = ""), x = unit(6, "mm"), y = unit(y_start-6, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 10))
    lgd_Starfish_points = Legend(at = c("CGR region"), type = "grid",
                                legend_gp = gpar(fill = colorConfig$starfish$CGR, col = NA),labels_gp = gpar(fontsize = 8))
    draw(lgd_Starfish_points, x = unit(8, "mm"), y = unit(y_start-12, "mm"), just = c("left", "bottom"))
    y_start = y_start-12
  }
  ################################################
  if (!is.null(args$gGnome)){
    trackcount=trackcount+1
    grid.text(paste("Track",trackcount,": gGnome",sep = ""), x = unit(6, "mm"), y = unit(y_start-6, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 10))
    lgd_gGnome_colors <- lapply(colorConfig$gGnome$gGnome_colors, function(color) {
      Legend(at = c(""), type = "grid", 
            legend_gp = gpar(fill = color, col = NA), 
            labels_gp = gpar(fontsize = 8))
    })
    lgd_gGnome_colors[[length(lgd_gGnome_colors) ]] <- Legend(at = c("complex events"), type = "grid",
                                                             legend_gp = gpar(fill = colorConfig$gGnome$gGnome_colors[length(lgd_gGnome_colors)], col = NA), 
                                                             labels_gp = gpar(fontsize = 8))
    lgd_gGnome_combined <- do.call(packLegend, c(lgd_gGnome_colors, list(direction = "horizontal", column_gap = unit(-0.5, "mm"))))                    
    draw(lgd_gGnome_combined, x = unit(8, "mm"), y = unit(y_start-12, "mm"), just = c("left", "bottom"))
    y_start = y_start-12
  }

  ################################################
  if (!is.null(args$aa)) {
    trackcount=trackcount+1
    grid.text(paste("Track",trackcount,": AmpliconArchitect",sep = ""), x = unit(6, "mm"), y = unit(y_start-6, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 10))
    lgd_AA_ecdna = Legend(at = c("ecDNA"), type = "grid",
                                        legend_gp = gpar(fill = colorConfig$AA$ecDNA, col = NA),labels_gp = gpar(fontsize = 8))
    lgd_AA_linear = Legend(at = c( "Linear"), type = "grid",
                                        legend_gp = gpar(fill = colorConfig$AA$Linear, col = NA),labels_gp = gpar(fontsize = 8))
    lgd_AA_bfb = Legend(at = c("BFB"), type = "grid",
                                        legend_gp = gpar(fill = colorConfig$AA$BFB, col = NA),labels_gp = gpar(fontsize = 8))
    lgd_AA_points_combined = packLegend(lgd_AA_ecdna, lgd_AA_linear,lgd_AA_bfb ,direction = "vertical", row_gap = unit(1, "mm"))
    draw(lgd_AA_points_combined, x = unit(8, "mm"), y = unit(y_start-22, "mm"), just = c("left", "bottom"))
    y_start = y_start-22
  }

  ################################################
  if (!is.null(args$sv)){
    trackcount=trackcount+1
    grid.text(paste("Track",trackcount,": SV",sep = ""), x = unit(6, "mm"), y = unit(y_start-6, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 10))
    if (args$show_related_sv && args$show_unrelated_sv) {
      lgd_sv_linel = Legend(at = c("", "","",""), type = "lines", 
                            legend_gp = gpar(col = c(colorConfig$sv$line_del,colorConfig$sv$line_dup,colorConfig$sv$line_inv,colorConfig$sv$line_tra),lwd = 2))
      lgd_sv_liner = Legend(at = c("DEL", "DUP","INV","TRA"), type = "lines", 
                            legend_gp = gpar(col = c(colorConfig$sv$nline_del,colorConfig$sv$nline_dup,colorConfig$sv$nline_inv,colorConfig$sv$nline_tra),lwd = 2))
      lgd_sv_combined = packLegend(lgd_sv_linel, lgd_sv_liner, direction = "horizontal", column_gap = unit(0, "mm"))      
      draw(lgd_sv_combined, x = unit(8, "mm"), y = unit(y_start-23, "mm"), just = c("left", "bottom"))
    }
    if (args$show_related_sv && !args$show_unrelated_sv) {
      lgd_sv_line = Legend(at = c("DEL", "DUP","INV","TRA"), type = "lines", 
                            legend_gp = gpar(col = c(colorConfig$sv$line_del,colorConfig$sv$line_dup,colorConfig$sv$line_inv,colorConfig$sv$line_tra),lwd = 2))      
      draw(lgd_sv_line, x = unit(8, "mm"), y = unit(y_start-23, "mm"), just = c("left", "bottom"))
    }
    if (!args$show_related_sv && args$show_unrelated_sv) {
      lgd_sv_line = Legend(at = c("DEL", "DUP","INV","TRA"), type = "lines", 
                            legend_gp = gpar(col = c(colorConfig$sv$nline_del,colorConfig$sv$nline_dup,colorConfig$sv$nline_inv,colorConfig$sv$nline_tra),lwd = 2))      
      draw(lgd_sv_line, x = unit(8, "mm"), y = unit(y_start-23, "mm"), just = c("left", "bottom"))
    }
  }
}, silent = FALSE)
dev.off()
