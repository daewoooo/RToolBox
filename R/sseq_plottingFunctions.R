
plotCompositeFileWCratio <- function(composite.files, region=NULL, chromosomes=paste0('chr', c(1:22, 'X')),  binsize = 100000, stepsize = 10000, min.reads = 0, alpha.cov = FALSE, annot.gr = NULL, highlight.gr = NULL, bsgenome = NULL, fa.index = NULL, plot.ideo = FALSE) {
  ## Get list of composite files
  #composite.files <- list.files(path = composite.files.path, pattern = "_th.2_minR10_back.02.RData", full.names = TRUE)
  ## Load BSgenome
  if (class(bsgenome) != 'BSgenome') {
    if (is.character(bsgenome)) {
      suppressPackageStartupMessages(library(bsgenome, character.only=TRUE))
      bsgenome <- as.object(bsgenome) # replacing string by object
    }
  } else if (is.null(bsgenome) & is.null(fa.index)) {
    stop("Parameter 'bsgenome' and 'fa.index' is not defined, aborting!!!")
  }
  
  ## Subselect annotation ranges
  if (!is.null(region) & class(region) == 'GRanges') {
    chromosomes <- as.character(unique(seqnames(region)))
    if (!is.null(annot.gr) & class(annot.gr) == 'GRanges') {
      annot.gr <- subsetByOverlaps(annot.gr, region)
    }
    if (!is.null(highlight.gr) & class(highlight.gr) == 'GRanges') {
      highlight.gr <- subsetByOverlaps(highlight.gr, region)
    }
  } 
  
  ## Create genomic bins
  if (class(bsgenome) == 'BSgenome') {
    binned.data <- primatR::makeBins(bsgenome = bsgenome, chromosomes = chromosomes, binsize = binsize, stepsize = stepsize)
  } else if (file.exists(fa.index)) {
    binned.data <- primatR::makeBins(fai = fa.index, chromosomes = chromosomes, binsize = binsize, stepsize = stepsize)
  }
  
  reads.grl <- GRangesList()
  for (i in seq_along(composite.files)) {
    composite.file <- composite.files[i]
    filename <- basename(composite.file)
    sample <- strsplit(filename, "_")[[1]][2]
    sample <- gsub(pattern = '\\.RData$', replacement = "", sample)
    #popul <- samples$superpop[grep(samples$name, pattern = sample.num)]
    #sample.col <- colors[[popul]]
    ## Load composite file reads
    directional.reads <- get(load(composite.file))
    ## Get reads from user defined region
    if (!is.null(region) & class(region) == 'GRanges') {
      directional.reads <- subsetByOverlaps(directional.reads, region)
      binned.data <- subsetByOverlaps(binned.data, region)
    }
    ## Get Watson and Crick read counts
    Watsonreads <- GenomicRanges::countOverlaps(binned.data, directional.reads[strand(directional.reads)=='-']) 
    Crickreads <- GenomicRanges::countOverlaps(binned.data, directional.reads[strand(directional.reads)=='+'])
    bothreads <- Watsonreads + Crickreads
    mcols(binned.data)$bothreads <- bothreads
    mcols(binned.data)$Watsonreads <- Watsonreads
    mcols(binned.data)$Crickreads <- Crickreads
    ## Get Watson and Crick ratios
    mcols(binned.data)$Watson.reads.frac <- mcols(binned.data)$Watsonreads / (mcols(binned.data)$Watsonreads + mcols(binned.data)$Crickreads)
    mcols(binned.data)$Crick.reads.frac <- mcols(binned.data)$Crickreads / (mcols(binned.data)$Watsonreads + mcols(binned.data)$Crickreads)
    binned.data$sample <- sample
    ## Set NaNs to zero
    mcols(binned.data)$Watson.reads.frac[is.nan(mcols(binned.data)$Watson.reads.frac)] <- 0
    mcols(binned.data)$Crick.reads.frac[is.nan(mcols(binned.data)$Crick.reads.frac)] <- 0
    ## Set bins with less than min.reads to NA
    mcols(binned.data)$Watson.reads.frac[mcols(binned.data)$bothreads < min.reads] <- NA
    mcols(binned.data)$Crick.reads.frac[mcols(binned.data)$bothreads < min.reads] <- NA
    ## Set alpha level based on read counts
    if (alpha.cov) {
      ## Get quantiles for specific values of the data distribution
      quantile.range <- quantile(binned.data$bothreads, probs = seq(0, 1, 0.2))
      ## prepare label text (use two adjacent values for range text)
      #label.text <- rollapply(round(quantile.range, 2), width = 2, by = 1, FUN = function(i) paste(i, collapse = " : "))
      ## Find category of predefined ranges defined by 'quantile.range' 
      categ <- findInterval(binned.data$bothreads, quantile.range, all.inside = TRUE)
      binned.data$cov.categ <- factor(categ, levels = 1:5)
    }
    #binned.region$superpop <- as.character(popul)
    reads.grl[[i]] <- binned.data
  }
  reads.gr <- unlist(reads.grl, use.names = FALSE)
  plt.df <- as.data.frame(reads.gr)
  plt.df$midpoint <- plt.df$start + ( (plt.df$end - plt.df$start) %/% 2 )
  ## Make sure samples are in the same order as the input composite files
  plt.df$sample <- factor(plt.df$sample, levels = unique(plt.df$sample))
  
  my.theme <- theme(legend.position="none",
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    strip.text.y.left = element_text(angle = 0))
  
  ## Construct ggplot
  if (alpha.cov) {
    plt <- ggplot(plt.df) + 
      geom_linerange(aes(ymin=1, ymax=Watson.reads.frac, x=midpoint, color='Crick', alpha=cov.categ), size=1) + 
      geom_linerange(aes(ymin=Watson.reads.frac, ymax=0, x=midpoint, color='Watson',alpha=cov.categ), size=1) +
      scale_alpha_manual(values = c('1' = 0.2, '2' = 0.4, '3'=0.6, '4'=0.8, '5'=1))
  } else {
    plt <- ggplot(plt.df) + 
      geom_linerange(aes(ymin=1, ymax=Watson.reads.frac, x=midpoint, color='Crick'), size=1) + 
      geom_linerange(aes(ymin=Watson.reads.frac, ymax=0, x=midpoint, color='Watson'), size=1)
  } 
  plt <- plt + geom_hline(yintercept = 0.5, color='white', linetype='dashed', size=0.25) +
    facet_grid(sample ~ ., switch='y') +
    scale_color_manual(values = c("paleturquoise4", "sandybrown")) +
    scale_x_continuous(expand = c(0,0), name = as.character(region)) +
    my.theme
  
  ## Highlight regions
  if (!is.null(highlight.gr) & class(highlight.gr) == 'GRanges') {
    highlight.gr <- subsetByOverlaps(highlight.gr, region)
    highlight.df <- as.data.frame(highlight.gr)
    plt <- plt + geom_rect(data=highlight.df, aes(xmin=start, xmax=end, ymin=0, ymax=Inf), fill='red', alpha=0.2)
  }
  ## Add annotation
  if (!is.null(annot.gr) & length(annot.gr) > 0 & class(annot.gr) == 'GRanges') {
    annot.df <- as.data.frame(annot.gr)
    #annot.df$ID <- 'annot'
    ## Set x axis range
    range.df <- as.data.frame(range(reads.gr))
    annot.plt <- ggplot() +
      geom_rect(data=range.df, aes(xmin=start, xmax=end, ymin=0, ymax=1), fill='white') +
      geom_rect(data=annot.df, aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=ID)) +
      scale_fill_manual(values = brewer.pal(n = 9, name = 'Set1')) +
      scale_x_continuous(expand = c(0,0), labels = comma) +
      facet_grid(ID ~ ., switch='y') +
      my.theme
    plt <- plot_grid(plt, annot.plt, ncol = 1, align = 'v', axis = 'lr', rel_heights = c(1, 0.15))
  }
  ## Add ideogram
  if (plot.ideo & !is.null(region) & class(region) == 'GRanges') {
    region.df <- as.data.frame(region)
    chr <- as.character(unique(seqnames(region)))
    if (is.character(chr) & length(chr) == 1) {
      hg38IdeogramCyto <- getIdeogram("hg38", cytobands = TRUE, subchr = chr)
      ideo.df <- as.data.frame(hg38IdeogramCyto)
      ideo.plt <- ggplot() + 
        geom_rect(ideo.df, aes(xmin=start, xmax=end, ymin=-0.1, ymax=0.1, fill=gieStain), color='black', show.legend = FALSE) + 
        geom_rect(region.df, aes(xmin=start, xmax=end, ymin=-0.1, ymax=0.1), fill='red', alpha=0.5) + 
        scale_fill_giemsa() +
        scale_x_continuous(labels = comma) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
      plt <- plot_grid(ideo.plt, plt, ncol = 1, align = 'v', axis = 'lr', rel_heights = c(0.1, 1))
    }
  }  
  return(plt) 
} 

plotBinnedReads <- function(reads.gr = NULL, highlight.gr = NULL, binsize = 100000, stepsize = 50000, title = NULL) {
  ## Set stepsize if not defined
  if (is.numeric(stepsize)) {
    if (stepsize < binsize) {
      stepsize == stepsize
    } else {
      stepsize == binsize
    }
  } else {
    stepsize == binsize
  }
  ## Make sure seqlevels reflect data
  reads.gr <- keepSeqlevels(reads.gr, value = unique(as.character(seqnames(reads.gr))), pruning.mode = 'coarse')
  
  ## Bin the genome
  gen.ranges <- range(reads.gr, ignore.strand=TRUE)
  bins <- GRangesList()
  for (i in seq_along(gen.ranges)) {
    range <- gen.ranges[i]
    chr.len <- end(range)
    
    bin.starts <- seq(from = 1, to = chr.len-binsize, by = stepsize)
    bin.ends <- seq(from = binsize, to = chr.len, by = stepsize)
    
    chr.bins <- GRanges(seqnames=as.character(seqnames(range)), 
                        ranges=IRanges(start=bin.starts, end=bin.ends))
    bins[[i]] <- chr.bins
  }
  binned.data <- unlist(bins, use.names = FALSE)
  binned.data <- subsetByOverlaps(binned.data, gen.ranges)
  
  ## Counts overlaps between bins and our reads
  Watsonreads <- GenomicRanges::countOverlaps(binned.data, reads.gr[strand(reads.gr)=='-']) 
  Crickreads <- GenomicRanges::countOverlaps(binned.data, reads.gr[strand(reads.gr)=='+'])
  bothreads <- Watsonreads + Crickreads
  mcols(binned.data)$bothreads <- bothreads
  mcols(binned.data)$Watsonreads <- Watsonreads
  mcols(binned.data)$Crickreads <- Crickreads
  
  ## Prepare data for plotting
  dfplot.reads <- as.data.frame(binned.data)
  
  my_theme <- theme(
    legend.position="top",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank()
  )
  
  ## Get midpoint values for each genomic bin
  dfplot.reads$midpoint <- dfplot.reads$start + ( (dfplot.reads$end - dfplot.reads$start) %/% 2 )
  
  ## Filter bins with extremely high amount of reads
  Crickreads.outlier <- stats::quantile(dfplot.reads$Crickreads, 0.999)
  Watsonreads.outlier <- stats::quantile(dfplot.reads$Watsonreads, 0.999)
  ## Set outlier bins to the limit
  dfplot.reads$Crickreads[dfplot.reads$Crickreads >= Crickreads.outlier] <- Crickreads.outlier
  dfplot.reads$Watsonreads[dfplot.reads$Watsonreads >= Watsonreads.outlier] <- Watsonreads.outlier
  dfplot.reads$mWatsonreads <- -dfplot.reads$Watsonreads
  
  ## Construct ggplot
  ggplt <- ggplot() + 
    geom_linerange(data=dfplot.reads, aes(ymin=0, ymax=Crickreads, x=midpoint, color='Crick'), size=0.2) +
    geom_linerange(data=dfplot.reads, aes(ymin=0, ymax=mWatsonreads, x=midpoint, color='Watson'), size=0.2) +
    ylab("Read counts") +
    xlab("Genomic region") +
    scale_x_continuous(expand = c(0,0)) +
    scale_color_manual(values = c("paleturquoise4", "sandybrown"), name='Strand') +
    scale_fill_manual(values = c("gray65", "gray85"), name='RegionState') +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    my_theme
  
  ## Prepare Genomic ranges to highlight for plotting
  if (!is.null(highlight.gr) & class(highlight.gr) == 'GRanges') {
    highlight.gr <- subsetByOverlaps(highlight.gr, range(reads.gr, ignore.strand=TRUE), maxgap = 1000000)
    if (length(highlight.gr) > 0) {
      highlight.df <- as.data.frame(highlight.gr)
      ggplt <- ggplt + geom_rect(data=highlight.df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill='gray', alpha=0.25)
    }  
  }
  
  ## Add title if defined
  if (!is.null(title)) {
    if (is.character(title) & nchar(title) > 0) {
      ggplt <- ggplt + ggtitle(title)
    }
  }  
  
  return(ggplt)
}


plotBinnedReads_multiple <- function(composite.files = NULL, region = NULL, highlight.gr = NULL, highlight.pos = NULL, sd.annot='', binsize = 100000, stepsize = 50000, min.mapq = 10, composite.file.field.id = 1, composite.file.field.sep = '_', title = NULL, read.counts.only = FALSE) {
  ## Set stepsize if not defined
  if (is.numeric(stepsize)) {
    if (stepsize < binsize) {
      stepsize == stepsize
    } else {
      stepsize == binsize
    }
  } else {
    stepsize == binsize
  }
  
  reads.plt <- list()
  for (i in seq_along(composite.files)) {
    composite.file <- composite.files[i]
    filename <- basename(composite.file)
    ## Get ID from the composite filename
    file.id <- strsplit(filename, composite.file.field.sep)[[1]][composite.file.field.id]
    ## Remove remaining RData suffix
    file.id <- gsub(file.id, pattern = '\\.RData', replacement = '')
    ## Load composite file reads
    reads.gr <- get(load(composite.file))
    ## Filter reads based on mapping quality
    if (min.mapq > 0) {
      reads.gr <- reads.gr[reads.gr$mapq >= min.mapq]
    }
    
    if (length(reads.gr) > 0) {
      ## Subset to user defined genomic region
      if (!is.null(region) & class(region) == 'GRanges') {
        reads.gr <- subsetByOverlaps(reads.gr, region)
      } else {
        region <- range(reads.gr, ignore.strand=TRUE)
      }
      
      ## Make sure seqlevels reflect data
      reads.gr <- keepSeqlevels(reads.gr, value = unique(as.character(seqnames(reads.gr))), pruning.mode = 'coarse')
      
      ## Bin the genome
      gen.ranges <- region
      bins <- GRangesList()
      for (j in seq_along(gen.ranges)) {
        range <- gen.ranges[j]
        chr.len <- end(range)
        
        bin.starts <- seq(from = 1, to = chr.len-binsize, by = stepsize)
        bin.ends <- seq(from = binsize, to = chr.len, by = stepsize)
        
        chr.bins <- GRanges(seqnames=as.character(seqnames(range)), 
                            ranges=IRanges(start=bin.starts, end=bin.ends))
        bins[[j]] <- chr.bins[width(chr.bins) == binsize]
      }
      binned.data <- unlist(bins, use.names = FALSE)
      binned.data <- subsetByOverlaps(binned.data, gen.ranges)
      
      ## Counts overlaps between bins and our reads
      Watsonreads <- GenomicRanges::countOverlaps(binned.data, reads.gr[strand(reads.gr)=='-']) 
      Crickreads <- GenomicRanges::countOverlaps(binned.data, reads.gr[strand(reads.gr)=='+'])
      bothreads <- Watsonreads + Crickreads
      mcols(binned.data)$bothreads <- bothreads
      mcols(binned.data)$Watsonreads <- Watsonreads
      mcols(binned.data)$Crickreads <- Crickreads
      
      ## Prepare data for plotting
      dfplot.reads <- as.data.frame(binned.data)
      
      ## Get midpoint values for each genomic bin
      dfplot.reads$midpoint <- dfplot.reads$start + ( (dfplot.reads$end - dfplot.reads$start) %/% 2 )
      
      ## Filter bins with extremely high amount of reads
      Crickreads.outlier <- stats::quantile(dfplot.reads$Crickreads, 0.999)
      Watsonreads.outlier <- stats::quantile(dfplot.reads$Watsonreads, 0.999)
      ## Set outlier bins to the limit
      dfplot.reads$Crickreads[dfplot.reads$Crickreads >= Crickreads.outlier] <- Crickreads.outlier
      dfplot.reads$Watsonreads[dfplot.reads$Watsonreads >= Watsonreads.outlier] <- Watsonreads.outlier
      dfplot.reads$mWatsonreads <- -dfplot.reads$Watsonreads
      ## Add file.id
      dfplot.reads$file.id <- file.id
      ## Store binned plotting data
      reads.plt[[i]] <- dfplot.reads
    } else {
      dfplot.reads <- data.frame('seqnames'='', 'start'=start(region), 'end'=end(region), 'width'=width(region), 'strand'='*', 'bothreads'=0, 'Watsonreads'=0, 'Crickreads'=0, 'midpoint'=start(region), 'mWatsonreads'=0, 'file.id'=file.id)
      reads.plt[[i]] <- dfplot.reads
    }
  }
  reads.plt.df <- do.call(rbind, reads.plt)
  ## Set file.id column as factor following order of input composite files
  reads.plt.df$file.id <- factor(reads.plt.df$file.id, levels = unique(reads.plt.df$file.id))
  
  if (read.counts.only) {
    return(reads.plt.df)
    stop('Reporting only binned read counts!!!')
  }
  
  ## Set plotting theme
  my_theme <- theme(
    legend.position="bottom",
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    strip.text.y.left = element_text(angle = 0),
    strip.background =element_blank(),
    panel.spacing.y=unit(0.1, "lines"),
    #plot.margin = unit(c(0, 0, 0, 0), "cm")
  )
  
  ## Construct ggplot
  ggplt <- ggplot() + 
    geom_linerange(data=reads.plt.df, aes(ymin=0, ymax=Crickreads, x=midpoint, color='Crick'), size=0.2) +
    geom_linerange(data=reads.plt.df, aes(ymin=0, ymax=mWatsonreads, x=midpoint, color='Watson'), size=0.2) +
    facet_grid(file.id ~ ., scales = 'free', switch='y') +
    ylab("Read counts: Watson|Crick") +
    xlab("Genomic region") +
    scale_x_continuous(expand = c(0,0), limits = c(start(region), end(region)), labels = comma) +
    scale_color_manual(values = c("paleturquoise4", "sandybrown"), name='Strand') +
    scale_fill_manual(values = c("gray65", "gray85"), name='RegionState') +
    my_theme
  
  ## Add title if defined as x-axis label
  if (!is.null(title)) {
    if (is.character(title) & nchar(title) > 0) {
      #ggplt <- ggplt + ggtitle(title)
      ggplt <- ggplt + xlab(title)
    } else {
      warning("Parameter 'title' can only be character string!!!")
    }  
  }
  
  ## Prepare Genomic ranges to highlight for plotting
  if (!is.null(highlight.gr) & class(highlight.gr) == 'GRanges') {
    highlight.gr <- subsetByOverlaps(highlight.gr, region, type = 'within')
    ## Subset to user defined genomic region
    #if (!is.null(region) & class(region) == 'GRanges') {
    #  highlight.gr <- subsetByOverlaps(highlight.gr, region)
    #}
    if (length(highlight.gr) > 0) {
      highlight.df <- as.data.frame(highlight.gr)
      ggplt <- ggplt + geom_rect(data=highlight.df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill='gray', alpha=0.25)
    }  
  }
  
  ## Prepare Genomic ranges to highlight for plotting
  if (!is.null(highlight.pos) & class(highlight.pos) == 'GRanges') {
    highlight.pos <- subsetByOverlaps(highlight.pos, region, type = 'within')
    ## Subset to user defined genomic region
    #if (!is.null(region) & class(region) == 'GRanges') {
    #  highlight.pos <- subsetByOverlaps(highlight.pos, region)
    #}
    if (length(highlight.pos) > 0) {
      highlight.pos <- getRegionBoundaries(highlight.pos)
      highlight.pos.df <- as.data.frame(highlight.pos)
      ggplt <- ggplt + geom_vline(xintercept = highlight.pos.df$start, linetype='dotted')
    }  
  }
  
  ## Add SD annotation
  if (!is.null(sd.annot) & class(sd.annot) == 'GRanges') {
    # sd.annot.df <- read.table(sd.annot)
    # #sd.annot.df <- sd.annot.df[,c(2,3,4,27)]
    # sd.annot.df <- sd.annot.df[,c(1,2,3,24)]
    # colnames(sd.annot.df) <- c('seqnames', 'start', 'end', 'fracMatch')
    # sd.annot.gr <- makeGRangesFromDataFrame(sd.annot.df, keep.extra.columns = TRUE)
    ## Subset to user defined genomic region
    sd.annot.gr <- subsetByOverlaps(sd.annot, region, type = 'within')
    if (length(sd.annot.gr) > 0) {
      ## Prepare data for plotting
      sd.annot.df <- as.data.frame(sd.annot.gr)
      ## Define SD colors
      sd.categ <- findInterval(sd.annot.df$fracMatch, vec = c(0.95, 0.98, 0.99))
      sd.categ <- dplyr::recode(sd.categ, '0' = '<95%', '1' = '95-98%', '2' = '98-99%', '3'='>=99%')
      sd.categ <- factor(sd.categ, levels=c('<95%', '95-98%', '98-99%', '>=99%'))
      sd.annot.df$sd.categ <- sd.categ
      ## Get colors
      pal <- wes_palette("Zissou1", 4, type = "continuous")
  
      ## Set x axis range
      range.df <- as.data.frame(region)
      annot.plt <- ggplot() +
        geom_rect(data=range.df, aes(xmin=start, xmax=end, ymin=0, ymax=1), fill='white') +
        geom_rect(data=sd.annot.df, aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=sd.categ)) +
        scale_fill_manual(values = pal, name='SDs') +
        scale_x_continuous(expand = c(0,0), labels = comma) +
        my_theme +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="top")
        sd.height = 1/length(composite.files) + 0.05
      ggplt <- plot_grid(annot.plt, ggplt, ncol = 1, align = 'v', axis = 'lr', rel_heights = c(sd.height, 1))
    }  
  }
  return(ggplt)
}
