## Supplemental Code for: Gaps and complex structurally variant loci in phased genome assemblies ##
###################################################################################################
## See below set of function used to process aligned assemblies to the reference genome (T2T-CHM13v1.1).
## See methods section of the above mentioned manuscript for more details
## Make sure to source all function before usage using R function source().

#' Load a de novo assembly aligned to the reference in PAF format into a \code{\link{GRanges-class}} object.
#' 
#' This function will take a PAF file of contig alignments to a reference genome and converts them 
#' into a \code{\link{GRanges-class}} object.
#'
#' @param paf.file Alignments reported in PAF file format.
#' @param index A user defined identifier for given PAF alignments. 
#' @param min.mapq Alignments of mapping quality lower then this threshold will be removed.
#' @param min.ctg.size A minimum length a final contig after gaps are collapsed.
#' @param report.ctg.ends Set to \code{TRUE} if aligned position of each contig ends should be reported
#' @param min.ctg.ends A minimum length of alignment to be considered when reporting contig end alignments.
#' @return A \code{\link{GRanges-class}} object.
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @author David Porubsky
#' @export
#' 
paf2ranges <- function(paf.file=NULL, index=NULL, min.mapq=10, min.aln.width=10000, min.ctg.size=500000, report.ctg.ends=FALSE, min.ctg.ends=50000) {
  ## Get total processing time
  #ptm <- proc.time()
  
  if (file.exists(paf.file)) {
    ptm <- startTimedMessage("\nLoading PAF file: ", paf.file)
    paf <- utils::read.table(paf.file, stringsAsFactors = FALSE, comment.char = '&')
    ## Keep only first 12 columns
    paf <- paf[,c(1:12)]
    ## Add header
    header <- c('q.name', 'q.len', 'q.start', 'q.end', 'strand', 't.name', 't.len', 't.start', 't.end', 'n.match', 'aln.len', 'mapq') 
    colnames(paf) <- header
    stopTimedMessage(ptm)
  } else {
    stop(paste0("PAF file ", bedfile, " doesn't exists !!!"))
  }  
  ## Filter by mapping quality
  if (min.mapq > 0) {
    ptm <- startTimedMessage("    Keeping alignments of min.mapq: ", min.mapq)
    paf <- paf[paf$mapq >= min.mapq,]
    stopTimedMessage(ptm)
  }
  if (nrow(paf) == 0) {
    stop("None of the PAF alignments reach user defined mapping quality (min.mapq) !!!")
  }
  ## Convert data.frame to GRanges object
  #paf.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, keep.extra.columns = FALSE, seqnames.field = 'q.name', start.field = 'q.start', end.field = 'q.end')
  paf.gr <- GenomicRanges::makeGRangesFromDataFrame(paf, keep.extra.columns = FALSE, seqnames.field = 't.name', start.field = 't.start', end.field = 't.end')
  mcols(paf.gr) <- paf[,c('q.len', 't.len', 'n.match', 'aln.len', 'mapq', 'q.name')]
  names(paf.gr) <- NULL
  #paf.gr$target.gr <- GenomicRanges::GRanges(seqnames=paf$t.name, ranges=IRanges::IRanges(start=paf$t.start, end=paf$t.end), 'q.name'= paf$q.name)
  paf.gr$query.gr <- GenomicRanges::GRanges(seqnames=paf$q.name, ranges=IRanges::IRanges(start=paf$q.start, end=paf$q.end), 't.name'= paf$t.name)
  ## Add index if defined
  if (!is.null(index) & is.character(index)) {
    paf.gr$ID <- index
  }
  ## Ignore strand
  #GenomicRanges::strand(paf.gr) <- '*'
  ## Filter out small alignments
  if (min.aln.width > 0) {
    ptm <- startTimedMessage("    Keeping alignments of min width: ", min.aln.width, 'bp')
    #paf.gr <- paf.gr[width(paf.gr) >= min.aln.width]
    paf.gr <- paf.gr[paf.gr$aln.len >= min.aln.width]
    stopTimedMessage(ptm)
  }
  if (length(paf.gr) == 0) {
    stop("None of the PAF alignments reach user defined alignment size (min.aln.width !!!")
  }
  ## Filter out small contig sizes
  if (min.ctg.size > 0) {
    ptm <- startTimedMessage("    Keeping contigs of min size: ", min.ctg.size, 'bp')
    paf.gr <- paf.gr[paf.gr$q.len >= min.ctg.size]
    stopTimedMessage(ptm)
  }
  ## Keep only seqlevels that remained after data filtering
  if (length(paf.gr) > 0) {
    paf.gr <- GenomeInfoDb::keepSeqlevels(paf.gr, value = unique(as.character(GenomeInfoDb::seqnames(paf.gr))))
    paf.gr$query.gr <- GenomeInfoDb::keepSeqlevels(paf.gr$query.gr, value = unique(as.character(GenomeInfoDb::seqnames(paf.gr$query.gr))))
  } else {
    stop("None of the PAF alignments reach user defined contig size (min.ctg.size) !!!")
  }
  
  if (report.ctg.ends == TRUE) {
    ptm <- startTimedMessage("    Reporting end positions of each contig")
    #paf.grl <- split(paf.gr, GenomeInfoDb::seqnames(paf.gr))
    paf.grl <- split(paf.gr, GenomeInfoDb::seqnames(paf.gr$query.gr))
    to.collapse <- which(lengths(paf.grl) > 1)
    #ctg.ends.grl <- S4Vectors::endoapply( paf.grl[to.collapse], function(x) range(sort(x, ignore.strand=TRUE)$target.gr) )
    
    ## Helper function
    gr2ranges <- function(gr, min.ctg.ends=0, allowed.ctg.length=0.05) {
      gr <- GenomeInfoDb::keepSeqlevels(gr, value = unique(as.character(seqnames(gr))))
      q.name <- unique(gr$q.name)
      if (length(q.name) > 1) {
        stop("Submitted GRanges object contains more than one 'q.name', not allowed !!!")
      }
      if (min.ctg.ends > 0) {
        gr <- gr[gr$aln.len >= min.ctg.ends]
      }
      if (length(gr) > 0) {
        ## Order by contig alignments
        gr <- gr[order(gr$query.gr)]
        if (length(gr) > 1) {
          ## Report target ranges corresponding to the contig ends
          gr.ends <- gr[c(1, length(gr))]
          ## Report genomic range for contig ends mapping
          gr.range <- range(gr.ends, ignore.strand=TRUE)
          ## Get best contiguous alignment ##
          ## Report target sequence with the highest number of continuously aligned bases
          target2keep <- names(which.max(sum(split(width(gr), seqnames(gr)))))
          gr.contig <- gr[seqnames(gr) == target2keep]
          gr.contig <- gr.contig[,0]
          strand(gr.contig) <- '*'
          gr.gaps <- gaps(gr.contig, start = min(start(gr.contig)))
          if (length(gr.gaps) > 0) {
            ## Remove gaps longer than query(contig) length
            q.len <- unique(gr$q.len)
            ## Remove gaps longer than total alignment length
            #aln.len <- sum(gr$aln.len)
            ## Adjust allowed contig length based on fraction allowed to be added to the contigs size
            if (!is.null(allowed.ctg.length)) {
              if (allowed.ctg.length > 0) {
                q.len <- q.len + (q.len * allowed.ctg.length)
                #aln.len <- aln.len + (aln.len * allowed.ctg.length)
              }
            }
            gr.gaps <- gr.gaps[width(gr.gaps) < q.len]
            #gr.gaps <- gr.gaps[width(gr.gaps) < aln.len]
            ## Keep only gaps that together with contig length are no longer than query(contig) length
            ctg.size <- sum(width(reduce(gr.contig)))
            gap.size <- sum(width(reduce(gr.gaps)))  
            while ((ctg.size + gap.size) > q.len & length(gr.gaps) > 0) {
              #while ((ctg.size + gap.size) > aln.len & length(gr.gaps) > 0) {  
              ## Remove longest gap
              gr.gaps <- gr.gaps[-which.max(width(gr.gaps))]
              gap.size <- sum(width(gr.gaps))
            }
          }  
          ## Allow only gaps that are no longer than smallest alignment
          #gr.gaps <- gr.gaps[width(gr.gaps) <= min(gr$aln.len)]
          ## Collapse contiguous ranges
          gr.contig <- reduce(c(gr.contig, gr.gaps))
          ## Add contig names and split alignment ids
          if (length(gr.contig) > 1) {
            gr.contig$q.name <- q.name
            aln <- paste0('p', 1:length(gr.contig))
            gr.contig$aln <- aln
          } else {
            gr.contig$q.name <- q.name
            gr.contig$aln <- 'single'
          }  
          #gr.contig <- gr.contig[which.max(width(gr.contig))]
          # gr.ranges <- range(gr[order(gr$query.gr)], ignore.strand=TRUE)
          # if (report.longest.range) {
          #   return(gr.ranges[which.max(width(gr.ranges))])
          # } else {
          #   return(gr.ranges)
          # }
          gr.contig$ends <- paste(as.character(gr.range), collapse = ';')
          #gr.range$longest.aln <- gr.contig
        } else {
          #gr.range <- gr[,0]
          #strand(gr.range) <- '*'
          #gr.range$longest.aln <- gr.range
          gr.contig <- gr[,0]
          strand(gr.contig) <- '*'
          gr.contig$q.name <- q.name
          gr.contig$aln <- 'single'
          gr.contig$ends <- paste(as.character(gr.contig), collapse = ';')
        }  
      } else {
        #gr.range <- GRanges()
        gr.contig <- GRanges()
      }  
      #return(gr.range)
      return(gr.contig)
    }
    
    #ctg.ends.grl <- S4Vectors::endoapply( paf.grl[to.collapse], function(x) range(x[order(x$query.gr)], ignore.strand=TRUE) )
    
    ## Report contigs alignment ranges only for alignments of a certain size (default: 50kb)
    if (min.ctg.ends > 0) {
      ctg.ends.grl <- S4Vectors::endoapply( paf.grl[to.collapse], 
                                            function(x) gr2ranges(gr = x, min.ctg.ends = min.ctg.ends) )  
    } else {
      ctg.ends.grl <- S4Vectors::endoapply( paf.grl[to.collapse], 
                                            function(x) gr2ranges(gr = x) )
    }
    
    simple.ends.gr <- unlist(paf.grl[-to.collapse], use.names = FALSE)
    if (length(simple.ends.gr) > 0) {
      #simple.ends <- as.character(simple.ends.gr$target.gr)
      simple.ends <- as.character(simple.ends.gr)
      #names(simple.ends) <- unique(as.character(GenomeInfoDb::seqnames(simple.ends.gr)))
      names(simple.ends) <- unique(as.character(GenomeInfoDb::seqnames(simple.ends.gr$query.gr)))
      
      simple.gr <- simple.ends.gr[,0]
      strand(simple.gr) <- '*'
      #simple.gr$longest.aln <- simple.gr
      simple.gr$q.name <- as.character(seqnames(simple.ends.gr$query.gr))
      simple.gr$aln <- 'single'
      simple.gr$ends <- as.character(simple.gr)
      names(simple.gr) <- NULL
      
      split.ends <- as.character(unlist(ctg.ends.grl))
      split.ends.l <- split(split.ends, names(split.ends))
      split.ends <- sapply(split.ends.l, function(x) paste(x, collapse = ';'))
    }  
    
    split.ends.gr <- unlist(ctg.ends.grl)
    #split.ends.gr$q.name <- names(split.ends.gr)
    names(split.ends.gr) <- NULL
    
    #ctg.ends <- c(simple.ends, split.ends)
    #paf.gr$ctg.end.pos <- ctg.ends[match(as.character(seqnames(paf.gr)), names(ctg.ends))]
    #paf.gr$ctg.end.pos <- ctg.ends[match(as.character(seqnames(paf.gr$query.gr)), names(ctg.ends))]
    
    ends.gr <- sort(c(simple.gr, split.ends.gr))
    ends.gr$ID <- unique(paf.gr$ID)
    
    stopTimedMessage(ptm)
    return(list('ctg.aln'=paf.gr, 'ctg.ends'=ends.gr))
  } else {
    return(list('ctg.aln'=paf.gr, 'ctg.ends'=NULL))
  }
  
  ## Report total processing time
  #time <- proc.time() - ptm
  #message("Total time: ", round(time[3],2), "s")
  
  #return(paf.gr)
} 

#' This function will take in a \code{\link{GRanges-class}} or \code{\link{GRangesList-class}} object of alignments of a 
#' single contig to a reference and collapses alignment gaps based user specified maximum allowed gap.
#'
#' @param ranges A \code{\link{GRanges-class}} or \code{\link{GRangesList-class}} object of regions of a single or multiple contigs aligned to a reference.
#' @param max.gap A maximum length of a gap within a single contig alignments to be collapsed.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
collapseGaps <- function(ranges, max.gap=100000) {
  ## Helper function definition
  fillGaps <- function(gr, max.gap=100000) {
    gr <- GenomeInfoDb::keepSeqlevels(gr, value = as.character(unique(GenomeInfoDb::seqnames(gr))), pruning.mode = 'coarse')
    gr <- GenomicRanges::sort(gr)
    gap.gr <- GenomicRanges::gaps(gr, start = min(start(gr)))
    gap.gr <- gap.gr[width(gap.gr) <= max.gap]
    if (length(gap.gr) > 0) {
      red.gr <- GenomicRanges::reduce(c(gr[,0], gap.gr))
      GenomicRanges::mcols(red.gr) <- GenomicRanges::mcols(gr)[length(gr),]
    } else {
      red.gr <- GenomicRanges::reduce(gr)
      GenomicRanges::mcols(red.gr) <- GenomicRanges::mcols(gr)[length(gr),]
    }
    return(red.gr)
  }
  
  if (class(ranges) == "CompressedGRangesList") {
    ## Process only contigs with split alignments
    to.fill <- which(lengths(ranges) > 1)
    red.ranges <-  suppressWarnings( S4Vectors::endoapply(ranges[to.fill], function(gr) fillGaps(gr=gr, max.gap=max.gap)) )
    red.ranges <- unlist(red.ranges, use.names = FALSE)
    ## Add ranges with no split alignment
    red.ranges <- c(unlist(ranges[-to.fill], use.names = FALSE), red.ranges)
  } else if (class(ranges) == "GRanges") {
    red.ranges <- fillGaps(gr=ranges, max.gap = max.gap)
  } else {
    stop("Only objects of class 'GRanges' or 'GRangesList' are allowed as input for parameter 'ranges' !!!")
  }
  return(red.ranges)
}

#' This function will take in a \code{\link{GRanges-class}} or \code{\link{GRangesList-class}} object of genomic ranges of a 
#' single a contig aligned to a reference and reports all alignment gaps.
#'
#' @param id.col A column number from the original \code{\link{GRanges-class}} object to be reported as an unique ID.
#' @inheritParams collapseGaps
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
reportGaps <- function(ranges, id.col=NULL) {
  ## Helper function definitions
  getGaps <- function(gr=NULL) {
    gap.gr <- GenomicRanges::gaps(gr, start = min(start(gr)))
    gap.gr <- gap.gr[GenomicRanges::strand(gap.gr) == '*']
    return(gap.gr)
  }  
  
  processGaps <- function(gr, id.col=NULL) {
    ## Make sure only seqlevels present in the submitted ranges are kept
    gr <- GenomeInfoDb::keepSeqlevels(gr, value = as.character(unique(GenomeInfoDb::seqnames(gr))), pruning.mode = 'coarse')
    ## Sort ranges by position
    GenomicRanges::strand(gr) <- '*'
    gr <- GenomicRanges::sort(gr)
    ## Make sure rows are not named
    names(gr) <- NULL
    ## Keep only ranges with the same seqnames
    if (length(GenomeInfoDb::seqlevels(gr)) > 1) {
      max.seqname <- GenomeInfoDb::seqlevels(gr)[which.max(S4Vectors::runLength(GenomeInfoDb::seqnames(gr)))]
      gr <- gr[GenomeInfoDb::seqnames(gr) == max.seqname]
      warning("Multiple 'seqlevels' present in submitted ranges, keeping only ranges for: ", max.seqname)
    }
    ## TODO for alignment landing on different chromosomes report gap == 0 and both alignments
    if (length(gr) > 1) {
      ## Calculate gaps
      #gap.gr <- GenomicRanges::gaps(gr, start = min(start(gr)))
      #gap.gr <- gap.gr[strand(gap.gr) == '*']
      gap.gr <- getGaps(gr)
      
      ## Add ID column from the original gr object if defined
      if (!is.null(id.col)) {
        if (id.col > 0 & ncol(GenomicRanges::mcols(gr)) >= id.col) {
          GenomicRanges::mcols(gap.gr) <- rep(unique(GenomicRanges::mcols(gr)[id.col]), length(gap.gr))
        } else {
          warning("User defined 'id.col' number is larger the then total number of columns in input 'gr', skipping adding id column ...")
        }
      }
            
      ## Report alignment on each side of the gap
      if (length(gap.gr) > 0) {
        ## Report upstream and downstream target ranges
        ## Upstream
        up.idx <- IRanges::follow(gap.gr, gr)
        up.idx.keep <- !is.na(up.idx)
        gap.gr$up.gr <- GenomicRanges::GRanges(seqnames=GenomeInfoDb::seqnames(gap.gr), 
                                               ranges=IRanges::IRanges(start=GenomicRanges::start(gap.gr), end=GenomicRanges::start(gap.gr)))
        gap.gr$up.gr[up.idx.keep] <- gr[up.idx[up.idx.keep]][,0]
        ## Downstream
        down.idx <- IRanges::precede(gap.gr, gr)
        down.idx.keep <- !is.na(down.idx)
        gap.gr$down.gr <- GenomicRanges::GRanges(seqnames=GenomeInfoDb::seqnames(gap.gr), 
                                                 ranges=IRanges::IRanges(start=GenomicRanges::end(gap.gr), end=GenomicRanges::end(gap.gr)))
        gap.gr$down.gr[down.idx.keep] <- gr[down.idx[down.idx.keep]][,0]
        
        ## Report upstream and downstream query ranges and gaps if defined
        if ('query.gr' %in% names(mcols(gr)) & class(gr$query.gr) == 'GRanges') {
          ## Upstream
          gap.gr$query.up.gr <- GenomicRanges::GRanges(seqnames=GenomeInfoDb::seqnames(gap.gr), 
                                                       ranges=IRanges::IRanges(start=GenomicRanges::start(gap.gr), end=GenomicRanges::start(gap.gr)))
          suppressWarnings( gap.gr$query.up.gr[up.idx.keep] <- gr$query.gr[up.idx[up.idx.keep]][,0] )
          ## Downstream
          gap.gr$query.down.gr <- GenomicRanges::GRanges(seqnames=GenomeInfoDb::seqnames(gap.gr), 
                                                         ranges=IRanges(start=GenomicRanges::end(gap.gr), end=GenomicRanges::end(gap.gr)))
          suppressWarnings( gap.gr$query.down.gr[down.idx.keep] <- gr$query.gr[down.idx[down.idx.keep]][,0] )
          ## Calculate gaps
          gap.grl <- split(gap.gr, 1:length(gap.gr))
          query.gap.grl <- S4Vectors::endoapply(gap.grl, function(x) getGaps(c(x$query.up.gr, x$query.down.gr)))
          query.gap.gr <- unlist(query.gap.grl, use.names = FALSE)
          ## Initialize gap ranges
          gap.gr$query.gap.gr <- GenomicRanges::GRanges(seqnames=rep('unknown', length(gap.gr)), 
                                                        ranges=IRanges::IRanges(start=1, end=1))
          if (length(query.gap.gr) > 0) {
            gap.idx <- which(lengths(query.gap.grl) > 0)
            gap.gr$query.gap.gr[gap.idx] <- query.gap.gr
          }
        }
        return(gap.gr)
      } else {
        ## Report upstream and downstream query ranges and gaps if defined
        if ('query.gr' %in% names(mcols(gr)) & class(gr$query.gr) == 'GRanges') {
          dummy.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1), ID='dummy')
          dummy.gr$up.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
          dummy.gr$down.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
          dummy.gr$query.up.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
          dummy.gr$query.down.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
          dummy.gr$query.gap.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
        } else {
          dummy.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1), ID='dummy')
          dummy.gr$up.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
          dummy.gr$down.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
        }  
        return(dummy.gr)
      }
    } else {
      ## Report upstream and downstream query ranges and gaps if defined
      if ('query.gr' %in% names(mcols(gr)) & class(gr$query.gr) == 'GRanges') {
        dummy.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1), ID='dummy')
        dummy.gr$up.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
        dummy.gr$down.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
        dummy.gr$query.up.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
        dummy.gr$query.down.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
        dummy.gr$query.gap.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
      } else {
        dummy.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1), ID='dummy')
        dummy.gr$up.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
        dummy.gr$down.gr <- GRanges(seqnames = 'dummy', ranges = IRanges(start=1, end=1))
      }  
      return(dummy.gr)
      #return(GRanges())
    }
  }  
  
  if (class(ranges) == "CompressedGRangesList") {
    ## Keep only contigs with split alignments
    ptm <- startTimedMessage("Reporting gaps")
    
    ranges <- ranges[lengths(ranges) > 1]
    gaps <-  suppressWarnings( S4Vectors::endoapply(ranges, function(gr) processGaps(gr=gr, id.col=id.col)) )
    gaps <- unlist(gaps, use.names = FALSE)
    ## Remove empty ranges
    gaps <- gaps[seqnames(gaps) != 'dummy']
    #gaps <-  suppressWarnings( S4Vectors::lapply(ranges, function(gr) processGaps(gr=gr, id.col=id.col)) )
    #do.call(c, gaps)
    
    stopTimedMessage(ptm)
  } else if (class(ranges) == "GRanges") {
    ptm <- startTimedMessage("Reporting gaps")
    
    gaps <- processGaps(gr=ranges, id.col=id.col)
    
    stopTimedMessage(ptm)
  } else {
    stop("Only objescts of class 'GRanges' or 'GRangesList' are allowed as input for parameter 'ranges' !!!")
  }
  return(gaps)
}


#' This function will take in a \code{\link{GRanges-class}} or \code{\link{GRangesList-class}} object of genomic ranges of a 
#' single a contig aligned to a reference and reports coverage of regions that overlaps each other. 
#'
#' @inheritParams collapseGaps
#' @inheritParams reportGaps
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
reportContigCoverage <- function(ranges, id.col=NULL) {
  ## Helper function definitions
  processContigCoverage <- function(gr=NULL, id.col=NULL) {
    ## Get disjoin ranges
    gr.disjoint <- GenomicRanges::disjoin(gr, ignore.strand=TRUE)
    ## Add ID column from the original gr object if defined
    if (!is.null(id.col)) {
      if (id.col > 0 & ncol(GenomicRanges::mcols(gr)) >= id.col) {
        GenomicRanges::mcols(gr.disjoint) <- rep(unique(GenomicRanges::mcols(gr)[id.col]), length(gr.disjoint))
      } else {
        warning("User defined 'id.col' number is larger the then total number of columns in input 'gr', skipping adding id column ...")
      }
    }
    ## Report coverage of each disjoint range
    gr.disjoint$cov <- IRanges::countOverlaps(gr.disjoint, gr)
    ## Export final ranges with coverage
    return(gr.disjoint)
  }
  
  if (class(ranges) == "CompressedGRangesList") {
    covs <-  suppressWarnings( S4Vectors::endoapply(ranges, function(gr) processContigCoverage(gr=gr, id.col=id.col)) )
    covs <- unlist(covs, use.names = FALSE)
  } else if (class(ranges) == "GRanges") {
    covs <- processContigCoverage(gr=ranges, id.col=id.col)
  } else {
    stop("Only objects of class 'GRanges' or 'GRangesList' are allowed as input for parameter 'ranges' !!!")
  }
  return(covs)
}

#' This function will take in a \code{\link{GRanges-class}} object of genomic ranges and enumerate the number of overlapping bases
#' with other set of user defined 'query' genomic ranges.
#'
#' @param ranges A \code{\link{GRanges-class}} object of genomic regions to get the number of overlapping bases with 'query.ranges'.
#' @param query.ranges A \code{\link{GRanges-class}} object of genomic regions for which one want to report overlaps. 
#' @param index A user defined name of a column where the number of overlapping bases between 'ranges' and 'query.ranges' will be reported.
#' @return A \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#'
reportOverlapBases <- function(ranges=NULL, query.ranges=NULL, index=NULL) {
  if (!is.null(query.ranges) & !is.null(ranges)) {
    if (class(ranges) == "GRanges" & class(query.ranges) == "GRanges") {
      ## Get disjoint ranges between query and user.defined set of ranges
      disj.gr <- suppressWarnings( GenomicRanges::disjoin(c(ranges[,0], query.ranges[,0])) )
      ## Get disjoint ranges that are overlapping with query.ranges
      disj.query.gr <- IRanges::subsetByOverlaps(disj.gr, query.ranges)
      ## Get disjoint ranges overlapping with ranges of interest
      disj.query.roi.gr <- IRanges::subsetByOverlaps(disj.query.gr, ranges)
      ## Split by regions of interest
      hits <- IRanges::findOverlaps(ranges, disj.query.roi.gr)
      disj.query.roi.grl <- split(disj.query.roi.gr[S4Vectors::subjectHits(hits)], S4Vectors::queryHits(hits))
      query.bases <- sapply(disj.query.roi.grl, function(gr) sum(width(reduce(gr))))
      ## Add overlapping bases counts
      if (!is.null(index) & nchar(index) > 0) {
        new.col.idx <- ncol(GenomicRanges::mcols(ranges)) + 1
        GenomicRanges::mcols(ranges)[new.col.idx] <- 0
        colnames(GenomicRanges::mcols(ranges))[new.col.idx] <- index
        GenomicRanges::mcols(ranges)[new.col.idx][unique(S4Vectors::queryHits(hits)),] <- query.bases
      } else {
        ranges$query.bases <- 0
        ranges$query.bases[unique(S4Vectors::queryHits(hits))] <- query.bases
      } 
    }
  }
  return(ranges)
}  


#' Summarize coverage of de novo assembly aligned to the reference.
#' 
#' This function will take a PAF file of contig alignments to a reference genome then summarize the coverage
#' of these alignment with respect to the reference and lastly reports a \code{\link{GRanges-class}} summary 
#' along with genome-wide visualisation.
#'
#' @param h1.paf Alignments reported in PAF file format for haplotype 1.
#' @param h2.paf Alignments reported in PAF file format for haplotype 1.
#' @inheritParams paf2ranges
#' @return A  \code{\link{GRanges-class}} object.
#' @author David Porubsky
#' @export
#' 
asm2referenceCoverage <- function(h1.paf = NULL, h2.paf = NULL, min.mapq = 0, min.aln.width = 0, min.ctg.size = 0, index = NULL) {
  ## Read in alignments to the reference in PAF format for haplotype1
  h1.gr <- paf2ranges(paf.file = h1.paf
                      , index = 'H1'
                      , min.mapq = min.mapq
                      , min.aln.width = min.aln.width
                      , min.ctg.size = min.ctg.size
                      , report.ctg.ends = FALSE)
  
  ## Read in alignments to the reference in PAF format for haplotype2
  h2.gr <- paf2ranges(paf.file = h2.paf
                      , index = 'H2'
                      , min.mapq = min.mapq
                      , min.aln.width = min.aln.width
                      , min.ctg.size = min.ctg.size
                      , report.ctg.ends = FALSE)
  
  gr <- suppressWarnings( c(h1.gr$ctg.aln, h2.gr$ctg.aln) )
  ## Calculate coverage
  gr.cov <- as( coverage(gr), 'GRanges')
  cov.hap <- sum(width(gr.cov[gr.cov$score == 1]))
  cov.dip <- sum(width(gr.cov[gr.cov$score == 2]))
  cov.multicov <- sum(width(gr.cov[gr.cov$score > 2]))
  
  if (is.null(index)) {
    index <- 'NA'
  }
  
  ## Export ploidy
  cov.regions <- gr.cov[gr.cov$score > 0]
  names(mcols(cov.regions)) <- 'ploidy'
  cov.regions$index <- index
  seqlengths(cov.regions) <- NA
  
  ## Make a ploidy plot
  plt.df <- as.data.frame(cov.regions)
  plt.df$ploidy.categ <- '2n' 
  plt.df$ploidy.categ[plt.df$ploidy == 1] <- '1n'
  plt.df$ploidy.categ[plt.df$ploidy > 2] <- '>2n'
  plt.df$seqnames <- factor(plt.df$seqnames, levels = gtools::mixedsort(as.character(unique(plt.df$seqnames))))
  my.theme <- theme(legend.position="bottom",
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    strip.text.y.left = element_text(angle = 0),
                    strip.background =element_blank())
  plt <- ggplot() +
    geom_rect(data=plt.df, aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=ploidy.categ)) +
    scale_fill_manual(values = c('1n' = 'deepskyblue2', '2n' = 'darkolivegreen3', '>2n' = 'brown3')) +
    ggnewscale::new_scale_fill() +
    #geom_rect(data=SDs.df, aes(xmin=start, xmax=end, ymin=1, ymax=1.5), fill='orange') +
    #geom_rect(data=cent.df , aes(xmin=start, xmax=end, ymin=1, ymax=1.5), fill='black') +
    scale_x_continuous(expand = c(0,0)) +
    facet_grid(seqnames ~ ., switch='y') +
    my.theme
  if (index != 'NA') {
    plt <- plt + ggtitle(index)
  }
  ## Return results
  return(list(plot = plt, ploidy.ranges = cov.regions))
}

      
## Time messaging functions ##
##############################
messageU <- function(..., underline='=', overline='=') {
  
  x <- paste0(..., collapse='')
  if (!is.null(overline)) {
    message(rep(overline, nchar(x)))
  }
  message(x)
  if (!is.null(underline)) {
    message(rep(underline, nchar(x)))
  }
  
}


startTimedMessage <- function(...) {
  
  x <- paste0(..., collapse='')
  message(x, appendLF=FALSE)
  ptm <- proc.time()
  return(ptm)
  
}


stopTimedMessage <- function(ptm) {
  
  time <- proc.time() - ptm
  message(" ", round(time[3],2), "s")
  
}
