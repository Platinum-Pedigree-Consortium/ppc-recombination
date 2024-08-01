#' Find meiotic recombination breakpoints using circular binary segmentation (CBS)
#'
#' This function takes as an input a \code{\link{VRanges-class}} object with extra meta-column that contains
#' for each variant parental identity. P1 -> parent1 & P2 -> parent2.
#'
#' @param comparison A \code{\link{VRanges-class}} object
#' @param minSeg Minimal length (number of variants) being reported as haplotype block (\code{fastseg} parameter).
#' @param minSegBp Minimal length of reported inherited parental segment.
#' @param smooth Number of consecutive variants being considered as a random error and so being corrected (flipped).
#' @param collapse.amb Set to \code{TRUE} if segments with ambiguous haplotype assignments should be collapsed.
#' @param mask.gr A \code{\link{GRanges-class}} object containing genomic regions that should be masked
#' @param mask.overlap.th Defines a maximum allowed overlap between reported haplotype segments and a set of masked regions
#' defined in `mask.gr` [default: 50%]
#' @importFrom fastseg fastseg
#' @importFrom dplyr recode
#' @author David Porubsky
#' @export

findRecomb <- function(comparison=NULL, minSeg=NULL, minSegBp=100000, smooth=3, collapse.amb=TRUE, mask.gr=NULL, mask.overlap.th = 50) {
  
  ## Helper function
  switchValue <- function(x) {
    if (x == 1) {
      x <- 0
    } else {
      x <- 1
    }
  }
  
  # collapseBins <- function(gr, id.field=0) {
  #   ind.last <- cumsum(runLength(Rle(mcols(gr)[,id.field]))) ##get indices of last range in a consecutive(RLE) run of the same value
  #   ind.first <- c(1,cumsum(runLength(Rle(mcols(gr)[,id.field]))) + 1) ##get indices of first range in a consecutive(RLE) run of the same value
  #   ind.first <- ind.first[-length(ind.first)]  ##erase last index from first range indices
  #   collapsed.gr <- GRanges(seqnames=seqnames(gr[ind.first]), ranges=IRanges(start=start(gr[ind.first]), end=end(gr[ind.last])), mcols=mcols(gr[ind.first]))
  #   names(mcols(collapsed.gr)) <- names(mcols(gr[ind.first]))
  #   return(collapsed.gr)
  # }
  
  ## Remove missing values
  comparison.filt <- comparison[mcols(comparison)[,1] != 'N']
  ## Define the minimum segment size for fastseq function as median variant count per minSegBp
  if (is.null(minSeg)) {
    #minSeg <- ceiling(length(comparison.filt) * 0.01)
    binned.gr <- tileGenome(seqlengths = seqlengths(comparison.filt), tilewidth = minSegBp, cut.last.tile.in.chrom = TRUE)
    binned.gr$n.vars <- countOverlaps(binned.gr, comparison.filt)
    minSeg <- median(binned.gr$n.vars)
  }  
  
  if (length(comparison.filt) >= 2*minSeg) {
    ## Recode comparison into in 0/1 vector
    comparison.filt$comp.vector <- dplyr::recode(mcols(comparison.filt)[,1], 'P1' = 0, 'P2' = 1)
    
    ## Run CBS segmentation on 0/1 vector
    segs <- fastseg(comparison.filt$comp.vector, minSeg=minSeg, segMedianT = c(0.8, 0.2))
    #segs <- fastseg(comparison.filt$comp.vector, minSeg=minSeg)
    
    while (any(segs$num.mark <= smooth)) {
      toSwitch <- which(segs$num.mark <= smooth)
      switch.segs <- segs[toSwitch]
      switch.pos <- mapply(function(x,y) {x:y}, x=switch.segs$startRow, y=switch.segs$endRow)
      switch.pos <- unlist(switch.pos)
      
      switched.vals <- sapply(comparison.filt$comp.vector[switch.pos], switchValue) #SWITCH
      #comparison$comp.vector <- comparison$comp.vector[-switch.pos]  #DELETE
      comparison.filt$comp.vector[switch.pos] <- switched.vals
      
      ## Rerun fastseg
      segs <- fastseg(comparison.filt$comp.vector, minSeg=minSeg, segMedianT = c(0.8, 0.2))
      segs <- fastseg(comparison.filt$comp.vector, minSeg=minSeg)
    }
    ## Remove ambiguous
    # toRemove <- which(segs$seg.mean > 0.2 & segs$seg.mean < 0.8)
    # if (length(toRemove) > 1) {
    #   remove.segs <- segs[toRemove]
    #   remove.pos <- mapply(function(x,y) {x:y}, x=remove.segs$startRow, y=remove.segs$endRow)
    #   remove.pos <- unlist(remove.pos)
    #   comparison.filt <- comparison.filt[-remove.pos]
    #   segs <- fastseg(comparison.filt$comp.vector, minSeg=minSeg)
    # } 
    
    ## Assign genomic positions
    #gen.ranges <- IRanges(start=start(comparison.filt)[segs$startRow], end=end(comparison.filt)[segs$endRow])
    gen.ranges <- IRanges(start=start(comparison.filt)[segs$startRow], end=start(comparison.filt)[segs$endRow])
    ranges(segs) <- gen.ranges
    ## Add chromosome name and length
    seqlevels(segs) <- unique(as.character(seqnames(comparison.filt)))
    seqlengths(segs) <- seqlengths(comparison.filt)
    
    ## Assign haplotype
    segs$match[segs$seg.mean <= 0.25] <- 'hap1'
    segs$match[segs$seg.mean >= 0.75] <- 'hap2'
    segs$match[segs$seg.mean > 0.25 & segs$seg.mean < 0.75] <- 'amb'
    # segs$match[segs$seg.mean <= 0.2] <- 'hap1'
    # segs$match[segs$seg.mean >= 0.8] <- 'hap2'
    # segs$match[segs$seg.mean > 0.2 & segs$seg.mean < 0.8] <- 'amb'
    # segs$match[segs$seg.mean <= 0.2] <- 'hap1'
    # segs$match[segs$seg.mean >= 0.8] <- 'hap2'
    # segs$match[segs$seg.mean > 0.2 & segs$seg.mean < 0.8] <- 'amb'
    
    ## Remove segments with mixed H1 and H2 signal
    if (collapse.amb) {
      segs <- segs[segs$match != 'amb']
    }
    
  } else {
    message("    Low density of informative SNVs, skipping ...")
    ## Report empty GRanges object
    segs <- GenomicRanges::GRanges()
  }
  segm <- segs
  
  ## Merge consecutive haplotype segments
  if (length(segs) > 1) {
    segm <- primatR::collapseBins(gr = segm, id.field = 6, measure.field = 2)
  } else {
    segm <- segm
  }
  
  ## Make sure first and last haplotype segment reaches start and end of the chromosome
  if (length(segm) > 0) {
    start(segm)[which.min(start(segm))] <- 1
    if (!is.na(seqlengths(segm))) {
      end(segm)[which.max(end(segm))] <- seqlengths(segm)
    }
  }
  
  ## Remove segments shorter than minSegBp
  if (!is.null(minSegBp) & minSegBp > 0) {
    segm <- segm[width(segm) >= minSegBp]
  }
  
  ## Remove segments with less variants than defined in minSeg
  if (!is.null(minSeg) & minSeg > 0) {
    segm <- segm[segm$num.mark >= minSeg]
  }
  
  ## Merge consecutive haploype segments
  if (length(segm) > 1) {
    segm <- primatR::collapseBins(gr = segm, id.field = 6, measure.field = 2)
  } else {
    segm <- segm
  }
  
  ## Remove segments that overlap masked region by a user defined overlap
  if (!is.null(mask.gr) & length(segm) > 0) {
    if (methods::is(mask.gr, 'list')) {
      for (i in seq_along(mask.gr)) {
        tmp <- primatR::getRangesOverlaps(query.gr = segm, subject.gr = mask.gr[[i]], index = 'MASK')
        segm <- segm[tmp$PercOverlap_MASK <= mask.overlap.th]
      }
    } else {
      tmp <- primatR::getRangesOverlaps(query.gr = segm, subject.gr = mask.gr, index = 'MASK')
      segm <- segm[tmp$PercOverlap_MASK <= mask.overlap.th]
    }
  }
  
  ## Merge consecutive haploype segments
  if (length(segm) > 1) {
    segm <- primatR::collapseBins(gr = segm, id.field = 6, measure.field = 2)
  } else {
    segm <- segm
  }
  
  ## Make sure first and last haplotype segment reaches start and end of the chromosome
  if (length(segm) > 0) {
    start(segm)[which.min(start(segm))] <- 1
    if (!is.na(seqlengths(segm))) {
      end(segm)[which.max(end(segm))] <- seqlengths(segm)
    }
  }
  
  ## Get meiotic recombination breakpoints
  if (length(segm) > 1) {
    # suppressWarnings( recomb.break <- gaps(segm, start = start(segm)) )
    # start(recomb.break) <- start(recomb.break) - 1
    # end(recomb.break) <- end(recomb.break) + 1
    
    ## Make sure segments are not overlapping [Report message!!!]
    shift.n <- c(pmin(0, start(segm)[-1] - end(segm)[-length(segm)]), 0)
    if (!all(shift.n == 0)) {
      message('    Adjusting overlapping genomic segments to continous scale!!!')
      segm <- resize(segm, width = width(segm) + shift.n, fix = 'start')
    }
    
    recomb.break <- GenomicRanges::GRanges(seqnames = unique(seqnames(segm)),
                                           ranges = IRanges(start = end(segm)[-length(segm)],
                                                            end = start(segm)[-1])
    )
  } else {
    ## Create a dummy breakpoint
    #recomb.break <- GenomicRanges::GRanges(seqnames = 'none', ranges = IRanges(start=1, end=1))
    ## Report empty GRanges object
    recomb.break <- GenomicRanges::GRanges()
  }
  return(list(hap.segm = segm, recomb.break = recomb.break, fastseg.minSeg = minSeg))
}