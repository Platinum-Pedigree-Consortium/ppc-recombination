#' Map meiotic recombination events in a single family trio [father-mother-child].
#'
#' This function takes as an input phased VCF files for each member of a family trio and maps breakpoints
#' of meiotic recombination in each inherited parental homolog.
#'
#' @param parent1 A path to a VCF file to be loaded for a parent 1.
#' @param parent2 A path to a VCF file to be loaded for a parent 2.
#' @param child A path to a VCF file to be loaded for a child.
#' @param method A user defined method to be used to map changes in haplotype blocks [default: CBS]
#' @param male.X Set to TRUE if processing chromosome X from male child [default: FALSE]
#' @inheritParams findRecomb
#' @inheritParams vcf2vranges
#' @importFrom fastseg fastseg
#' @importFrom dplyr recode
#' @return A \code{list} object that contains mapped meiotic breakpoints and inherited hapltype segments
#' for each homolog in a child.
#' @author David Porubsky
#' @export
#'
mapRecomb <- function(parent1=NULL, parent2=NULL, child=NULL, genome='hg38', method='CBS', minSeg=NULL, minSegBp=100000, smooth=3, collapse.amb=TRUE, mask.gr=NULL, mask.overlap.th = 50, male.X=FALSE) {
  
  ## Helper function
  switchValue <- function(x) {
    if (x == 1) {
      x <- 0
    } else {
      x <- 1
    }
  }
  
  ## Load VCF files
  if (!is.null(parent1) & file.exists(parent1)) {
    par1 <- vcf2vranges(vcfFile = parent1, genoField = 'GT', translateBases = TRUE, genome = genome)
    par1 <- keepSeqlevels(par1, value = unique(seqnames(par1)), pruning.mode = 'coarse')
  } else {
    ## Create a dummy VCF record
    par1 <- VRanges(seqnames = 'chr', ranges = IRanges(start=1, end=1),
                    ref = 'A', alt = 'N',
                    totalDepth = 0, refDepth = 0, altDepth = 0,
                    sampleNames = 'parent1', softFilterMatrix = matrix(FALSE))
  }
  if (!is.null(parent2) & file.exists(parent2)) {
    par2 <- vcf2vranges(vcfFile = parent2, genoField = 'GT', translateBases = TRUE, genome = genome)
    par2 <- keepSeqlevels(par2, value = unique(seqnames(par2)), pruning.mode = 'coarse')
  } else {
    ## Create a dummy VCF record
    par2 <- VRanges(seqnames = 'chr', ranges = IRanges(start=1, end=1),
                    ref = 'A', alt = 'N',
                    totalDepth = 0, refDepth = 0, altDepth = 0,
                    sampleNames = 'parent2', softFilterMatrix = matrix(FALSE))
  }
  if (!is.null(child) & file.exists(child)) {
    child <- vcf2vranges(vcfFile = child, genoField = 'GT', translateBases = TRUE, genome = genome)
    child <- keepSeqlevels(child, value = unique(seqnames(child)), pruning.mode = 'coarse')
  }
  ## Keep only SNVs
  #child <- child[isSNV(child)]
  ## Filter out variants from child with identical start position
  #dup.vars <- start(child)[duplicated(start(child))]
  #child <- child[!start(child) %in% dup.vars]
  child <- child[!duplicated(start(child))]
  
  ## Remove child's missing values (Expects diploid chromosomes)
  ## Switch off chromosome X in males that carry only one copy !!!
  if (male.X == FALSE) {
    mask <- child$H1 == 'N' | child$H2 == 'N'
    child <- child[!mask]
  }  
  ## Keep only child's HET SNVs
  mask <- child$H1 == child$H2
  child <- child[!mask]
  
  ## Get only shared SNVs between datasets
  comparison.obj <- child
  comparison.obj$par1.H1 <- 'N'
  comparison.obj$par1.H2 <- 'N'
  comparison.obj$par2.H1 <- 'N'
  comparison.obj$par2.H2 <- 'N'
  if (!is.null(parent1) & file.exists(parent1)) {
    shared.par1 <- findOverlaps(par1, child)
    comparison.obj$par1.H1[subjectHits(shared.par1)] <- par1$H1[queryHits(shared.par1)]
    comparison.obj$par1.H2[subjectHits(shared.par1)] <- par1$H2[queryHits(shared.par1)]
  }
  if (!is.null(parent2) & file.exists(parent2)) {
    shared.par2 <- findOverlaps(par2, child)
    comparison.obj$par2.H1[subjectHits(shared.par2)] <- par2$H1[queryHits(shared.par2)]
    comparison.obj$par2.H2[subjectHits(shared.par2)] <- par2$H2[queryHits(shared.par2)]
  }
  
  ## Keep parental HOM variants only to be able to assign variants inherited from a given parent
  ## Compare child.H1 to par1
  mask <- comparison.obj$H1 != 'N' & comparison.obj$par1.H1 != 'N' & comparison.obj$par1.H2 != 'N' & comparison.obj$par1.H1 == comparison.obj$par1.H2
  c1_to_par1.match <- length(comparison.obj$H1[mask][comparison.obj$H1[mask] == comparison.obj$par1.H1[mask]])
  ## Compare child.H2 to par1
  mask <- comparison.obj$H2 != 'N' & comparison.obj$par1.H1 != 'N' & comparison.obj$par1.H2 != 'N' & comparison.obj$par1.H1 == comparison.obj$par1.H2
  c2_to_par1.match <- length(comparison.obj$H2[mask][comparison.obj$H2[mask] == comparison.obj$par1.H1[mask]])
  ## Compare child.H1 to par2
  mask <- comparison.obj$H1 != 'N' & comparison.obj$par2.H1 != 'N' & comparison.obj$par2.H2 != 'N' & comparison.obj$par2.H1 == comparison.obj$par2.H2
  c1_to_par2.match <- length(comparison.obj$H1[mask][comparison.obj$H1[mask] == comparison.obj$par2.H1[mask]])
  ## Compare child.H2 to par2
  mask <- comparison.obj$H2 != 'N' & comparison.obj$par2.H1 != 'N' & comparison.obj$par2.H2 != 'N' & comparison.obj$par2.H1 == comparison.obj$par2.H2
  c2_to_par2.match <- length(comparison.obj$H2[mask][comparison.obj$H2[mask] == comparison.obj$par2.H1[mask]])
  
  ## Get parental homologs ##
  if ((c1_to_par1.match + c2_to_par2.match) > (c2_to_par1.match + c1_to_par2.match)) {
    ## Get comparison of c1 to par2
    comparison.obj$c1_to_par1 <- 'N'
    ## Filter only non-missing alleles and those that are bi-allelic
    mask <- comparison.obj$H1 != 'N' & comparison.obj$par1.H1 != 'N' & comparison.obj$par1.H2 != 'N' & comparison.obj$par1.H1 != comparison.obj$par1.H2 &
      (comparison.obj$H1 == comparison.obj$par1.H1 | comparison.obj$H1 == comparison.obj$par1.H2 | comparison.obj$H1 == comparison.obj$par2.H1 | comparison.obj$H1 == comparison.obj$par2.H2) &
      (comparison.obj$H2 == comparison.obj$par1.H1 | comparison.obj$H2 == comparison.obj$par1.H2 | comparison.obj$H2 == comparison.obj$par2.H1 | comparison.obj$H2 == comparison.obj$par2.H2)
    comparison.obj$c1_to_par1[mask] <- ifelse(comparison.obj$H1[mask] == comparison.obj$par1.H1[mask], 'P1', 'P2')
    ## Get comparison of c2 to par1
    comparison.obj$c2_to_par2 <- 'N'
    ## Filter only non-missing alleles and those that are bi-allelic
    mask <- comparison.obj$H2 != 'N' & comparison.obj$par2.H1 != 'N' & comparison.obj$par2.H2 != 'N' & comparison.obj$par2.H1 != comparison.obj$par2.H2 &
      (comparison.obj$H1 == comparison.obj$par1.H1 | comparison.obj$H1 == comparison.obj$par1.H2 | comparison.obj$H1 == comparison.obj$par2.H1 | comparison.obj$H1 == comparison.obj$par2.H2) &
      (comparison.obj$H2 == comparison.obj$par1.H1 | comparison.obj$H2 == comparison.obj$par1.H2 | comparison.obj$H2 == comparison.obj$par2.H1 | comparison.obj$H2 == comparison.obj$par2.H2)
    comparison.obj$c2_to_par2[mask] <- ifelse(comparison.obj$H2[mask] == comparison.obj$par2.H1[mask], 'P1', 'P2')
    comparisons <- comparison.obj
    inherited.homologs <- c(paste0('H1.', runValue(sampleNames(child)), " <= ",  runValue(sampleNames(par1))),
                            paste0('H2.', runValue(sampleNames(child)), " <= ",  runValue(sampleNames(par2)))
    )
    H1.inherited <- as.character(runValue(sampleNames(par1)))
    H2.inherited <- as.character(runValue(sampleNames(par2)))
  } else {
    ## Get comparison of c1 to par2
    comparison.obj$c1_to_par2 <- 'N'
    ## Filter only non-missing alleles and those that are bi-allelic
    mask <- comparison.obj$H1 != 'N' & comparison.obj$par2.H1 != 'N' & comparison.obj$par2.H2 != 'N' & comparison.obj$par2.H1 != comparison.obj$par2.H2 &
      (comparison.obj$H1 == comparison.obj$par1.H1 | comparison.obj$H1 == comparison.obj$par1.H2 |  comparison.obj$H1 == comparison.obj$par2.H1 | comparison.obj$H1 == comparison.obj$par2.H2) &
      (comparison.obj$H2 == comparison.obj$par1.H1 | comparison.obj$H2 == comparison.obj$par1.H2 |  comparison.obj$H2 == comparison.obj$par2.H1 | comparison.obj$H2 == comparison.obj$par2.H2)
    comparison.obj$c1_to_par2[mask] <- ifelse(comparison.obj$H1[mask] == comparison.obj$par2.H1[mask], 'P1', 'P2')
    ## Get comparison of c2 to par1
    comparison.obj$c2_to_par1 <- 'N'
    ## Filter only non-missing alleles and those that are bi-allelic
    mask <- comparison.obj$H2 != 'N' & comparison.obj$par1.H1 != 'N' & comparison.obj$par1.H2 != 'N' & comparison.obj$par1.H1 != comparison.obj$par1.H2 &
      (comparison.obj$H1 == comparison.obj$par1.H1 | comparison.obj$H1 == comparison.obj$par1.H2 |  comparison.obj$H1 == comparison.obj$par2.H1 | comparison.obj$H1 == comparison.obj$par2.H2) &
      (comparison.obj$H2 == comparison.obj$par1.H1 | comparison.obj$H2 == comparison.obj$par1.H2 |  comparison.obj$H2 == comparison.obj$par2.H1 | comparison.obj$H2 == comparison.obj$par2.H2)
    comparison.obj$c2_to_par1[mask] <- ifelse(comparison.obj$H2[mask] == comparison.obj$par1.H1[mask], 'P1', 'P2')
    comparisons <- comparison.obj
    inherited.homologs <- c(paste0('H1.', runValue(sampleNames(child)), " <= ",  runValue(sampleNames(par2))),
                            paste0('H2.', runValue(sampleNames(child)), " <= ",  runValue(sampleNames(par1)))
    )
    H1.inherited <- as.character(runValue(sampleNames(par2)))
    H2.inherited <- as.character(runValue(sampleNames(par1)))
  }
  ## Keep only comparisons columns
  ncols <- length(mcols(comparisons))
  comparisons <- comparisons[,c(ncols-1, ncols)]
  names(mcols(comparisons)) <- c('H1.comp', 'H2.comp')
  
  if (method == 'CBS') {
    H1.recomb <- findRecomb(comparison = comparisons[,'H1.comp'], minSeg = minSeg, minSegBp=minSegBp, smooth = smooth, collapse.amb = collapse.amb, mask.gr = mask.gr, mask.overlap.th = mask.overlap.th)
    H2.recomb <- findRecomb(comparison = comparisons[,'H2.comp'], minSeg = minSeg, minSegBp=minSegBp, smooth = smooth, collapse.amb = collapse.amb, mask.gr = mask.gr, mask.overlap.th = mask.overlap.th)
  }
  
  ## Add field of inherited homologs
  if (length(H1.recomb$hap.segm) > 0) {
    H1.recomb$hap.segm$inherited <- H1.inherited
  }
  if (length(H1.recomb$recomb.break) > 0) {
    H1.recomb$recomb.break$inherited <- H1.inherited
  }
  if (length(H2.recomb$hap.segm) > 0) {
    H2.recomb$hap.segm$inherited <- H2.inherited
  }
  if (length(H2.recomb$recomb.break) > 0) {
    H2.recomb$recomb.break$inherited <- H2.inherited
  }
  
  ## Return results object
  return(list(comparisons = comparisons,
              H1.segments = H1.recomb$hap.segm,
              H1.recomb.breaks = H1.recomb$recomb.break,
              H2.segments = H2.recomb$hap.segm,
              H2.recomb.breaks = H2.recomb$recomb.break,
              inherited.homologs = inherited.homologs)
  )
}