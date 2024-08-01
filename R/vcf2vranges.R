#' Load a VCF file
#'
#' This function loads raw VCF file into a \code{\link{VRanges-class}} object.
#'
#' @param vcfFile A path to a VCF file to be loaded.
#' @param genoField A vector of genotype IDs to be loaded from VCF [e.g. 'GT']
#' @param translateBases Set to \code{TRUE} if REF and ALT alleles should be reported as A,C,G or T.
#' @param genome A reference genome used by \code{readVcfAsVRanges} function. [e.g. 'hg38' - human]
#' @importFrom tidyr separate
#' @importFrom VariantAnnotation readVcfAsVRanges
#' @return A \code{\link{VRanges-class}} object.
#' @author David Porubsky
#' @export
#'
vcf2vranges <- function(vcfFile=NULL, genoField=NULL, translateBases=TRUE, genome='hg38') {
  ## Load vcf file into vranges object
  if (all(is.character(genoField) & nchar(genoField) > 0)) {
    suppressWarnings( vcf.vranges <- VariantAnnotation::readVcfAsVRanges(x = vcfFile, genome = genome,
                                                                         param = ScanVcfParam(fixed=c('ALT'), info = NA, geno = genoField)) )
  } else {
    suppressWarnings( vcf.vranges <- VariantAnnotation::readVcfAsVRanges(x = vcfFile, genome = genome,
                                                                         param = ScanVcfParam(fixed=c('ALT'), info = NA)) )
  }
  ## Split genotype field into hapltypes
  gen.field <-  tidyr::separate(data = as(mcols(vcf.vranges), 'data.frame'), col = GT, into = c("H1", "H2"))
  ## Translate 0/1 haplotypes into nucleotides
  if (translateBases) {
    allele1 <- rep(".", length(vcf.vranges))
    allele2 <- rep(".", length(vcf.vranges))
    allele1[gen.field$H1 == 0] <- ref(vcf.vranges)[gen.field$H1 == 0]
    allele1[gen.field$H1 == 1] <- alt(vcf.vranges)[gen.field$H1 == 1]
    allele1[allele1 == "."] <- 'N'
    allele2[gen.field$H2 == 0] <- ref(vcf.vranges)[gen.field$H2 == 0]
    allele2[gen.field$H2 == 1] <- alt(vcf.vranges)[gen.field$H2 == 1]
    allele2[allele2 == "."] <- 'N'
    gen.field$H1 <- allele1
    gen.field$H2 <- allele2
  }
  ## Construct final VRanges object
  mcols(vcf.vranges) <- gen.field
  return(vcf.vranges)
}