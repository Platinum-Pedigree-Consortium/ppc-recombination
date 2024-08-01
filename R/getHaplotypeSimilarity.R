#' Calculate binned similarity between haplotypes.
#'
#' This function takes a single haplotype of a child and compares it to both parental haplotypes from whom
#' the child's haplotype was inherited from.
#'
#' @param df A \code{data.frame} with 5 columns; `start` and `end` position of a variant (allele), child's haplotype (`child`) and
#' two parental haplotypes (`p1` and `p2`) from whom the child's haplotype was inherited from.
#'
#' @return A \code{matrix}
#' @author David Porubsky
#' @export
#'
getHaplotypeSimilarity <- function(df) {
  start <- as.numeric(min(df[,1]))
  end <- as.numeric(max(df[,2]))
  child <- df[,3]
  p1 <- df[,4]
  p2 <- df[,5]
  agree.hap1 <- length(child[child == p1])
  agree.hap1.perc <- as.numeric(agree.hap1/length(child))
  agree.hap2 <- length(child[child == p2])
  agree.hap2.perc <- as.numeric(agree.hap2/length(child))
  
  df.new <- data.frame(start=start,
                       end=end,
                       p1.simil=agree.hap1.perc,
                       p2.simil=agree.hap2.perc)
  return(df.new)
}