#' Function to find a changepoint where child haplotype switches matching allele 1 to allele 2 (or vice versa) of inherited parental haplotype.
#'
#' @param rec.mat A two-row matrix composed of 1's and 2's assigning allele 1 or allele 2 match of child haplotype to parental haplotype.
#' @param bin.width A a size of a sliding to window to scan for changepoints in 'rec.mat' vectors.
#' @return A \code{data.frame} containing start and end position of a changepoint between allele 1 and 2. 
#' @author David Porubsky
#'
getAlleleChangepoints <- function(rec.mat=NULL, bin.width = 50) {
  ## Get binned deltas per matrix row
  deltas.l <- list()
  for (i in 1:nrow(rec.mat)) {
    ID <- rownames(rec.mat)[i]
    #num.vector <- as.numeric(rec.mat[i,])
    num.vector <- factor(rec.mat[i,], levels = c('1', '2', '0'))
    #slide through the numeric vector and count the number of occurrences of each allele in a given window
    counts <- zoo::rollapply(data = as.factor(num.vector), width=bin.width, by=1, FUN=function(x) table(x, useNA='always'))
    #get allele counts for left and right window
    counts.left <- counts[1:(nrow(counts) - bin.width),]
    counts.right <- counts[(bin.width + 1):nrow(counts),]
    #calculate difference between counts for allele 1 in the left and right window [perhaps allele 1 is not always informative!!!] 
    diff <- as.numeric( abs(counts.left[,'Freq1'] - counts.right[,'Freq1']) )
    #get indices of window boundaries for left and right window
    left.bin.idx <- seq(bin.width, length(num.vector) - bin.width, by=1)
    right.bin.idx <- seq(bin.width + 1, length(num.vector) - (bin.width - 1), by=1)
    #translate window boundaries into genomic location
    left.bin.genCoord <- as.numeric(names(num.vector)[left.bin.idx])
    right.bin.genCoord <- as.numeric(names(num.vector)[right.bin.idx])
    #export data into data.frame and store in a list
    diff.df <- data.frame(ID=ID, start.idx=left.bin.idx, end.idx=right.bin.idx, start=left.bin.genCoord, end=right.bin.genCoord, delta=diff)
    deltas.l[[i]] <- diff.df
  }
  #deltas.all <- do.call(rbind, deltas.l)
  ## Sum up deltas
  delta.sum <- deltas.l[[1]]$delta + deltas.l[[2]]$delta
  
  ## Export position of the max peak
  max.peak <- deltas.l[[1]][which.max(delta.sum), c('start.idx','end.idx','start','end')]
  
  return(max.peak)
}  
