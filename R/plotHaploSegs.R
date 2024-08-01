#' Plot inherited haplotype segments as ideogram.
#'
#' This function takes as na input a \code{\link{GRanges-class}} object with genomic positions of each
#' inherited haplotype block and project them as whole-genome ideogram.
#'
#' @param hapSegm.gr A \code{\link{GRanges-class}} object with regions of each inherited haplotype block
#' @param layout Set to 'vertical' or 'horizonal' based on desired layout of an resultant ideogram.
#' @param color.by A metadata columns in `hapSegm.gr` to be used to define color scheme.
#' @param color.palette A discrete color palette defined as named character vector (elements = colors, names = discrete levels).
#' @param mask.gr A \code{\link{GRanges-class}} object with regions that will plotted over the final ideogram.
#' @param mask.color A user defined color to be used for regions defined in `mask.gr`.
#' @importFrom dplyr "%>%"
#' @return A \code{ggplot} object.
#' @author David Porubsky
#' @export
#'
plotHaploSegs <- function(hapSegm.gr=NULL, layout='vertical', color.by=NULL, color.palette=NULL, mask.gr=NULL, mask.color=NULL) {
  ## Remove ranges with NA values
  hapSegm.gr <- hapSegm.gr[!is.na(hapSegm.gr$match)]
  ## Prepare data for plotting
  plt.df <- as.data.frame(hapSegm.gr)
  ## Set levels for homolog plotting
  homolog.ids <- unique(plt.df$inherited)
  plt.df$ymin <- 0
  plt.df$ymax <- 0
  plt.df$ymin[plt.df$inherited == homolog.ids[1]] <- 0
  plt.df$ymax[plt.df$inherited == homolog.ids[1]] <- 2
  plt.df$ymin[plt.df$inherited == homolog.ids[2]] <- 3
  plt.df$ymax[plt.df$inherited == homolog.ids[2]] <- 5
  
  ## Create unique ID
  if (!color.by %in% colnames(plt.df)) {
    plt.df$ID <- paste0(plt.df$match, ".", plt.df$inherited)
    color.by <- 'ID'
  }
  ## Define color palette
  if (is.null(color.palette)) {
    color.palette <- c("cadetblue3","coral", "dodgerblue4", "firebrick")
  }
  
  ## Set plotting themes
  theme_horizontal <- theme(legend.position ="top",
                            axis.text.y=element_blank(),
                            strip.text.y.left = element_text(angle = 0),
                            axis.ticks.y = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank())
  
  theme_vertical <- theme(legend.position ="top",
                          axis.line = element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          strip.text.y = element_text(angle = 180),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank())
  
  ## Get chromosome breaks and labels
  max.len <- signif(max(plt.df$end), digits = 2)
  breaks <- seq(from = 0, to = max.len, length.out = 6)
  labels <- breaks / 1000000
  labels <- paste0(labels, 'Mb')
  chr.num <- length(unique(plt.df$seqnames))
  
  ## Count number of breakpoints per homolog
  break.count <- plt.df %>%
    group_by(inherited, seqnames) %>%
    summarise(count=n() - 1) %>%
    group_by(inherited) %>%
    summarise(count=sum(count))
  ## Prepare breakpoint annotation
  break.count.annot <- paste0(break.count$inherited, ": ", break.count$count)
  #break.count.annot <- as.data.frame(break.count.annot)
  #break.count.annot$level <- 1:nrow(break.count.annot)
  break.count.annot <- paste0(break.count.annot, collapse = "\n")
  
  if (layout == 'horizontal') {
    plt <- ggplot(plt.df) +
      geom_rect(aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=.data[[color.by]])) +
      facet_grid(seqnames ~ ., switch = 'y') +
      scale_fill_manual(values=color.palette, name=break.count.annot) +
      scale_x_continuous(breaks = breaks, labels = labels, expand = c(0,0)) +
      theme_horizontal
    ## Mask user defined regions
    if (methods::is(mask.gr, 'GRanges')) {
      mask.df <- as.data.frame(mask.gr)
      if (is.null(mask.color)) {
        mask.color <- 'white'
      }
      plt <- plt + geom_rect(data=mask.df, aes(ymin=0, ymax=5, xmin=start, xmax=end), fill=mask.color)
    }
  } else if (layout == 'vertical') {
    plt <- ggplot(plt.df) +
      geom_rect(aes(xmin=ymin, xmax=ymax, ymin=start, ymax=end, fill=.data[[color.by]])) +
      facet_grid(. ~ seqnames, switch = 'x') +
      scale_fill_manual(values=color.palette, name=break.count.annot) +
      scale_y_continuous(breaks = breaks, labels = labels, expand = c(0,0)) +
      theme_vertical
    ## Mask user defined regions
    if (methods::is(mask.gr, 'GRanges')) {
      mask.df <- as.data.frame(mask.gr)
      if (is.null(mask.color)) {
        mask.color <- 'white'
      }
      plt <- plt + geom_rect(data=mask.df, aes(xmin=0, xmax=5, ymin=start, ymax=end), fill=mask.color, alpha=0.5)
    }
  } else {
    message("Unsupported layout!!! Please choose from 'horizontal' or 'vertical'.")
  }
  return(plt)
}