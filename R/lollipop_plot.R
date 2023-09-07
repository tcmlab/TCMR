#' Lollipop chart display of KEGG/GO results
#'
#' @param data R clusterprofiler package for KEGG and GO results
#' @param top According to the order of p adjust value from small to large
#' the number of categories to show
#' @param color color color see "RColorBrewer::display.brewer.all()"
#' @param title title
#' @param text.size text size
#' @param text.width text width
#' @param ... additional parameters
#'
#' @return lollipop plot
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stringr str_wrap
#' @importFrom stats reorder
#'
#' @examples
#' \dontrun{
#' data(xfbdf)
#' library(org.Hs.eg.db)
#' library(clusterProfiler)
#' eg <- bitr(unique(xfbdf$target),
#'   fromType = "SYMBOL",
#'   toType = "ENTREZID",
#'   OrgDb = "org.Hs.eg.db"
#' )
#' KK <- enrichKEGG(
#'   gene = eg$ENTREZID,
#'   organism = "hsa",
#'   pvalueCutoff = 0.05
#' )
#' lollipop_plot(KK)
#' }
lollipop_plot <- function(data,
                          top = 15,
                          color = "RdBu",
                          title = NULL,
                          text.size = 10,
                          text.width = 35, ...) {
  # data processing
  data2 <- mutate(data,
    richFactor = Count / as.numeric(sub(
      "/\\d+",
      "", BgRatio
    ))
  )

  # ggplot2 plotting
  p <- ggplot(data2,
    showCategory = top,
    aes(richFactor, stats::reorder(Description, richFactor))
  ) +
    geom_segment(aes(xend = 0, yend = Description)) +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_gradientn(
      colours = RColorBrewer::brewer.pal(8, color),
      trans = "log10",
      guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    scale_size_continuous(range = c(2, 10)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        size = text.size,
        colour = "black",
        vjust = 1
      ),
      axis.text.y = element_text(
        size = text.size,
        colour = "black",
        hjust = 1
      )
    ) +
    theme(axis.title = element_text(
      margin = margin(10, 5, 0, 0),
      color = "black",
      size = text.size
    )) +
    xlab("Rich Factor") +
    ylab(NULL) +
    ggtitle(title) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = text.width))
  return(p)
}
