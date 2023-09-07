#' Bar graphs showed the results of GO and KEGG analysis
#'
#' @param data R clusterprofiler package for KEGG and GO results
#' @param top according to the order of p adjust value from small to large
#' the number of categories to show
#' @param color color see "RColorBrewer::display.brewer.all()"
#' @param text.size text size
#' @param text.width y-axis label length
#' @param text.bar The size of the text at the top of the bar graph
#' @param title  title
#' Kyoto Encyclopedia of Genes and Genomes, KEGG
#' cellular component, CC
#' molecular function, MF
#' biological process, BP
#' @param ... additional parameters
#'
#' @return barplot
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stats reorder
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stringr str_wrap
#'
#' @examples
#' \dontrun{
#' data(xfbdf, package = "TCMR")
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
#' KK <- setReadable(KK, "org.Hs.eg.db", keyType = "ENTREZID")
#' BP <- enrichGO(
#'   gene = eg$ENTREZID,
#'   "org.Hs.eg.db",
#'   ont = "BP",
#'   pvalueCutoff = 0.05,
#'   readable = TRUE
#' )
#' bar_plot(KK,title = "KEGG")
#' bar_plot(BP,title = "biological process")
#' }
bar_plot <- function(data,
                     top = 10,
                     color = "RdBu",
                     text.size = 10,
                     text.width = 35,
                     text.bar = 4,
                     title = NULL, ...) {
  # data processing
  data2 <- mutate(data,
    richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))
  )
  # ggplot2 plotting
  p <- ggplot(data2, showCategory = top) +
    geom_col(aes(x = Count, y = reorder(Description, Count), fill = p.adjust),
      color = "black",
      width = 0.6
    ) +
    geom_text(aes(x = Count, y = reorder(Description, Count), label = Count),
      hjust = -0.3,
      vjust = 0.5,
      size = text.bar
    ) + # Count value
    scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(8, color),
      trans = "log10",
      guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5,
        size = text.size,
        colour = "black"
      ),
      axis.text.y = element_text(
        angle = 0,
        size = text.size,
        face = "plain",
        colour = "black"
      )
    ) +
    scale_y_discrete(expand = c(0, -1)) +
    theme(legend.position = "right") +
    ylab(NULL) +
    ggtitle(title) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = text.width))
  return(p)
}
