#' Bubble chart display of KEGG/GO results
#'
#' @param data R clusterprofiler package for KEGG or GO results
#' @param padjust such as padjust = 0.01, mean that p.adjust <= 0.01
#' @param top According to the order of P adjust value from small to large
#' the number of categories to show
#' @param text.size text size
#' @param color color see "RColorBrewer::display.brewer.all()"
#' @param ... additional parameters
#'
#' @return bubble plot
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stats reorder
#' @importFrom ggrepel geom_label_repel
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' \dontrun{
#' data(xfbdf)
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
#' bubble_plot(KK)
#' }
bubble_plot <- function(data,
                        padjust = 0.01,
                        top = 30,
                        text.size = 4,
                        color = "RdBu", ...) {
  # data processing
  kk.se <- subset(data@result, p.adjust <= padjust)
  kk.se$logP <- -log10(kk.se$p.adjust)
  kk.se <- kk.se %>% arrange(desc(logP))
  col_bar2 <- colorRampPalette(brewer.pal(8, color))(top)
  names(col_bar2) <- kk.se[1:top, ]$Description

  # ggplot2 plotting
  p <- ggplot(kk.se[1:top, ]) +
    geom_point(aes(
      x = logP,
      y = Description,
      size = Count,
      color = Description
    )) +
    scale_color_manual(values = col_bar2) +
    theme_test() +
    guides(colour = "none") +
    ggrepel::geom_label_repel(
      aes(
        x = logP,
        y = Description,
        label = Description
      ),
      size = text.size,
      nudge_y = 0.1
    ) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    xlim(c(min(kk.se$logP) - 2, max(kk.se$logP) + 2)) +
    xlab(bquote(-Log[10] ~ italic("Padjust"))) +
    ylab(bquote("Pathway"))
  return(p)
}
