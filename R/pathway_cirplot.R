#' Display the relationship between genes and pathways
#' enriched by KEGG or GO in the form of a chord diagram
#'
#' @param data data.frame
#' @param top According to the order of p adjust value from small to large
#' the number of categories to show
#' @param start.degree the angle at which the circle plot starts
#' @param label.name "ID" or "Description"
#' @param text.width label length
#' @param text.size the font size of each layer circle
#' @param color color see "RColorBrewer::display.brewer.all()"
#' @param ... additional parameters
#'
#' @return chord diagram
#' @export
#'
#' @import circlize
#' @import dplyr
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom cowplot ggdraw
#' @importFrom tidyr separate_rows
#'
#' @examples
#' \dontrun{
#' data(xfbdf, package = "TCMR")
#' library(org.Hs.eg.db)
#' library(clusterProfiler)
#' library(DOSE)
#' library(tidyverse)
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
#' pathway_cirplot(KK)
#' }
pathway_cirplot <- function(data,
                            top = 5,
                            start.degree = -90,
                            label.name = "Description",
                            text.width = 15,
                            text.size = 0.6,
                            color = "RdBu", ...) {
  df <- separate_rows(data@result[1:top, ], geneID, sep = "/") %>%
    mutate(col = "grey")
  gene_col <- df$col %>% head(length(unique(df$geneID)))
  names(gene_col) <- unique(df$geneID)
  circlize::circos.par(
    canvas.xlim = c(-1, 1),
    canvas.ylim = c(-1, 1),
    start.degree = start.degree
  )
  if (label.name == "ID") {
    df <- df %>% dplyr::select(ID, geneID)
    cols_number <- length(unique(df$ID))
    cols <- brewer.pal(ifelse(cols_number > 8, 8, cols_number), color)
    pathway_col <- colorRampPalette(cols)(length(unique(df$ID)))
    names(pathway_col) <- unique(df$ID)
  } else if (label.name == "Description") {
    df <- df %>% dplyr::select(Description, geneID)
    df$Description <- lapply(
      strwrap(df$Description,
        width = text.width,
        simplify = FALSE
      ),
      paste,
      collapse = "\n"
    )
    cols_number <- length(unique(df$Description))
    cols <- brewer.pal(ifelse(cols_number > 8, 8, cols_number), color)
    pathway_col <- colorRampPalette(cols)(length(unique(df$Description)))
    names(pathway_col) <- unique(df$Description)
  } else {
    print("The label.name is ID' or 'Description'. ")
  }
  grid.col <- c(pathway_col, gene_col)
  set.seed(1234)
  chordDiagram(df,
    grid.col = grid.col,
    link.decreasing = TRUE,
    transparency = 0.1,
    big.gap = 1,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = .1)
  )

  for (si in get.all.sector.index()) {
    xlim <- get.cell.meta.data("xlim",
      sector.index = si,
      track.index = 1
    )
    ylim <- get.cell.meta.data("ylim",
      sector.index = si,
      track.index = 1
    )
    circos.text(mean(xlim),
      ylim[1],
      labels = si,
      sector.index = si,
      track.index = 1,
      facing = "clockwise",
      cex = 0.6,
      adj = c(0, 0.5),
      niceFacing = TRUE
    )
  }
  circos.clear()
}
