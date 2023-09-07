#' Chord diagram display of KEGG or GO results
#'
#' @param data data.frame
#' @param top number of pathways
#' @param id.color The color of the first layer ID of the path
#' @param id.cex The text size of the first layer ID of the path
#' @param gene.color The color of the number of all genes in this pathway
#' @param gene.cex The font size of the number of all genes in this pathway
#' @param ratio.color
#' （1）The color of the number of genes enriched this time on this pathway
#' （2）Other genes involved in this enrichment analysis
#' @param ratio.cex
#' The text size of the numbers involved in this enrichment analysis
#' @param p.color padjust color
#' @param p.cex The text size of padjust
#' @param richFactor.color richFactor color
#' @param legend.position legend position
#' @param text.width legend text width
#' @param ... additional parameters
#'
#' @return chord diagram
#' @export
#'
#' @import circlize
#' @import dplyr
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stringr str_split
#' @importFrom ComplexHeatmap Legend
#' @importFrom ComplexHeatmap packLegend
#' @importFrom stringr str_wrap
#' @importFrom grid gpar
#' @importFrom grid pushViewport
#' @importFrom grid grid.draw
#' @importFrom grid viewport
#'
#' @examples
#' \dontrun{
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#' library(DOSE)
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
#' cir_plot(KK)
#' }
cir_plot <- function(data,
                     top = 20,
                     id.color = "RdBu",
                     id.cex = 0.8,
                     gene.color = "#CFB0D4",
                     gene.cex = 1,
                     ratio.color = c("#6291bd", "#B1D1E7"),
                     ratio.cex = 1,
                     p.color = c("#eda9aa", "#c25254"),
                     p.cex = 1,
                     richFactor.color = "RdBu",
                     legend.position = c(0.8, 0.45),
                     text.width = 35,
                     ...) {
  # data processing
  enrich.res <- data@result
  enrich.res$left <- 0
  enrich.res$this_pathway_gene <- enrich.res$BgRatio %>%
    sapply(function(x) {
      str_split(x, "/")[[1]][1]
    }) %>%
    as.numeric()
  enrich.res$all_pathway_gene <- enrich.res$BgRatio %>%
    sapply(function(x) {
      str_split(x, "/")[[1]][2]
    }) %>%
    as.numeric()
  enrich.res$DEGnum <- enrich.res$GeneRatio %>%
    sapply(function(x) {
      str_split(x, "/")[[1]][2]
    }) %>%
    as.numeric()
  enrich.res$right <- max(enrich.res$this_pathway_gene)
  enrich.res <- mutate(enrich.res,
    richFactor = Count / as.numeric(sub(
      "/\\d+",
      "", BgRatio
    ))
  )
  enrich.res <- enrich.res %>% arrange(desc(richFactor))
  rownames(enrich.res) <- enrich.res$ID
  enrich.res <- enrich.res[1:top, ]

  # start drawing
  plotdata <- enrich.res
  circos.clear()

  circlize::circos.par(
    canvas.xlim = c(-0.2, 1.2),
    canvas.ylim = c(-1.0, 1.2),
    points.overflow.warning = FALSE,
    clock.wise = TRUE,
    start.degree = 90,
    gap.degree = 0.8,
    xaxis.clock.wise = TRUE
  )

  # The first circle: classification information
  col <- brewer.pal(8, name = id.color)
  ID_col <- colorRampPalette(col)(length(unique(enrich.res$ID)))
  names(ID_col) <- unique(enrich.res$ID)
  circlize::circos.genomicInitialize(plotdata[, c("ID", "left", "right")],
    plotType = NULL
  )
  circlize::circos.track(
    ylim = c(0, 1), track.height = 0.08, bg.border = NA, bg.col = ID_col,
    panel.fun = function(x, y) {
      ylim <- get.cell.meta.data("ycenter")
      xlim <- get.cell.meta.data("xcenter")
      sector.name <- get.cell.meta.data("sector.index")
      circos.text(xlim, ylim, sector.name, cex = id.cex, niceFacing = TRUE)
    }
  )

  # The second circle: how many genes are there in total in this pathway
  circlize::circos.genomicTrackPlotRegion(
    plotdata[, c("ID", "left", "this_pathway_gene")],
    track.height = 0.08, bg.border = NA,
    stack = TRUE,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value,
        col = gene.color,
        border = NA, ...
      )
      ylim <- get.cell.meta.data("ycenter")
      xlim <- plotdata[, c(
        "ID",
        "left",
        "this_pathway_gene"
      )][get.cell.meta.data("sector.index"), 3] / 2
      sector.name <- plotdata[, c(
        "ID",
        "left",
        "this_pathway_gene"
      )][get.cell.meta.data("sector.index"), 3]
      circos.text(xlim, ylim,
        sector.name,
        cex = gene.cex,
        niceFacing = TRUE
      )
    }
  )

  # The third circle: generatio related
  plotdata3 <- plotdata[, c("ID", "left", "Count", "DEGnum", "right")]
  plotdata3$ratio <- plotdata3$Count / plotdata3$DEGnum
  plotdata3$len <- plotdata3$ratio * plotdata3$right
  plotdata3$len2 <- plotdata3$right - plotdata3$len

  tmpdf1 <- plotdata3[, c("ID", "left", "len")]
  colnames(tmpdf1) <- c("ID", "start", "end")
  tmpdf1$type <- 1
  tmpdf2 <- plotdata3[, c("ID", "len", "right")]
  colnames(tmpdf2) <- c("ID", "start", "end")
  tmpdf2$type <- 2
  tmpdata <- tmpdf1 %>% rbind(tmpdf2)
  color_assign <- colorRamp2(breaks = c(1, 2), colors = ratio.color)
  circlize::circos.genomicTrackPlotRegion(
    tmpdata,
    track.height = 0.08, bg.border = NA, stack = TRUE,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value,
                         col = color_assign(value),
                         border = NA, ...)

      ylim <- get.cell.meta.data("ycenter")
      xlim <- plotdata3[get.cell.meta.data("sector.index"), "len"] / 2
      sector.name <- plotdata3[get.cell.meta.data("sector.index"), "Count"]
      circos.text(xlim, ylim, sector.name, cex = 1.0, niceFacing = TRUE)

      xlim <- plotdata3[get.cell.meta.data("sector.index"), "len"] +
        plotdata3[get.cell.meta.data("sector.index"), "len2"] / 2
      sector.name <- plotdata3[get.cell.meta.data("sector.index"), "DEGnum"] -
        plotdata3[get.cell.meta.data("sector.index"), "Count"]
      circos.text(xlim, ylim, sector.name,
        cex = ratio.cex, niceFacing = TRUE
      )
    }
  )

  # The fourth circle: p.adjust value
  plotdata4 <- plotdata[, c("ID", "left", "right", "p.adjust")]
  total.len <- unique(plotdata4$right)
  plotdata4$p.adjust_neg <- -log10(plotdata4$p.adjust)
  plotdata4$relative_value <-
    plotdata4$p.adjust_neg / max(plotdata4$p.adjust_neg) * total.len

  p_max <- max(plotdata4$p.adjust_neg) %>% ceiling()
  colorsChoice <- grDevices::colorRampPalette(p.color)
  color_assign <- colorRamp2(breaks = 0:p_max,
                             colors = colorsChoice(p_max + 1))

  plotdata4 <- plotdata4[, c(
    "ID",
    "left",
    "relative_value",
    "p.adjust_neg"
  )]

  circlize::circos.genomicTrackPlotRegion(
    plotdata4,
    track.height = 0.08, bg.border = NA, stack = TRUE,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value,
        col = color_assign(value),
        border = NA, ...
      )

      ylim <- get.cell.meta.data("ycenter")
      xlim <- plotdata4[
        get.cell.meta.data("sector.index"),
        "relative_value"
      ] / 2
      sector.name <- plotdata4[
        get.cell.meta.data("sector.index"),
        "p.adjust_neg"
      ] %>%
        round(2)

      tmpx <- -log10(0.01) / max(plotdata4$p.adjust_neg) * total.len
      circos.lines(c(tmpx, tmpx),
        c(get.cell.meta.data("ylim")[1], get.cell.meta.data("ylim")[2]),
        col = "gray", lwd = 0.45
      )

      circos.text(xlim, ylim, sector.name, cex = p.cex, niceFacing = TRUE)
    }
  )

  # Fifth Circle: Enrichment Factors
  plotdata5 <- plotdata[, c(
    "ID", "left", "right",
    "richFactor", "Description"
  )]
  label_data <- plotdata5["Description"]
  col <- RColorBrewer::brewer.pal(8, name = richFactor.color)
  color_assign <-
    grDevices::colorRampPalette(col)(length(unique(enrich.res$Description)))

  names(color_assign) <- plotdata5$Description
  circlize::circos.genomicTrack(
    plotdata5,
    ylim = c(0, max(plotdata5$richFactor)),
    track.height = max(plotdata5$richFactor),
    bg.col = "gray95",
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      sector.name <- get.cell.meta.data("sector.index")
      circos.genomicRect(region, value,
        col = color_assign[label_data[sector.name, 1]],
        border = NA, ytop.column = 1, ybottom = 0, ...
      )
      circos.lines(c(0, max(region)), c(0.5, 0.5), col = "gray", lwd = 0.3)
    }
  )


  # draw legend
  label.text <- as.character(unique(paste0(
    plotdata$ID, ":",
    plotdata$Description
  )))
  category_legend <- ComplexHeatmap::Legend(
    labels = stringr::str_wrap(label.text, width = text.width),
    background = color_assign,
    type = "points",
    pch = NA,
    labels_gp = grid::gpar(fontsize = 8),
    grid_height = unit(0.5, "cm"),
    grid_width = unit(0.5, "cm")
  )

  thispathway_gene_legend <- ComplexHeatmap::Legend(
    labels = c("Gene number of all this pathway"),
    background = gene.color,
    type = "points", pch = NA,
    labels_gp = grid::gpar(fontsize = 8),
    grid_height = unit(0.5, "cm"),
    grid_width = unit(0.5, "cm")
  )

  generatio_legend <- ComplexHeatmap::Legend(
    labels = c("Genes in this pathway", "other Genes"),
    background = ratio.color,
    type = "points", pch = NA,
    labels_gp = grid::gpar(fontsize = 8),
    grid_height = unit(0.5, "cm"),
    grid_width = unit(0.5, "cm")
  )

  pvalue_legend <- ComplexHeatmap::Legend(
    col_fun = colorRamp2(
      breaks = 0:p_max,
      colors = grDevices::colorRampPalette(p.color)(p_max + 1)
    ),
    legend_height = unit(3, "cm"),
    labels_gp = grid::gpar(fontsize = 8),
    title = "-log10(p.adjust)",
    title_gp = grid::gpar(fontsize = 9),
    title_position = "lefttop",
    direction = "horizontal"
  )

  pack_Legend <- ComplexHeatmap::packLegend(
    category_legend,
    thispathway_gene_legend,
    generatio_legend,
    pvalue_legend,
    row_gap = unit(0.2, "cm")
  )
  grid::pushViewport(grid::viewport(
    x = legend.position[1],
    y = legend.position[2]
  ))
  grid::grid.draw(pack_Legend)
  circos.clear()
}
