#' circle chart showed the results of GO and KEGG analysis
#'
#' @param data R clusterprofiler package for KEGG and GO results
#' @param top According to the order of p adjust value from small to large
#' the number of categories to show
#' @param label.name "ID" or "Description"
#' @param root root name
#' @param color.node node color
#' @param color.alpha color alpha
#' @param text.size text size
#' (1) first circle text size
#' (2) Second circle text size
#' @param ... additional parameters
#'
#' @return figure
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import ggraph
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr separate_rows
#' @importFrom ccgraph gather_graph_node
#' @importFrom ccgraph gather_graph_edge
#' @importFrom tidygraph tbl_graph
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
#' pathway_ccplot(KK)
#' }
#'
pathway_ccplot <- function(data,
                           top = 5,
                           label.name = "Description",
                           root = "KEGG",
                           color.node = "Paired",
                           color.alpha = 0.5,
                           text.size = c(3, 4),
                           ...) {
  path <- separate_rows(data@result[1:top, ], geneID, sep = "/")
  if (label.name == "ID") {
    kegg.df <- path %>%
      select(ID, geneID) %>%
      mutate(frequen = as.numeric("1"))
  } else if (label.name == "Description") {
    kegg.df <- path %>%
      select(Description, geneID) %>%
      mutate(frequen = as.numeric("1"))
  } else {
    print("The label.name is 'ID' or 'Description'. ")
  }

  # use the selected terms to build the data format
  se_index <- c(ifelse(label.name == "ID", "ID", "Description"), "geneID")
  nodes_kegg <- ccgraph::gather_graph_node(kegg.df,
    index = se_index,
    value = "frequen",
    root = root
  )
  edges_kegg <- ccgraph::gather_graph_edge(kegg.df,
    index = se_index,
    root = root
  )
  graph <- tbl_graph(nodes_kegg, edges_kegg)
  ggraph(graph, layout = "dendrogram", circular = TRUE) +
    geom_edge_diagonal(aes(color = node1.node.branch),
      alpha = color.alpha
    ) +
    geom_node_point(aes(size = node.size, color = node.branch),
      alpha = color.alpha
    ) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none") +
    scale_size(range = c(0.5, 30)) +
    geom_node_text(
      aes(
        x = 1.0175 * x,
        y = 1.0175 * y,
        label = node.short_name,
        angle = -((-node_angle(x, y) + 90) %% 180) + 90,
        filter = leaf, color = node.branch
      ),
      size = text.size[1], hjust = "outward"
    ) +
    scale_colour_manual(values = rep(
      RColorBrewer::brewer.pal(8, color.node),
      ifelse(label.name == "ID",
        length(unique(path$ID)),
        length(unique(path$"Description"))
      )
    )) +
    # add inner circle text label
    geom_node_text(
      show.legend = FALSE,
      aes(
        label = node.short_name,
        filter = !leaf,
        color = node.branch
      ),
      fontface = "bold",
      size = text.size[2]
    ) +
    coord_cartesian(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2))
}
