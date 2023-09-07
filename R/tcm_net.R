#' Network interaction diagram of herbs, ingredients, and targets
#'
#' @param network.data data frame
#' must contain herb, molecule, target three columns of data
#' @param node.color node color
#' see "RColorBrewer::display.brewer.all()"
#' @param node.size  node size
#' @param label.size label size
#' @param label.degree
#' the node degree is the number of connections that
#' the node has with the other nodes.
#' Nodes with connections greater than or
#' equal to degree will be displayed.
#' @param edge.color edge color
#' @param edge.width edge width
#' @param graph.layout etwork Diagram Layout:
#' "kk", nicely", "grid","sphere", "randomly",
#' "bipartite", "star","tree", nicely", "grid",
#' "sphere", "randomly","gem", "graphopt","lgl",
#' "mds", "sugiyama"
#' @param rem.dis.inter remove single free unconnected nodes
#' @param label.repel label repel
#'
#' @return Network Diagram
#' @export
#' @import dplyr
#' @import ggplot2
#' @import ggraph
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom igraph V
#' @importFrom igraph E
#' @importFrom igraph graph_from_data_frame
#' @examples
#' \dontrun{
#' data("xfbdf", package = "TCMR")
#' network.data <- xfbdf %>%
#'   dplyr::select(herb, molecule, target) %>%
#'   sample_n(100, replace = FALSE) %>%
#'   as.data.frame()
#' tcm_net(network.data, rem.dis.inter = TRUE)
#' }
tcm_net <- function(network.data,
                    node.color = "RdBu",
                    node.size = c(2, 8),
                    label.size = 4,
                    label.degree = 3,
                    label.repel = TRUE,
                    edge.color = "lightgrey",
                    edge.width = c(0.2, 2),
                    graph.layout = "kk",
                    rem.dis.inter = FALSE) {
  links <- rbind(
    network.data %>%
      dplyr::select(herb, molecule) %>%
      dplyr::rename(from = herb, to = molecule) %>%
      dplyr::mutate(weight = 1),
    network.data %>%
      dplyr::select(molecule, target) %>%
      dplyr::rename(from = molecule, to = target) %>%
      dplyr::mutate(weight = 1)
  ) %>%
    dplyr::distinct()
  # Create a network graph from links and nodes
  if (rem.dis.inter == FALSE && label.repel == TRUE) {
    nodes <- links %>%
      {
        data.frame(node = c(.$from, .$to))
      } %>%
      distinct()
    net <- igraph::graph_from_data_frame(
      d = links,
      vertices = nodes,
      directed = FALSE
    )
    # V and E are functions of the igraph package,
    # which are used to modify the nodes (nodes) and
    # connections (links) of the network graph respectively.
    igraph::V(net)$degree <- igraph::degree(net)
    igraph::V(net)$class <- c(
      rep(
        colnames(network.data)[1],
        length(intersect(network.data[, 1], nodes$node))
      ),
      rep(
        colnames(network.data)[2],
        length(intersect(network.data[, 2], nodes$node))
      ),
      rep(
        colnames(network.data)[3],
        length(intersect(network.data[, 3], nodes$node))
      )
    )
    igraph::V(net)$size <- igraph::degree(net)
    igraph::E(net)$score <- igraph::E(net)$weight

    # Plot with ggraph
    # ggraph is a package based on ggplot2,
    # the syntax is similar to regular ggplot2
    col <- RColorBrewer::brewer.pal(8, name = node.color)
    colorN <- grDevices::colorRampPalette(colors = col)(nrow(nodes))
    names(colorN) <- nodes$node

    ggraph(net, layout = graph.layout) +
      geom_edge_link0(aes(edge_linewidth = weight), edge_colour = edge.color) +
      geom_edge_fan(aes(edge_width = score, colour = score),
        color = edge.color, show.legend = TRUE
      ) +
      ggraph::geom_node_point(aes(
        color = degree,
        size = degree,
        shape = class
      ), alpha = 1.0) +
      scale_color_gradientn(colours = rev(colorN)) +
      geom_node_text(aes(filter = degree >= label.degree,
                         label = name), size = label.size, repel = TRUE) +
      scale_edge_width(range = edge.width) +
      scale_size_continuous(name = "degree", range = node.size) +
      theme_graph(base_family = "sans")
  } else if (rem.dis.inter == TRUE && label.repel == TRUE) {
    links <- links %>%
      mutate(from_c = count(., from)$n[match(from, count(., from)$from)]) %>%
      mutate(to_c = count(., to)$n[match(to, count(., to)$to)]) %>%
      filter(!(from_c == 1 & to_c == 1)) %>%
      dplyr::select(1, 2, 3)
    nodes <- links %>%
      {
        data.frame(node = c(.$from, .$to))
      } %>%
      distinct()
    col <- RColorBrewer::brewer.pal(8, name = node.color)
    colorN <- grDevices::colorRampPalette(colors = col)(nrow(nodes))
    names(colorN) <- nodes$node

    # Create a network graph from links and nodes
    net <- igraph::graph_from_data_frame(
      d = links,
      vertices = nodes,
      directed = FALSE
    )
    # V and E are functions of the igraph package,
    # which are used to modify the nodes (nodes) and
    # connections (links) of the network graph respectively.
    igraph::V(net)$degree <- igraph::degree(net)
    igraph::V(net)$class <- c(
      rep(
        colnames(network.data)[1],
        length(intersect(network.data[, 1], nodes$node))
      ),
      rep(
        colnames(network.data)[2],
        length(intersect(network.data[, 2], nodes$node))
      ),
      rep(
        colnames(network.data)[3],
        length(intersect(network.data[, 3], nodes$node))
      )
    )
    igraph::V(net)$size <- igraph::degree(net)
    igraph::E(net)$score <- igraph::E(net)$weight

    # Plot with ggraph
    ggraph(net, layout = graph.layout) +
      geom_edge_link0(aes(edge_linewidth = weight),
        edge_colour = edge.color
      ) +
      geom_edge_fan(aes(edge_width = score, colour = score),
        color = edge.color, show.legend = TRUE
      ) +
      ggraph::geom_node_point(aes(
        color = degree, size = degree,
        shape = class
      ), alpha = 1.0) +
      scale_color_gradientn(colours = rev(colorN)) +
      geom_node_text(
        aes(
          filter = degree >= label.degree,
          label = name
        ),
        size = label.size,
        repel = TRUE
      ) +
      scale_edge_width(range = edge.width) +
      scale_size_continuous(
        name = "degree",
        range = node.size
      ) +
      theme_graph(base_family = "sans")
  } else if (rem.dis.inter == FALSE && label.repel == FALSE) {
    nodes <- links %>%
      {
        data.frame(node = c(.$from, .$to))
      } %>%
      distinct()
    net <- igraph::graph_from_data_frame(
      d = links,
      vertices = nodes, directed = FALSE
    )
    # V and E are functions of the igraph package,
    # which are used to modify the nodes (nodes) and
    # connections (links) of the network graph respectively.
    igraph::V(net)$degree <- igraph::degree(net)
    igraph::V(net)$degree <- igraph::degree(net)
    igraph::V(net)$class <- c(
      rep(
        colnames(network.data)[1],
        length(intersect(network.data[, 1], nodes$node))
      ),
      rep(
        colnames(network.data)[2],
        length(intersect(network.data[, 2], nodes$node))
      ),
      rep(
        colnames(network.data)[3],
        length(intersect(network.data[, 3], nodes$node))
      )
    )
    igraph::V(net)$size <- igraph::degree(net)
    igraph::E(net)$score <- igraph::E(net)$weight
    # Plot with ggraph
    # ggraph is a package based on ggplot2,
    # the syntax is similar to regular ggplot2
    col <- RColorBrewer::brewer.pal(8, name = node.color)
    colorN <- grDevices::colorRampPalette(colors = col)(nrow(nodes))
    names(colorN) <- nodes$node

    ggraph(net, layout = graph.layout) +
      geom_edge_link0(aes(edge_linewidth = weight),
        edge_colour = edge.color
      ) +
      geom_edge_fan(aes(edge_width = score, colour = score),
        color = edge.color, show.legend = TRUE
      ) +
      ggraph::geom_node_point(aes(
        color = degree, size = degree,
        shape = class
      ), alpha = 1.0) +
      scale_color_gradientn(colours = rev(colorN)) +
      geom_node_text(
        aes(
          filter = degree >= label.degree,
          label = name
        ),
        size = label.size,
        repel = FALSE
      ) +
      scale_edge_width(range = edge.width) +
      scale_size_continuous(name = "degree", range = node.size) +
      theme_graph(base_family = "sans")
  } else if (rem.dis.inter == TRUE && label.repel == FALSE) {
    links <- links %>%
      mutate(from_c = count(., from)$n[match(from, count(., from)$from)]) %>%
      mutate(to_c = count(., to)$n[match(to, count(., to)$to)]) %>%
      filter(!(from_c == 1 & to_c == 1)) %>%
      dplyr::select(1, 2, 3)
    nodes <- links %>%
      {
        data.frame(node = c(.$from, .$to))
      } %>%
      distinct()
    col <- RColorBrewer::brewer.pal(8, name = node.color)
    colorN <- grDevices::colorRampPalette(colors = col)(nrow(nodes))
    names(colorN) <- nodes$node

    # Create a network graph from links and nodes
    net <- igraph::graph_from_data_frame(
      d = links,
      vertices = nodes,
      directed = FALSE
    )
    # # V and E are functions of the igraph package,
    # which are used to modify the nodes (nodes) and
    # connections (links) of the network graph respectively.
    igraph::V(net)$degree <- igraph::degree(net) # 每个节点连接的节点数
    igraph::V(net)$class <- c(
      rep(
        colnames(network.data)[1],
        length(intersect(network.data[, 1], nodes$node))
      ),
      rep(
        colnames(network.data)[2],
        length(intersect(network.data[, 2], nodes$node))
      ),
      rep(
        colnames(network.data)[3],
        length(intersect(network.data[, 3], nodes$node))
      )
    )
    igraph::V(net)$size <- igraph::degree(net)
    igraph::E(net)$score <- igraph::E(net)$weight

    # Plot with ggraph
    p <- ggraph(net, layout = graph.layout) +
      geom_edge_link0(aes(edge_linewidth = weight),
        edge_colour = edge.color
      ) +
      geom_edge_fan(aes(edge_width = score, colour = score),
        color = edge.color, show.legend = TRUE
      ) +
      ggraph::geom_node_point(aes(
        color = degree,
        size = degree,
        shape = class
      ), alpha = 1.0) +
      scale_color_gradientn(colours = rev(colorN)) +
      ggraph::geom_node_text(aes(filter = degree >= label.degree, label = name),
        size = label.size, repel = FALSE
      ) +
      scale_edge_width(range = edge.width) +
      scale_size_continuous(name = "degree", range = node.size) +
      theme_graph(base_family = "sans")
    return(p)
  }
}
