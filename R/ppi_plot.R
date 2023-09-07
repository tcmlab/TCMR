#' Protein-protein interaction analysis
#'
#' @param data Human protein-protein interaction (PPI)
#' data were downloaded from the STRING database.
#' @param nodes.color  node color: see "RColorBrewer::display.brewer.all()"
#' @param node.size  node size
#' @param label.size markup text size
#' @param label.degree The node degree is the number of connections that
#' the node has with the other nodes.
#' Nodes with connections greater than or
#' equal to degree will be displayed.
#' @param edge.color edge color
#' @param edge.width edge width
#' @param rem.dis.inter remove single free unconnected nodes
#' @param graph.layout Network Diagram Layout:
#' "kk", nicely", "grid","sphere", "randomly",
#' "bipartite", "star","tree", nicely", "grid",
#' "sphere", "randomly","gem", "graphopt","lgl",
#' "mds", "sugiyama"
#' @param label.repel label repel

#' @return figure: network diagrams
#' @export
#' @import dplyr
#' @import ggplot2
#' @import ggraph
#' @importFrom RColorBrewer brewer.pal
#' @importFrom igraph V
#' @importFrom igraph E
#' @importFrom igraph graph_from_data_frame
#' @examples
#' \dontrun{
#' data(string, package = "TCMR")
#' data <- string %>%
#'   dplyr::rename(
#'     from = X.node1,
#'     to = node2,
#'     weight = combined_score
#'   ) %>%
#'   dplyr::select(from, to, weight) %>%
#'   dplyr::distinct()
#'
#' ppi_plot(data, label.repel = FALSE)
#' ppi_plot(data,
#'   label.degree = 1,
#'   nodes.color = "Spectral",
#'   label.repel = TRUE
#' )
#' }
ppi_plot <- function(data,
                     nodes.color = "RdBu",
                     node.size = c(1, 10),
                     label.size = 4,
                     label.degree = 5,
                     label.repel = TRUE,
                     edge.color = "lightgrey",
                     edge.width = c(0.2, 2),
                     rem.dis.inter = FALSE,
                     graph.layout = "kk") {
  links <- data
  # node data
  nodes <- links %>%
    {
      data.frame(gene = c(.$from, .$to))
    } %>%
    dplyr::distinct()

  # Create a network graph from links and nodes
  if (rem.dis.inter == FALSE && label.repel == TRUE) {
    net <- igraph::graph_from_data_frame(
      d = links,
      vertices = nodes,
      directed = FALSE
    )
    # V and E are functions of the igraph package,
    # which are used to modify the nodes (nodes) and
    # connections (links) of the network graph respectively.
    igraph::V(net)$degree <- igraph::degree(net)
    igraph::V(net)$size <- igraph::degree(net)
    igraph::E(net)$score <- igraph::E(net)$weight
    # Plot with ggraph
    ggraph(net, layout = graph.layout) +
      geom_edge_fan(
        aes(
          edge_width = score,
          colour = score
        ),
        color = edge.color,
        show.legend = TRUE
      ) +
      geom_node_point(aes(color = degree, size = size), alpha = 1.0) +
      scale_color_gradientn(
        colours =
          rev(RColorBrewer::brewer.pal(8, nodes.color))
      ) +
      geom_node_text(aes(filter = degree >= label.degree, label = name),
        size = label.size,
        repel = TRUE
      ) +
      ggraph::scale_edge_width(range = edge.width) +
      scale_size_continuous(name = "degree", range = node.size) +
      theme_graph(base_family = "sans")
  } else if (rem.dis.inter == TRUE && label.repel == TRUE) {
    # remove dissociative interactions
    # If the from of a link in the links data frame
    # only appears once, and the to appears only once, remove it
    links_2 <- links %>%
      mutate(from_c = count(., from)$n[match(from, count(., from)$from)]) %>%
      mutate(to_c = count(., to)$n[match(to, count(., to)$to)]) %>%
      filter(!(from_c == 1 & to_c == 1)) %>%
      dplyr::select(1, 2, 3)
    # new node data
    nodes_2 <- links_2 %>%
      {
        data.frame(gene = c(.$from, .$to))
      } %>%
      distinct()
    # Create a network diagram
    net_2 <- igraph::graph_from_data_frame(
      d = links_2,
      vertices = nodes_2,
      directed = FALSE
    )
    # Add the necessary parameters
    igraph::V(net_2)$degree <- igraph::degree(net_2)
    igraph::V(net_2)$size <- igraph::degree(net_2)
    igraph::E(net_2)$score <- igraph::E(net_2)$weight
    ggraph(net_2, layout = graph.layout) +
      geom_edge_fan(aes(edge_width = score, colour = score),
        color = edge.color, show.legend = TRUE
      ) +
      geom_node_point(aes(color = degree, size = size), alpha = 1.0) +
      scale_color_gradientn(colours = rev(
        RColorBrewer::brewer.pal(8, nodes.color)
      )) +
      geom_node_text(aes(filter = degree >= label.degree, label = name),
        size = label.size, repel = TRUE
      ) +
      ggraph::scale_edge_width(range = edge.width) +
      scale_size_continuous(name = "degree", range = node.size) +
      theme_graph(base_family = "sans")
  } else if (rem.dis.inter == FALSE && label.repel == FALSE) {
    net <- igraph::graph_from_data_frame(
      d = links,
      vertices = nodes,
      directed = FALSE
    )
    igraph::V(net)$degree <- igraph::degree(net)
    igraph::V(net)$size <- igraph::degree(net)
    igraph::E(net)$score <- igraph::E(net)$weight

    # Plot with ggraph
    ggraph(net, layout = graph.layout) +
      geom_edge_fan(aes(edge_width = score, colour = score),
        color = edge.color,
        show.legend = TRUE
      ) +
      geom_node_point(aes(color = degree, size = size), alpha = 1.0) +
      scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(
        8, nodes.color
      ))) +
      geom_node_text(aes(filter = degree >= label.degree, label = name),
        size = label.size,
        repel = FALSE
      ) +
      ggraph::scale_edge_width(range = edge.width) +
      scale_size_continuous(name = "degree", range = node.size) +
      theme_graph(base_family = "sans")
  } else if (rem.dis.inter == TRUE && label.repel == FALSE) {
    # remove dissociative interactions
    links_2 <- links %>%
      mutate(from_c = count(., from)$n[match(from, count(., from)$from)]) %>%
      mutate(to_c = count(., to)$n[match(to, count(., to)$to)]) %>%
      filter(!(from_c == 1 & to_c == 1)) %>%
      dplyr::select(1, 2, 3)
    # new node data
    nodes_2 <- links_2 %>%
      {
        data.frame(gene = c(.$from, .$to))
      } %>%
      distinct()
    # Create a network diagram
    net_2 <- igraph::graph_from_data_frame(
      d = links_2,
      vertices = nodes_2,
      directed = FALSE
    )
    # Add the necessary parameters
    igraph::V(net_2)$degree <- igraph::degree(net_2)
    igraph::V(net_2)$size <- igraph::degree(net_2)
    igraph::E(net_2)$score <- igraph::E(net_2)$weight
    ggraph(net_2, layout = graph.layout) +
      geom_edge_fan(aes(edge_width = score, colour = score),
        color = edge.color,
        show.legend = TRUE
      ) +
      geom_node_point(aes(color = degree, size = size), alpha = 1.0) +
      scale_color_gradientn(
        colours =
          rev(RColorBrewer::brewer.pal(8, nodes.color))
      ) +
      geom_node_text(aes(filter = degree >= label.degree, label = name),
        size = label.size,
        repel = FALSE
      ) +
      ggraph::scale_edge_width(range = edge.width) +
      scale_size_continuous(name = "degree", range = node.size) +
      theme_graph(base_family = "sans")
  }
}
