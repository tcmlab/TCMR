#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom grDevices recordPlot
#' @importFrom rlang .data
#' @importFrom rlang .env
#' @importFrom stats aggregate
#' @importFrom utils head
## usethis namespace: end
#'
#' @docType package
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "Herb_cn_name",
    "HERB_NAME_PIN_YIN",
    "Herb_name_pin_yin",
    "MOL_ID",
    "Molecule_name",
    "OB",
    "DL",
    "Target_full_name",
    "Target",
    "Disease_name",
    "X.node1",
    "node2",
    "combined_score",
    "herb",
    "molecule_id",
    "molecule",
    "target",
    "D",
    "Description",
    "GeneRatio",
    "BgRatio",
    "pvalue",
    "p.adjust",
    "qvalue",
    "geneID",
    "Count",
    "ONTOLOGY",
    "name",
    ".,",
    "from_c",
    "from",
    "to",
    "to_c",
    "score",
    "weight",
    ".",
    "ID",
    "logP",
    "next_node",
    "next_x",
    "node",
    "richFactor",
    "RichFactor",
    "rowname",
    "size",
    "tcmsp",
    "type",
    "x",
    "Candidate.Target.genes",
    "Chemical.Component",
    "Candidate Target genes",
    "Chemical Component",
    "QED",
    "source",
    "gene",
    ".groups",
    "n",
    "names_from",
    "values_from",
    "query",
    "params",
    "intersects",
    "leaf",
    "node.branch",
    "node.short_name",
    "node.size",
    "node1.node.branch",
    "y",
    "ratio"
  ))
}


