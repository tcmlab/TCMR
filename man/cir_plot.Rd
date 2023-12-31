% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cir_plot.R
\name{cir_plot}
\alias{cir_plot}
\title{Chord diagram display of KEGG or GO results}
\usage{
cir_plot(
  data,
  top = 15,
  id.color = "RdBu",
  id.cex = 0.8,
  gene.color = "#CFB0D4",
  gene.cex = 0.8,
  ratio.color = c("#6291bd", "#B1D1E7"),
  ratio.cex = 0.8,
  p.color = c("#eda9aa", "#c25254"),
  p.cex = 0.8,
  richFactor.color = "RdBu",
  richFactor.cex = 0.7,
  legend.position = c(0.8, 0.45),
  text.width = 35,
  ...
)
}
\arguments{
\item{data}{data.frame}

\item{top}{number of pathways}

\item{id.color}{The color of the first layer ID of the path}

\item{id.cex}{The text size of the first layer ID of the path}

\item{gene.color}{The color of the number of all genes in this pathway}

\item{gene.cex}{The font size of the number of all genes in this pathway}

\item{ratio.color}{（1）The color of the number of genes enriched this time on this pathway
（2）Other genes involved in this enrichment analysis}

\item{ratio.cex}{The text size of the numbers involved in this enrichment analysis}

\item{p.color}{padjust color}

\item{p.cex}{The text size of padjust}

\item{richFactor.color}{richFactor color}

\item{richFactor.cex}{The text size of richFactor}

\item{legend.position}{legend position}

\item{text.width}{legend text width}

\item{...}{additional parameters}
}
\value{
chord diagram
}
\description{
Chord diagram display of KEGG or GO results
}
\examples{
\dontrun{
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
eg <- bitr(unique(xfbdf$target),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)
KK <- enrichKEGG(
  gene = eg$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)
KK <- setReadable(KK, "org.Hs.eg.db", keyType = "ENTREZID")
cir_plot(KK)
}
}
