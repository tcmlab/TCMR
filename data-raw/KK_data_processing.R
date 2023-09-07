library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
data(xfbdf, package = "TCMR")
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
