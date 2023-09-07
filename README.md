---
title: "TCMR"
description: a visualization tool for traditional Chinese medicine network pharmacology data analysis
author: Jinkun Liu
time: 2023-7-10
---

# 0. loaded R package 
```{r}
library(TCMR, quietly= TRUE)
library(clusterProfiler, quietly= TRUE)
library(org.Hs.eg.db, quietly= TRUE)
library(DOSE, quietly= TRUE)
library(tidyr)
library(stringr)
library(grDevices)
library(ggsankey)
library(ggvenn)
library(venn)
library(UpSetR)
```

# 1. tcm_comp
```{r}
xfbdf.compostion <- data.frame(
  herb = c(
    "mahuang", "kuxingren", "shengshigao",
    "shengyiren", "maocangzhu", "guanghuoxiang",
    "qinghaocao", "mabiancao", "ganlugen", "tinglizi",
    "huajuhong", "shenggancao", "huzhang"
  ),
  weight = c(6, 15, 30, 30, 10, 15, 12, 30, 30, 15, 15, 10, 20)
)
tcm_comp(xfbdf.compostion)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/tcm_comp.png" style="zoom: 25%;" />

# 2. herb_target

```{r}
herbs<-c('麻黄', '甘草','苦杏仁','石膏',
          '薏苡仁', '苍术', '青蒿', '猪苓',
          '马鞭草', '葶苈子','化橘红', 
          '虎杖', '广藿香','芦根')
xfbdf <- herb_target(herbs, type = "Herb_cn_name")
head(xfbdf)
```

```R
herb      molecule_id molecule target
cang zhu   MOL000173  wogonin   NOS2
cang zhu   MOL000173  wogonin  PTGS1
cang zhu   MOL000173  wogonin   ESR1
cang zhu   MOL000173  wogonin     AR
cang zhu   MOL000173  wogonin  SCN5A
cang zhu   MOL000173  wogonin  PPARG
```

```r
herbs2 <- c("ma huang", "ku xing ren")
fufang2 <- herb_target(herbs2, type ="Herb_name_pin_yin")
head(fufang2)
```

```r
herb         molecule_id molecule target
ku xing ren   MOL010921  estrone    D1R
ku xing ren   MOL010921  estrone   DRD1
ku xing ren   MOL010921  estrone  CHRM3
ku xing ren   MOL010921  estrone     F2
ku xing ren   MOL010921  estrone  CHRM1
ku xing ren   MOL010921  estrone  CHRM5
```



# 3. target_herb

```{r}
gene <- c("MAPK1", "JUN", "FOS", "RAC1", "IL1", "IL6")
herb.data <- target_herb(gene)
head(herb.data)
```

```R
   herb    molecule_id      molecule   target
ai di cha   MOL000422      kaempferol    JUN
ai di cha   MOL000098       quercetin    FOS
ai di cha   MOL000098       quercetin  MAPK1
ai di cha   MOL000098       quercetin    JUN
ai di cha   MOL000098       quercetin    IL6
    ai ye   MOL000358 beta-sitosterol    JUN
```



# 4. tcm_net

```{r}
data("xfbdf", package = "TCMR")
network.data <- xfbdf %>%
  dplyr::select(herb, molecule, target) %>%
  sample_n(100, replace = FALSE) %>%
  as.data.frame()
tcm_net(network.data, 
        label.degree = 3,
        rem.dis.inter = TRUE)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/tcm_net.png" style="zoom:25%;" />

# 5. degree_plot

```{r}
degree_plot(xfbdf,plot.set='horizontal')
```
<img src="/Users/ljk/Desktop/TCMR/文章图表/degree_plot2.png" style="zoom:25%;" />

# 6. tcm_sankey

```{r}
sankey.data <- xfbdf[sample(nrow(xfbdf), 30), ]
 tcm_sankey(sankey.data,
   text.size = 3,
   text.position = 1
 )
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/tcm_sankey.png" alt="tcm_sankey" style="zoom:25%;" />

# 7. tcm_alluvial

```{r}
alluvial.data <- xfbdf[sample(nrow(xfbdf), 30), ]
 tcm_alluvial(alluvial.data,
   text.size = 3,
   text.position = 1
 )
 
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/tcm_alluvial.png" style="zoom:25%;" />

# 8. bar_plot

```{r}
data(xfbdf, package = "TCMR")
eg <- bitr(unique(xfbdf$target), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
KK <- enrichKEGG(
  gene = eg$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)
KK <- setReadable(KK, "org.Hs.eg.db", keyType = "ENTREZID")
bar_plot(KK,title = "KEGG")
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/barplot_kk.png" style="zoom:25%;" />



```r
BP <- enrichGO(
  gene = eg$ENTREZID,
  "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE
)
bar_plot(BP,title = "biological process")
```



<img src="/Users/ljk/Desktop/TCMR/文章图表/barplot_bp.png" style="zoom:25%;" />

# 9. dot_plot

```{r}
dot_plot(KK, title = "KEGG")
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/dot_plot_KK.png" style="zoom:25%;" />

# 10. lollipop_plot

```{r}
lollipop_plot(KK, title = "KEGG")
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/lollipop_plot_KK.png" style="zoom:25%;" />

# 11. cir_plot

```{r}
cir_plot(KK)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/cir_plot_kk.png" style="zoom:25%;" />



# 12. pathway_cirplot

```{r}
#Filter the pathways to display
newdata<-KK %>% clusterProfiler.dplyr::slice(11, 15, 17, 33, 53)
pathway_cirplot(newdata)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/pathway_cirplot.png" style="zoom:25%;" />





# 13. pathway_ccplot

```{r}
pathway_ccplot(newdata)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/pathway_ccplot.png" style="zoom:25%;" />

# 14. bubble_plot

```{r}
bubble_plot(KK)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/bubble_plot.png" style="zoom:25%;" />

# 15. dot_sankey and dot_sankey2

```{r}
dot_sankey(newdata, dot.x = 0.35, dot.y = 0.25)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/dot_sankey.png" style="zoom:25%;" />



```r
data(kegg.filter, package = "TCMR")
#Because not all pathways and genes are what we want to display, and displaying too many genes at the same time will cause text overlap.
#Therefore, we built the dot.sankey function. 
#"kegg.filter" was derived from the data frame after kegg-enriched data is screened for pathways and genes.
dot_sankey2(kegg.filter)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/原始文件/dot_sankey2.png" style="zoom:25%;" />



# 16. go_barplot

```{r}
eg <- bitr(unique(xfbdf$target), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
go.diff <- enrichGO(
  gene = eg$ENTREZID,
  org.Hs.eg.db,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  ont = "all",
  readable = T
)
go_barplot(go.diff)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/go_barplot(go.diff).png" style="zoom:25%;" />



# 17. go_dotplot

```{r}
go_dotplot(go.diff)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/go_dotplot(go.diff).png" style="zoom:25%;" />

# 18. go_lollipop

```{r}
go_lollipop(go.diff)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/go_lollipop(go.diff).png" style="zoom:25%;" />

# 19. go_cir

```{r}
go_cir(go.diff)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/go_cir.png" style="zoom:80%;" />



# 20.  tcm_sankey_dot

```R
 KK2 <- KK %>% mutate(richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
 path <- separate_rows(KK2@result, geneID, sep = "/")
 data_sankey <- left_join(xfbdf, path,
   by = c("target" = "geneID"),
   relationship = "many-to-many"
 ) %>%
   distinct() %>%
   drop_na() %>%
   sample_n(30, replace = FALSE) %>%
   as.data.frame()
 tcm_sankey_dot(data_sankey)

```

<img src="/Users/ljk/Desktop/TCMR/文章图表/tcm_sankey_dot.png" style="zoom:25%;" />

```R
data_sankey2<- data_sankey %>% 
               dplyr::select(herb, molecule, target, Description)
tcm_sankey(data_sankey2,text.position = 1)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/tcm_sankey2.png" style="zoom:25%;" />



# 21. tcm_alluvial_dot

```R
tcm_alluvial_dot(data_sankey)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/tcm_alluvial_dot.png" style="zoom:25%;" />

```R
tcm_alluvial(data_sankey2, text.position = 1)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/tcm_alluvial2.png" style="zoom:25%;" />



# 22. ppi_plot

```{r}
data(string, package = "TCMR")
data <- string %>%
  dplyr::rename(
    from = X.node1,
    to = node2,
    weight = combined_score
  ) %>%
  dplyr::select(from, to, weight) %>%
  dplyr::distinct()

ppi_plot(string,
    label.degree = 1,
    nodes.color = "Spectral",
    label.repel = TRUE)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/ppi_plot2.png" style="zoom:80%;" />

# 23. dock_plot

```{r}
data <- matrix(rnorm(81), 9, 9)
data[1:9, seq(1, 9, 2)] <- data[1:9, seq(1, 9, 2)] - 4
colnames(data) <- paste("molecule", 1:9, sep = "")
rownames(data) <- paste("target", 1:9, sep = "")
data <- round(data, 2)
dock_plot(data)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/dock_plot.png" style="zoom:25%;" />

```R
dock_plot(data, shape = "circle", legend.height = 3)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/dock_plot2.png" style="zoom:25%;" />

# 24. venn_plot

```{r}
data(venn.data, package = "TCMR")
venn_plot(venn.data, type = 1)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/venn_plot(venn, type=1).png" style="zoom:25%;" />



```R
venn_plot(venn.data, type = 3)
```

<img src="/Users/ljk/Desktop/TCMR/文章图表/venn_plot(venn, type=3).png" style="zoom:30%;" />

# 25. etcm

```{r}
data("mahuang", package = "TCMR")
data <- etcm(mahuang, herb = "ma huang")
head(data)
```

```R
     herb   molecule target    QED
 ma huang Kaempferol   ACTB 0.9643
 ma huang Kaempferol    AHR 0.9643
 ma huang Kaempferol AKR1C1 0.8621
 ma huang Kaempferol   AKT1  0.963
 ma huang Kaempferol ATP5A1 0.9643
 ma huang Kaempferol  ATP5B 0.9643
```

