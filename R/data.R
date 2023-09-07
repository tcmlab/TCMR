#' Xuanfei Baidu granule (XFBD)
#'
#' @description
#' ephedra 6 g, bitter almond 15 g, raw gypsum 30 g,
#' raw almond 30 g, raw plaster 30 g,
#' Atractylodes macrocephala 10 g, patchouli 15 g,
#' artemisia grass 12 g, Polygonum cuspidatum 20 g,
#' verbena 30 g, dried Reed root 30 g,
#' Fructus thunbergii 15 g, tangerine red 15 g,
#' raw licorice 10 g.
#'
#' A data frame with 5024 rows and 4 variables
#' @source \url{http://www.nhc.gov.cn/cms-search/downFiles/a449a3e2e2c94d9a856d5faea2ff0f94.pdf}
#' @examples
#' data("xfbdf", package = "TCMR")
#' head(xfbdf)
"xfbdf"


#' string
#'
#' Protein-protein interaction analysis was performed using STRING database.
#' A data frame with 2234 rows and 3 variables
#' @source \url{https://string-db.org/}
#' @examples
#' string.data <- data(string)
#' head(string.data)
#'
"string"

#' ma huang
#'
#' The datasets was from ETCM database.
#' A data frame with 28 rows and 2 variables
#' @source \url{http://www.tcmip.cn/ETCM/index.php/Home/Index/yc_details.html?id=212}
#' @examples
#' head(mahuang)
#'
"mahuang"

#' KK
#'
#' The results of R clusterprofiler package for KEGG results
#' S4 object of class enrich results
#' @source R clusterprofiler package for KEGG results
#' @examples
#' head(KK@result)
#'
"KK"

#' kegg.filter
#'
#' The results after KEGG/GO screening of pathways and genes
#' A data frame with 8 rows and 9 variables
#' @source example data from the KEGG results
#' @examples
#' head(kegg.filter)
#'
"kegg.filter"

#' Venn diagram example data
#'
#' The datasets was from xfbdf and https://doi.org/10.1038/s41556-021-00796-6.
#' A data frame with 6439 rows and 2 variables
#' @source \url{https://www.nature.com/articles/s41556-021-00796-6#Sec47}
#' @examples
#' head(venn.data)
"venn.data"

#' TRRUST database
#' @description TRRUST contains 8,444 and 6,552 TF-target regulatory
#' relationships of 800 human TFs and 828 mouse TFs, respectively.
#'
#' The datasets was from TRRUST database
#' A data frame with 16543 rows and 6 variables
#' @source \url{https://www.grnpedia.org/trrust/}
#' @examples
#' head(trrust)
"trrust"

#' Herbal Properties and Meridians
#' @description Chinese Pharmacopoeia commonly used herbs property and flavor, tropism
#' The datasets was from Chinese Pharmacopoeia
#' A data frame with 409 rows and 5 variables
#' @source \url{https://db.ouryao.com/yd2020/index.php?cid=1}
#' @examples
#' head(herb_pm)
"herb_pm"
