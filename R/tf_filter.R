#' Finding transcription factors and their target genes
#'
#' @param data The data must be character as human or mouse gene symbol
#'
#' @return data.frame
#' @export
#' @examples
#' data(xfbdf, package = "TCMR")
#' newdata <- tf_filter(xfbdf$target)
#' head(newdata)
tf_filter <- function(data) {
  if (is.character(data)) {
    tf.data <- trrust[trrust$TF %in% data, ] %>% as.data.frame()
    return(tf.data)
  } else {
    print("The data must be character as human or mouse gene symbol.")
  }
}
