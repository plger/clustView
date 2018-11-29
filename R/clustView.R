#' clustView
#'
#' @param seurat An object of class `seurat`
#' @param enableClustree Logical; whether or not to enable the
#' clustering tree (default TRUE). Turn off to reduce loading
#' time.
#'
#' @export
clustView <- function(seurat, enableClustree=T){
  library(shiny)
  library(Seurat)
  shinyApp( ui=clustView.ui(),
            server=clustView.server(seurat, enableClustree=enableClustree), 
            options=list(launch.browser=T)
  )
}