#' Integrated scRNA-seq reference for breast tissue (subset)
#'
#' Integrated scRNA-seq reference using random subset of two scRNA-seq datasets for breast tissue. One sample is from Komen tissue bank.
#' The other sample is from Gtex snRNA-seq data.
#'
#' @format
#' a \code{\link[Seurat]{Seurat-class}} object with several metadata included:
#' \describe{
#'  \item{cellid}{cell barcodes}
#'  \item{celltypes}{annotated cell types}
#'  \item{subjectid}{subject id}
#'  \item{cohort}{data cohort, either "komentissuebank" or "gtex"}
#' }
#'
#' @source
#' \url{https://gtexportal.org/home/singleCellOverviewPage}
#' \url{https://pubmed.ncbi.nlm.nih.gov/33763657/}
"refdata"



#' Reference scRNA-seq data for mouse cortical cell taxonomy from the Allen Institute (subset)
#'
#' The data was generated using SMART-Seq2 protocol. 20 cells were randomly selected for each cell type. 15,000 highly variable
#' genes were selected.
#'
#' @format
#' a \code{\link[Seurat]{Seurat-class}} object with several metadata included:
#' \describe{
#'  \item{cellid}{cell barcodes}
#'  \item{celltypes}{annotated cell types}
#'  \item{subjectid}{subject id}
#' }
#'
#' @source
#' \url{https://www.nature.com/articles/nn.4216}
"refdata_brain"

#' relationship between average tpm and ffpe artifacts based on a gam model
#'
#' The model was learnt using seven paired FFPE-FFzn samples from BBD patients
#'
#' @format
#' a \code{\link[gamObject]{gamObject}} object
#'
"ffpemodel"

