#' Integrated scRNA-seq reference for breast tissue (subset)
#' 
#' Integrated scRNA-seq reference using random subset of two scRNA-seq datasets for breast tissue. One sample is from Komen tissue bank. 
#' The other sample is from Gtex snRNA-seq data.
#' 
#' @format `refdata`
#' a \code{\link[Seurat]{Seurat-class}} object with several metadata included:
#' \describe{
#'  \item{cellid}{cell barcodes}
#'  \item{celltypes}{annotated cell types}
#'  \item{subjectid}{subject id}
#'  \item{cohort}{data cohort, either "komentissuebank" or "gtex"}
#' }
#' 
#' @source 
#' <https://gtexportal.org/home/singleCellOverviewPage>
#' <A single-cell atlas of the healthy breast tissues reveals clinically relevant clusters of breast epithelial cells>
"refdata"