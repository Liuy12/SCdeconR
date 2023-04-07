
#' Integration of single-cell/nuclei RNA-seq data as reference
#' 
#' Integration of scRNA-seq or snRNA-seq data using either harmony or seurat.
#' 
#' @param ref_list a character vector of data paths to scRNA-seq/snRNA-seq. See \code{data_type} for accepted data types.
#' @param phenodata_list a character vector of data paths to metadata for elements in ref_list. All metadata within phenodata_list should have consistent column names. Columns represent
#' cell attributes, such as cell type, rows represent cells. Each element in \code{phenodata_list} should at least contain the first two columns as:
#' \enumerate{
#'  \item cell barcodes
#'  \item cell types
#' }
#' @param data_type data type of the input scRNA-seq/snRNA-seq data. Could be either a single character value from "cellranger", "h5", "matrix", or 
#' a vector/list of values with the same length as ref_list indicating the data type for each element.
#' @param method character value specifying the method to use. Has to be one of "harmony" or "seurat". See details for more information.
#' @param group_var a vector of character values indicating which variables within phenodata_list metadata to use for integration. Only applicable when method
#' is set to "harmony".
#' @param nfeature_rna minimum # of features with non-zero UMIs. Cells with # of features lower than nfeature_rna will be removed. Default to 200.
#' @param percent_mt maximum percentage of mitochondria (MT) mapped UMIs. Cells with MT percentage higher than percent_mt will be removed. Default to 40.
#' @param vars_to_regress a list of character values indicating the variables to regress for SCTransform normalization step. Default is to regress
#' out MT percentage ("percent_mt") & cell cycle effects ("phase")
#' @param ex_features a vector of character values indicating genes to exclude from anchor features. Those genes will not be considered as 
#' anchor features for integration, but will still be present in the integrated data.
#' @param cluster logical value indicating whether to perform clustering on the integrated data. If TRUE, unsupervised clustering will be performed, and 
#' the results will be saved in "seurat_clusters" metadata in the output Seurat object. 
#' @param resolution numeric value specifying resolution to use when cluster is set to TRUE. 
#' @param verbose logical value indicating whether to print messages. 
#' @param ... additional parameters passed to \code{\link[Seurat]{SCTransfrom}}.
#' 
#' @return a \code{\link[Seurat]{Seurat-class}} object.
#' 
#' @details 
#' data_type can be chosen from:
#' \describe{
#'  \item{cellranger}{path to a directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv files outputted by 10x's cell-ranger}
#'  \item{h5}{path to .h5 file outputted by 10x's cell-ranger}
#'  \item{matrix}{path to a matrix-like file, with rows representing genes, columns representing cells.}
#' }
#' 
#' SCTransform with \code{vst.flavor = "v2"} is used for normalization of individual data. Integration methods can be chosen from either "harmony" 
#' or "seurat". Harmony typically is more memory efficient and, recommended if you have large # of cells for integration. 
#' 
#' @export 
#' 
#' @examples 
#' \dontrun{
#' ## random subset of two scRNA-seq datasets for breast tissue
#' ref_list <- c(paste0(system.file("extdata", package = "SCdeconR"), "/refdata/sample1"),
#'               paste0(system.file("extdata", package = "SCdeconR"), "/refdata/sample2"))
#' phenodata_list <- c(paste0(system.file("extdata", package = "SCdeconR"), "/refdata/phenodata_sample1.txt"),
#'                     paste0(system.file("extdata", package = "SCdeconR"), "/refdata/phenodata_sample2.txt"))
#' 
#' ## Register backend for parallel processing
#' registerDoFuture()
#' plan("multisession", workers = 4)
#' 
#' ## construct integrated reference data
#' refdata <- construct_ref(ref_list = ref_list, 
#'                          phenodata_list = phenodata_list, 
#'                          data_type = "cellranger", 
#'                          method = "harmony", 
#'                          group_var = "subjectid", 
#'                          nfeature_rna = 50, 
#'                          vars_to_regress = "percent_mt")
#' }


construct_ref <- function(
    ref_list,
    phenodata_list,
    data_type = c("cellranger", "h5", "matrix"),
    method = c("harmony", "seurat"),
    group_var,
    nfeature_rna = 200,
    percent_mt = 40,
    vars_to_regress = c("percent_mt", "phase"),
    ex_features = NULL,
    cluster = TRUE,
    resolution = 0.8,
    verbose = TRUE, ...) {
  phenodata_list <- lapply(phenodata_list, function(i) data.table::fread(file = i, check.names = FALSE, header = TRUE, data.table = FALSE))
  phenodata_list_names <- lapply(phenodata_list, function(i) colnames(i))
  if (any(Reduce(intersect, phenodata_list_names) != phenodata_list_names[[which.max(sapply(phenodata_list_names, function(i) length(i)))]])) {
    stop("phenodata_list contains inconsistent colnames")
  }
  if (!group_var %in% colnames(phenodata_list[[1]])) {
    stop(paste0("groupvar not found in column names of phenodata_list: "), paste0(phenodata_list_names[[1]], collapse = ","))
  }
  if (!all(sapply(data_type, function(i) i %in% c("cellranger", "h5", "matrix")))) stop(paste0("data_type has to be one of: ", paste0(c("cellranger", "h5", "matrix"), collapse = ",")))
  if (length(data_type) > 1 && length(data_type) != length(ref_list)) stop("data_type could be either a single character value from 'cellranger', 'h5', 'matrix', or a vector of values with the same length as ref_list indicating the data type for each element")
  if (length(data_type) ==1) data_type <- rep(data_type, times = length(ref_list))
  object_list <- foreach(i = seq_along(ref_list), .packages = "Seurat") %dopar%
    load_scdata(
      ref = ref_list[[i]],
      data_type = data_type[[i]],
      meta_info = phenodata_list[[i]],
      nfeature_rna = nfeature_rna,
      percent_mt = percent_mt,
      vars_to_regress = vars_to_regress,
      id = i,
      verbose = verbose, ...
    )
  if (verbose) message("select integration features")
  object_features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 3000, verbose = verbose)
  object_features <- setdiff(object_features, ex_features)
  if (method == "seurat") {
    object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = object_features, verbose = verbose, assay = "SCT")
    merge_anchor <- FindIntegrationAnchors(object.list = object_list, dims = 1:30, normalization.method = "SCT", anchor.features = object_features, verbose = verbose)
    rm(list = c("object_list"))
    gc(verbose = verbose)
    message("perform integration")
    merge_data <- IntegrateData(anchorset = merge_anchor, normalization.method = "SCT", verbose = verbose, dims = 1:30)
    rm(list = c("merge_anchor"))
    gc(verbose = verbose)
  } else if (method == "harmony") {
    if (verbose) message("Generating merged data for harmony")
    merge_data <- merge(object_list[[1]], y = object_list[2:length(object_list)], merge.data = TRUE)
    rm(list = c("object_list"))
    gc(verbose = verbose)
    VariableFeatures(merge_data) <- object_features
    merge_data <- RunPCA(object = merge_data, assay = "SCT", verbose = verbose, npcs = 30)
    if (verbose) message("Running harmony")
    merge_data <- harmony::RunHarmony(
      object = merge_data,
      assay.use = "SCT",
      reduction = "pca",
      dims.use = 1:30,
      group.by.vars = group_var,
      plot_convergence = FALSE,
      verbose = verbose
    )
  }
  if (cluster) {
    reduction <- ifelse(method == "harmony", "harmony", "pca")
    merge_data <- RunPCA(merge_data, npcs = 30, verbose = verbose, features = VariableFeatures(merge_data))
    merge_data <- RunTSNE(merge_data, reduction = reduction, dims = 1:30)
    merge_data <- RunUMAP(merge_data, reduction = reduction, dims = 1:30)
    merge_data <- FindNeighbors(merge_data, reduction = reduction, dims = 1:30)
    merge_data <- FindClusters(merge_data, resolution = resolution, verbose = verbose)
  }
  if (verbose) message("normalizing median UMI counts across SCT models")
  merge_data <- PrepSCTFindMarkers(merge_data, verbose = verbose)
  return(merge_data)
}

#' Load, filter and normalize scRNA-seq/snRNA-seq data
#' 
#' Load and preprocess scRNA-seq/snRNA-seq data using seurat SCTransform workflow.
#' 
#' @param ref path to scRNA-seq/snRNA-seq data.
#' @param data_type a character value specifying data type of the input scRNA-seq/snRNA-seq data, should be one of "cellranger", "h5", "matrix".
#' @param meta_info a data.frame with rows representing cells, columns representing cell attributes. 
#' @param nfeature_rna minimum # of features with non-zero UMIs. Cells with # of features lower than nfeature_rna will be removed. Default to 200.
#' @param percent_mt maximum percentage of mitochondria (MT) mapped UMIs. Cells with MT percentage higher than percent_mt will be removed. Default to 40.
#' @param vars_to_regress a list of character values indicating the variables to regress for SCTransform normalization step. Default is to regress
#' out MT percentage ("percent_mt") & cell cycle effects ("phase")
#' @param id a character value specifying project or sample id. Only used for printing purposes.
#' @param verbose logical value indicating whether to print messages. 
#' @param ... additional parameters passed to \code{\link[Seurat]{SCTransfrom}}. 
#' 
#' @return a \code{\link[Seurat]{Seurat-class}} object.
#' 
#' @details 
#' For more details, refer to \code{\link{construct_ref}}
#' 
#' @export 
#' 
#' @examples 
#' \dontrun{
#' ref_list <- c(paste0(system.file("extdata", package = "SCdeconR"), "/refdata/sample1"),
#'               paste0(system.file("extdata", package = "SCdeconR"), "/refdata/sample2"))
#' phenodata_list <- c(paste0(system.file("extdata", package = "SCdeconR"), "/refdata/phenodata_sample1.txt"),
#'                     paste0(system.file("extdata", package = "SCdeconR"), "/refdata/phenodata_sample2.txt"))
#' tmp <- load_scdata(
#'   ref = ref_list[[1]], 
#'   data_type = c("cellranger"), 
#'   meta_info = fread(file = phenodata_list[[1]], check.names = FALSE, header = TRUE), 
#'   nfeature_rna = 50, 
#'   vars_to_regress = c("percent_mt"), 
#'   id = 1, 
#'   verbose = TRUE)
#' }


load_scdata <- function(
    ref,
    data_type = c("cellranger", "h5", "matrix"),
    meta_info,
    nfeature_rna = 200,
    percent_mt = 40,
    vars_to_regress = c("percent_mt", "phase"),
    id,
    verbose, ...) {
  if (data_type == "cellranger") {
    data <- Read10X(data.dir = ref)
  } else if (data_type == "h5") {
    data <- Read10X_h5(filename = ref)
  } else {
    data <- data.table::fread(file = ref, data.table = FALSE)
  }
  if (verbose) message(paste0("Loading sample completed for sample ", id))
  rownames(meta_info) <- meta_info[[1]]
  if (length(intersect(rownames(meta_info), colnames(data))) != length(union(rownames(meta_info), colnames(data)))) {
    stop(paste0("cell barcodes mismatch betweeen expression data and meta data for sample ", id))
  }
  meta_info <- meta_info[match(colnames(data), rownames(meta_info)), ]
  if (verbose) message(paste0("create seurat object for sample ", id))
  data <- CreateSeuratObject(counts = data, meta.data = as.data.frame(meta_info))
  if (verbose) message(paste0("calculate mt percent for sample "), id)
  data[["percent_mt"]] <- PercentageFeatureSet(object = data, pattern = "^MT-|^mt-")
  cellids <- colnames(data)[data$nFeature_RNA >= nfeature_rna & data$percent_mt <= percent_mt]
  data <- subset(data, cells = cellids)
  if ("phase" %in% vars_to_regress) {
    if (verbose) message(paste0("predict cell cycle phase for sample ", id))
    data <- NormalizeData(data, assay = "RNA", verbose = verbose)
    data("cc.genes")
    data <- CellCycleScoring(data, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE, verbose = verbose)
  }
  if (verbose) message(paste0("perform sctransform for sample ", id))
  data <- SCTransform(data, variable.features.n = 3000, vst.flavor = "v2", vars.to.regress = vars_to_regress, verbose = verbose, return.only.var.genes = FALSE, ...)
  return(data)
}
