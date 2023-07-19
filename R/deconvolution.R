#' Deconvolution of bulk RNA-seq data
#'
#' Deconvolution of bulk RNA-seq data based on single-cell reference data. Eight bulk deconvolution methods, along with eight normalization
#' methods and four transformation methods are available.
#'
#' @param bulk a matrix or data.frame of unnormalizaed & untransformed bulk RNA-seq gene expression values with rows representing genes, columns
#' representing samples
#' @param ref a matrix or data.frame of untransformed scRNA-seq gene expression counts with rows representing genes, columns representing cells.
#' This data will be used to deconvolute provided bulk RNA-seq data.
#' @param phenodata a data.frame with rows representing cells, columns representing cell attributes. It should at least contain the first three
#' columns as:
#' \enumerate{
#'  \item cell barcodes
#'  \item cell types
#'  \item subject ids
#' }
#' @param filter_ref logical value indicating whether outlier genes & cells should be removed from the provided reference data. Defaults to TRUE
#' @param marker_genes a data.frame of two columns. First column represents cell types in ref; second column represents gene names of marker genes. If specified,
#' those genes will be used to construct signature matrix for mark-gene based deconvolution methods, such as CIBERSORT, OLS, nnls, FARDEEP and RLR. Default to NULL,
#' carry out differential analysis to identify marker genes for each cell type in ref.
#' @param genes_to_remove a vector of gene names to remove from the reference scRNA-seq data. Default to NULL.
#' @param min_pct_ct a numeric value indicating the minimum required proportion of expressing cells per cell type for marker gene identification. Only applicable when marker_genes
#' is NULL. Default to 0.05.
#' @param decon_method character value specifying the deconvolution method to use. Has to be one of "scaden", "CIBERSORT", "OLS", "nnls", "FARDEEP", "RLR",
#' "MuSiC", "SCDC". See details for more information.
#' @param norm_method_sc character value specifying the normalization method to use for reference data. Has to be one of "none","LogNormalize",
#' "SCTransform", "scran", "scater", "Linnorm". See details for more information.
#' @param norm_method_bulk character value specifying the normalization method to use for bulk data. Has to be one of "none", "TMM",
#' "median_ratios", "TPM". See details for more information.
#' @param trans_method_sc character value specifying the transformation method to use for both bulk & reference data. Has to be one of "none", "log2", "sqrt",
#' "vst". See details for more information.
#' @param trans_method_bulk character value specifying the transformation method to use for both bulk & reference data. Has to be one of "none", "log2", "sqrt",
#' "vst". See details for more information.
#' @param gene_length a data.frame with two columns. The first column represents gene names that match with provided bulk data. The second column
#' represents length of each gene. Only applicable when norm_method is selected as "TPM".
#' @param lfc_markers log2 fold change cutoff used to identify marker genes for deconvolution. The option only applicable to marker-gene based
#' approaches, such as CIBERSORT, OLS, nnls, FARDEEP and RLR. Only applicable when marker_genes is NULL.
#' @param marker_strategy further strategy in selecting marker genes besides applying the log2 fold change cutoff. Can be chosen from: "all", "pos_fc",
#' "top_50p_logFC" or "top_50p_AveExpr". See details for more information. Only applicable when marker_genes is NULL.
#' @param to_remove character value representing the cell type to remove from reference data. Only applicable to simulation experiments in evaluating
#' effect of cell type removal from reference.
#' @param ffpe_artifacts logical value indicating whether to add simulated ffpe artifacts in the bulk data. Only applicable to simulation experiments in
#' evaluating the effect of FFPE artifacts.
#' @param prop a matrix or data.frame of simulated cell proportion values with rows representing cell types, columns representing samples. Only applicable to simulation
#' experiments in evaluating the effect of cell type removal from reference.
#' @param cibersortpath full path to CIBERSORT.R script.
#' @param pythonpath full path to python binary where scaden was installed with.
#' @param tmpdir temporary processing directory for scaden.
#' @param remove_tmpdir a logical value indicating whether to remove tmpdir once scaden is completed. Default to TRUE.
#' @param seed random seed used for simulating FFPE artifacts. Only applicable when ffpe_artifacts is set to TRUE.
#' @param nsamples number of artificial bulk samples to simulate for scaden. Default to 1000.
#' @param verbose a logical value indicating whether to print messages. Default to FALSE.
#'
#' @return a list containing two or four elements.
#' \describe{
#'  \item{first element}{a data.frame of predicted cell-type proportions, with rows representing cell types, columns representing samples.}
#'  \item{second element}{a data.frame of fitting errors of the algorithm; first column represents sample names, second column represents RMSEs.}
#'  \item{optional third element}{a data.frame of simulated cell proportion after removing the specified cell_type. Only applicable to simulation experiments.}
#'  \item{optional fourth element}{a data.frame of marker genes used for deconvolution. Only applicable to marker-gene based deconvolution methods.}
#' }
#'
#' @details
#' decon_method should be one of the following:
#' \describe{
#'  \item{\href{https://github.com/KevinMenden/scaden}{scaden}}{a deep learning based method using three multi-layer deep neural nets. To use scaden,
#' you need to firstly install scaden via \code{pip} or \code{conda}, the provide the python path to \code{pythonpath} option.}
#'  \item{CIBERSORT}{a marker gene based support vectors regression approach. CIBERSOR does not allow redistribution. To use CIBERSORT,
#' you need to request the source code from the authors & provide the path of CIBERSORT.R script to \code{cibersortpath} option.}
#'  \item{OLS}{ordinary least squares.}
#'  \item{\href{https://cran.r-project.org/web/packages/nnls/index.html}{nnls}}{non-negative least squares.}
#'  \item{\href{https://cran.r-project.org/web/packages/FARDEEP/index.html}{FARDEEP}}{robust regression using least trimmed squares}
#'  \item{\href{https://cran.r-project.org/web/packages/MASS/index.html}{RLR}}{robust regression using an M estimator}
#'  \item{\href{https://github.com/xuranw/MuSiC}{MuSiC}}{multi-subject single-cell deconvolution}
#'  \item{\href{https://github.com/meichendong/SCDC}{SCDC}}{an ENSEMBLE method to integrate deconvolution results from different scRNA-seq datasets}
#' }
#'
#' norm_method should be one of the following:
#' \describe{
#'  \item{none}{no normalization is performed.}
#'  \item{LogNormalize}{\code{\link[Seurat]{LogNormalize}} method from seurat.}
#'  \item{TMM}{TMM method from \code{\link[edgeR]{calcNormFactors}} function from edgeR.}
#'  \item{median_ratios}{median ratio method from \code{\link[DESeq2]{estimateSizeFactors,DESeqDataSet-method}} function from DESeq2.}
#'  \item{TPM}{Transcript per million. TPM has to be chosen if ffpe_artifacts is set to TRUE.}
#'  \item{SCTransform}{\code{\link[Seurat]{SCTransfrom}} method from Seurat.}
#'  \item{scran}{\code{\link[scran]{computeSumFactors}} method from scran.}
#'  \item{scater}{\code{\link[scater]{librarySizeFactors}} method from scater.}
#'  \item{Linnorm}{\code{\link[Linnorm]{Linnorm}} method from Linnorm.}
#' }
#'
#' trans_method should be one of the following:
#' \describe{
#'  \item{none}{no transformation is performed.}
#'  \item{log2}{log2 transformation. 0.1 is added to the data to avoid logarithm of 0s.}
#'  \item{sqrt}{square root transformation.}
#'  \item{vst}{\code{\link[DESeq2]{varianceStabilizingTransformation}} method from DESeq2.}
#' }
#'
#' marker_strategy should be one of the following:
#' \describe{
#'  \item{all}{all genes passed fold change threshold will be used.}
#'  \item{pos_fc}{only genes with positive fold changes will be used.}
#'  \item{top_50p_logFC}{only genes with top 50 percent positive fold changes will be used.}
#'  \item{top_50p_AveExpr}{only genes with top 50 percent average expression will be used.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' ## generate artificial bulk samples
#' data(refdata)
#' phenodata <- data.frame(cellid = colnames(refdata),
#'                         celltypes = refdata$celltype,
#'                         subjectid = refdata$subjectid)
#' bulk_sim <- bulk_generator(ref = GetAssayData(refdata, slot = "data", assay = "SCT"),
#'                            phenodata = phenodata,
#'                            num_mixtures = 20,
#'                            num_mixtures_sprop = 1)
#'
#' ## perform deconvolution based on "OLS" algorithm
#' decon_res <- scdecon(bulk = bulk_sim[[1]],
#'                      ref = GetAssayData(refdata, slot = "data", assay = "SCT"),
#'                      phenodata = phenodata,
#'                      filter_ref = TRUE,
#'                      decon_method = "OLS",
#'                      norm_method = "none",
#'                      trans_method = "none",
#'                      marker_strategy = "all")
#' }


scdecon <- function(
    bulk,
    ref,
    phenodata,
    filter_ref = TRUE,
    marker_genes = NULL,
    genes_to_remove = NULL,
    min_pct_ct = 0.05,
    decon_method = c("scaden", "CIBERSORT", "OLS", "nnls", "FARDEEP", "RLR", "MuSiC", "SCDC"),
    norm_method_sc = c("none","LogNormalize", "SCTransform", "scran", "scater", "Linnorm"),
    norm_method_bulk = c("none", "TMM", "median_ratios", "TPM"),
    trans_method_sc = c("none", "log2", "sqrt", "vst"),
    trans_method_bulk = c("none", "log2", "sqrt", "vst"),
    gene_length = NULL,
    lfc_markers = log2(1.5),
    marker_strategy = c("all", "pos_fc", "top_50p_logFC", "top_50p_AveExpr"),
    to_remove = NULL,
    ffpe_artifacts = FALSE,
    prop = NULL,
    cibersortpath = NULL,
    pythonpath = NULL,
    tmpdir = NULL,
    remove_tmpdir = TRUE,
    seed = 1234,
    nsamples = 1000,
    verbose = FALSE) {
  if (!decon_method %in% c("CIBERSORT", "OLS", "nnls", "FARDEEP", "RLR", "MuSiC", "SCDC", "scaden")) stop(paste0("decon_method must be one of ", paste0(c("CIBERSORT", "OLS", "nnls", "FARDEEP", "RLR", "MuSiC", "SCDC", "scaden"), collapse = ",")))
  if (!norm_method_sc %in% c("none","LogNormalize", "SCTransform", "scran", "scater", "Linnorm")) stop(paste0("norm_method_sc must be one of ", paste0(c("none","LogNormalize", "SCTransform", "scran", "scater", "Linnorm"), collapse = ",")))
  if (!norm_method_bulk %in% c("none", "TMM", "median_ratios", "TPM")) stop(paste0("norm_method_bulk must be one of ", paste0(c("none", "TMM", "median_ratios", "TPM"), collapse = ",")))
  if (!trans_method_sc %in% c("none", "log2", "sqrt", "vst")) stop(paste0("trans_method_sc must be one of ", paste0(c("none", "log2", "sqrt", "vst"), collapse = ",")))
  if (!trans_method_bulk %in% c("none", "log2", "sqrt", "vst")) stop(paste0("trans_method_bulk must be one of ", paste0(c("none", "log2", "sqrt", "vst"), collapse = ",")))
  if (decon_method %in% c("CIBERSORT", "OLS", "nnls", "FARDEEP", "RLR") && (!marker_strategy %in% c("all", "pos_fc", "top_50p_logFC", "top_50p_AveExpr"))) stop(paste0("marker_strategy must be one of ", paste0(c("all", "pos_fc", "top_50p_logFC", "top_50p_AveExpr"), collapse = ",")))
  if (decon_method == "scaden"){
    if (is.null(pythonpath) && reticulate::py_available()) {
      warning(paste0("pythonpath not specified, will use the python found here: ", reticulate::py_config()$python))
      pythonpath <- reticulate::py_config()$python
    } else if(is.null(pythonpath) && (!reticulate::py_available())) stop("pythonpath not specified, and python cannot be found via py_available()")
    reticulate::use_python(pythonpath)
    if(!reticulate::py_module_available("scaden")) stop("scaden not installed. You can install scaden via pip or conda. See installation page for details")
  }
  if( decon_method == "CIBERSORT" && (is.null(cibersortpath) || (!file.exists(cibersortpath)))) stop("cibersortpath not provided.")
  if (ncol(phenodata) < 3) stop("phenodata should contain at least the first three columns: cellid, celltype and subjectid.")
  colnames(phenodata)[1:3] <- c("cellid", "celltype", "subjectid")
  if (length(unique(phenodata$cellid)) != nrow(phenodata)) stop("values of cellid in phenodata not unique.")
  if (length(intersect(colnames(ref), phenodata$cellid)) != length(union(colnames(ref), phenodata$cellid))) stop("column names of reference data do not match with cellid of the phenodata.")
  if ((!is.null(to_remove)) && (!to_remove %in% phenodata$celltype)) stop("to_remove not present in celltype of phenodata")
  if (decon_method %in% c("CIBERSORT", "OLS", "nnls", "FARDEEP", "RLR")) decon_type <- "bulk" else decon_type <- "sc"
  if (norm_method_bulk == "TPM") {
    if(is.null(gene_length)) {
      stop("norm_method is specified as TPM, but no gene_length is provided")
    } else if(ncol(gene_length) != 2) {
      stop("gene_length need to have two columns, name of genes & length of genes")
    } else if(length(intersect(rownames(bulk), gene_length[[1]])) != length(union(rownames(bulk), gene_length[[1]]))){
      stop("gene names for gene_length are not consistent with row names of bulk")
    }
    colnames(gene_length) <- c('GeneName', 'Length')
  }
  if(!length(intersect(rownames(ref), rownames(bulk)))) stop("no overlapping genes/rownames between ref and bulk")
  if((!is.null(marker_genes)) && (ncol(marker_genes) !=2)) stop("marker_genes need to be a data.frame with two columns (cell type, gene name)")
  if((!is.null(marker_genes)) && (!length(Reduce(intersect, list(marker_genes[,2], rownames(bulk), rownames(ref)))))) stop("no overlapping genes between supplied marker gene list and ref/bulk")
  phenodata <- phenodata[match(colnames(ref), phenodata$cellid), ]
  rownames(phenodata) <- phenodata$cellid
  if (filter_ref) {
    if(verbose) message("filtering reference data.")
    lib_sizes <- colSums(ref)
    gene_names <- rownames(ref)
    mt_id <- grepl("^MT-|_MT-", gene_names, ignore.case = TRUE)
    cellstoremove <- c()
    if(length(which(mt_id))){
      mt_percent <- colSums(ref[mt_id, ]) / lib_sizes
      cellstoremove <- filter_cells(mt_percent)
    }
    rb_id <- grepl("^RPL|^RPS|_RPL|_RPS", gene_names, ignore.case = TRUE)
    if(length(which(rb_id))){
    rb_percent <- colSums(ref[rb_id, ]) / lib_sizes
    cellstoremove <- c(cellstoremove, filter_cells(rb_percent))
    }
    cellstoremove <- unique(c(cellstoremove,  filter_cells(lib_sizes)))
    if (length(cellstoremove) != 0) {
      ref <- ref[, -cellstoremove]
      phenodata <- phenodata[-cellstoremove, ]
    }
    keep <- sapply(unique(phenodata$celltype), function(i) {
      ct_hits <- which(phenodata$celltype == i)
      size <- ceiling(min_pct_ct * length(ct_hits))
      rowSums(ref[, ct_hits, drop = FALSE] != 0) >= size
    })
    ref <- ref[rowSums(keep) > 0, ]
  }
  if(!is.null(genes_to_remove)) ref <- ref[!rownames(ref) %in% genes_to_remove,]
  if(decon_type == "bulk"){
    celltypes <- phenodata$celltype
    avgexp_ct <- lapply(unique(celltypes), function(i) rowMeans(ref[, celltypes == i]))
    avgexp_ct <- do.call(cbind.data.frame, avgexp_ct)
    colnames(avgexp_ct) <- unique(celltypes)
    if(is.null(marker_genes)){
      if(verbose) message("find cell type marker genes using limma")
      keep <- sapply(unique(celltypes), function(i) {
        ct_hits <- which(celltypes == i)
        size <- ceiling(min_pct_ct * length(ct_hits))
        rowSums(ref[, ct_hits, drop = FALSE] != 0) >= size
      })
      ref_sel <- ref[rowSums(keep) > 0, ]
      ref_sel <- Normalization(ref_sel)
      annotation <- factor(celltypes)
      design <- model.matrix(~ 0 + annotation)
      colnames(design) <- gsub("annotation", "", colnames(design))
      cont_matrix <- matrix((-1 / ncol(design)), nrow = ncol(design), ncol = ncol(design))
      colnames(cont_matrix) <- colnames(design)
      diag(cont_matrix) <- (ncol(design) - 1) / ncol(design)
      tmp <- voom(ref_sel, design = design, plot = FALSE)
      fit <- lmFit(tmp, design)
      fit2 <- contrasts.fit(fit, cont_matrix)
      fit2 <- eBayes(fit2, trend = TRUE)
      markers <- markerfc(fit2, log2_threshold = lfc_markers)
      if(length(unique(markers$CT)) != length(unique(phenodata$celltype))) warning(paste0("some cell types don't have marker genes: ", paste0(setdiff(unique(phenodata$celltype), unique(markers$CT)), collapse = ',')))
    } else {
      colnames(marker_genes) <- c("CT", "gene")
      if(length(unique(marker_genes$CT)) != length(unique(phenodata$celltype))) warning(paste0("some cell types don't have marker genes: ", paste0(setdiff(unique(phenodata$celltype), unique(marker_genes$CT)), collapse = ',')))
      }
    }
  if (decon_type == "bulk") {
    if(verbose) message(paste0("perform transformation and normalization."))
    bulk <- scaling(bulk, norm_method_bulk, ffpe_artifacts = ffpe_artifacts, gene_length = gene_length)
    avgexp_ct <- scaling(avgexp_ct, norm_method_sc, ffpe_artifacts = FALSE, gene_length = gene_length)
    bulk <- transformation(bulk, trans_method_bulk)
    avgexp_ct <- transformation(avgexp_ct, trans_method_sc)
    if(is.null(marker_genes)) marker_distrib <- marker_strategies(markers, marker_strategy) else marker_distrib <- marker_genes
    if(!is.null(to_remove)){
      if(verbose) message(paste0("remove specified cell type: ", to_remove, "."))
      bulk <- bulk[,prop[to_remove,] != 0]
      avgexp_ct <- avgexp_ct[, colnames(avgexp_ct) %in% rownames(prop) & (!colnames(avgexp_ct) %in% to_remove)]
      prop <- prop[!rownames(prop) %in% to_remove, colnames(bulk)]
      marker_distrib <- marker_distrib[marker_distrib$CT %in% rownames(prop) & (marker_distrib$CT != to_remove),]
    }
    if(verbose) message(paste0("perform deconvolution analysis using ", decon_method, "."))
    results <- deconvolution(bulk = bulk, ref = avgexp_ct, decon_method = decon_method, marker_distrib = marker_distrib, cibersortpath = cibersortpath, verbose = verbose)
  } else if (decon_type == "sc") {
    if(verbose) message(paste0("perform transformation and normalization"))
    bulk <- scaling(bulk, norm_method_bulk, ffpe_artifacts = ffpe_artifacts, gene_length = gene_length)
    ref <- scaling(ref, norm_method_sc, ffpe_artifacts = FALSE, gene_length = gene_length)
    bulk <- transformation(bulk, trans_method_bulk)
    ref <- transformation(ref, trans_method_sc)
    if(!is.null(to_remove)){
      if(verbose) message(paste0("remove specified cell type: ", to_remove, "."))
      bulk <- bulk[,prop[to_remove,] != 0]
      ref <- ref[,phenodata$celltype != to_remove]
      prop <- prop[!rownames(prop) %in% to_remove, colnames(bulk)]
      phenodata <- phenodata[phenodata$celltype != to_remove,]
    }
    if(verbose) message(paste0("perform deconvolution analysis using ", decon_method, "."))
    results <- deconvolution(bulk = bulk, ref = ref, decon_method = decon_method, phenodata = phenodata, pythonpath = pythonpath, tmpdir = tmpdir, remove_tmpdir = remove_tmpdir, verbose = verbose, nsamples = nsamples)
  }
  output <- vector(mode = "list", length = 4)
  output[[1]] <- results[[1]]
  output[[2]] <- results[[2]]
  if(!is.null(prop)) output[[3]] <- prop
  if(decon_type == "bulk") output[[4]] <- marker_distrib
  names(output) <- c("prediction", "fit_error", "proportion", "markers")
  return(output)
}


#' Bar plot of cell type proportions across samples
#'
#' Bar plot of cell type proportions across samples
#'
#' @param prop a matrix or data.frame of cell proportion values with rows representing cell types, columns representing samples.
#' @param sort a logical value indicating whether to sort the samples based on cell type with highest median cell proportion across samples. Default to TRUE.
#' @param interactive a logical value indicating whether to generate interactive plot. Default to FALSE.
#' @export
#'

prop_barplot <- function(prop, sort = TRUE, interactive = FALSE){
  prop <- as.data.frame(t(prop))
  prop <-  prop %>% mutate(sampleid = rownames(prop)) %>%
  reshape2::melt(id.var = "sampleid", value.name = "ct_prop", variable.name = "ct") %>% as.data.frame()
 if(sort) {
  prop_median <- prop %>% group_by(ct) %>% summarise(median_prop = median(ct_prop))
  ct_sel <- prop_median$ct[which.max(prop_median$median_prop)]
  prop$sampleid <- factor(prop$sampleid, levels = prop$sampleid[prop$ct == ct_sel][order(prop$ct_prop[prop$ct == ct_sel],decreasing = T)])
  prop$ct <- factor(prop$ct, levels = prop_median$ct[order(prop_median$median_prop, decreasing = FALSE)])
  }
  gp <- ggplot() + geom_col(aes(x = sampleid,y = ct_prop, fill = ct), data = prop) +
  theme_classic() + labs(y='Predicted proportion', fill = "Cell types") +
  theme(axis.title = element_text(size=20,face='bold'),axis.text = element_text(size=15),axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),legend.title = element_text(size=15,face='bold'),legend.text = element_text(size=12))
  if(interactive) return(ggplotly(gp)) else return(gp)
}


deconvolution <- function(bulk, ref, decon_method, phenodata, marker_distrib, pythonpath = NULL, cibersortpath = NULL, tmpdir = NULL, remove_tmpdir = TRUE, verbose = FALSE, nsamples = 1000) {
  bulk_methods <- c("CIBERSORT", "OLS", "nnls", "FARDEEP", "RLR")
  sc_methods <- c("MuSiC", "SCDC", "scaden")

  ########## Using marker information for bulk_methods
  if (decon_method %in% bulk_methods) {
    ref <- ref[rownames(ref) %in% marker_distrib$gene, ]
    bulk <- bulk[rownames(bulk) %in% marker_distrib$gene, ]
  } else if (decon_method %in% sc_methods) {
    if (ncol(phenodata) < 3) stop("phenodata should contain at least the first three columns: cellid, celltype and subjectid.")
    colnames(phenodata)[1:3] <- c("cellid", "celltype", "subjectid")
    if (length(unique(phenodata$cellid)) != nrow(phenodata)) stop("values of cellid in phenodata not unique.")
    if (length(intersect(colnames(ref), phenodata$cellid)) != length(union(colnames(ref), phenodata$cellid))) stop("column names of reference data do not match with cellid of the phenodata.")
    phenodata <- phenodata[match(colnames(ref), phenodata$cellid), ]
    rownames(phenodata) <- phenodata$cellid
  }
  keep <- intersect(rownames(ref), rownames(bulk))
  if(verbose) message(length(keep), " number of genes overlapping between ref and bulk data.")
  ref <- ref[keep, ]
  bulk <- bulk[keep, ]
  if(decon_method %in% c("SCDC", "MuSiC")) {
    if(verbose) message("construct eset object.")
    ref_eset <- Biobase::ExpressionSet(assayData = as.matrix(ref), phenoData = Biobase::AnnotatedDataFrame(phenodata))
    bulk_eset <- Biobase::ExpressionSet(assayData = as.matrix(bulk))
  }
  ###################################
  if (decon_method == "CIBERSORT") { # without QN. By default, CIBERSORT performed QN (only) on the mixture.
    source(cibersortpath)
    results <- CIBERSORT(sig_matrix = ref, mixture_file = bulk, QN = FALSE)
    fiterror <- data.frame(sample = rownames(results), RMSE = results[, ncol(results)])
    results <- t(results[, 1:(ncol(results) - 3)])
  } else if (decon_method == "OLS") {
    results <- apply(bulk, 2, function(x) lm(x ~ as.matrix(ref))$coefficients[-1])
    results <- apply(results, 2, function(x) ifelse(x < 0, 0, x)) # explicit non-negativity constraint
    results <- apply(results, 2, function(x) x / sum(x)) # explicit STO constraint
    fiterror <- data.frame(
      sample = colnames(results),
      RMSE = sapply(1:ncol(bulk), function(i) {
        u <- sweep(ref, MARGIN = 2, results[, i], "*")
        k <- apply(u, 1, sum)
        sqrt((mean((k - bulk[, i])^2)))
      })
    )
    rownames(results) <- unlist(lapply(strsplit(rownames(results), ")"), function(x) x[2]))
  } else if (decon_method == "nnls") {
    results <- do.call(cbind.data.frame, lapply(apply(bulk, 2, function(x) nnls::nnls(as.matrix(ref), x)), function(y) y$x))
    results <- apply(results, 2, function(x) x / sum(x)) # explicit STO constraint
    rownames(results) <- colnames(ref)
    fiterror <- data.frame(sample = colnames(results), RMSE = sapply(1:ncol(bulk), function(i) {
      u <- sweep(ref, MARGIN = 2, results[, i], "*")
      k <- apply(u, 1, sum)
      sqrt((mean((k - bulk[, i])^2)))
    }))
  } else if (decon_method == "FARDEEP") {
    results <- t(FARDEEP::fardeep(ref, bulk, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
    results <- apply(results, 2, function(x) x / sum(x)) # explicit STO constraint
    fiterror <- data.frame(sample = colnames(results), RMSE = sapply(1:ncol(bulk), function(i) {
      u <- sweep(ref, MARGIN = 2, results[, i], "*")
      k <- apply(u, 1, sum)
      sqrt((mean((k - bulk[, i])^2)))
    }))
  } else if (decon_method == "RLR") { # RLR = robust linear regression
    results <- do.call(cbind.data.frame, lapply(apply(bulk, 2, function(x) MASS::rlm(x ~ as.matrix(ref), maxit = 100)), function(y) y$coefficients[-1]))
    results <- apply(results, 2, function(x) ifelse(x < 0, 0, x)) # explicit non-negativity constraint
    results <- apply(results, 2, function(x) x / sum(x)) # explicit STO constraint
    rownames(results) <- unlist(lapply(strsplit(rownames(results), ")"), function(x) x[2]))
    fiterror <- data.frame(sample = colnames(results), RMSE = sapply(1:ncol(bulk), function(i) {
      u <- sweep(ref, MARGIN = 2, results[, i], "*")
      k <- apply(u, 1, sum)
      sqrt((mean((k - bulk[, i])^2)))
    }))
  } else if (decon_method == "MuSiC") {
    results <- t(music_prop(
      bulk_eset, ref_eset, clusters = "celltype",
      markers = NULL, normalize = FALSE, samples = "subjectid",
      verbose = FALSE
    )$Est.prop.weighted)
    fiterror <- data.frame(sample = colnames(results), RMSE = sapply(1:ncol(bulk), function(i) {
      u <- sweep(ref, MARGIN = 2, results[, i], "*")
      k <- apply(u, 1, sum)
      sqrt((mean((k - bulk[, i])^2)))
    }))
  } else if (decon_method == "SCDC") {
    results <- t(SCDC::SCDC_prop(bulk.eset = bulk_eset, sc.eset = ref_eset, ct.varname = "celltype", sample = "subjectid", ct.sub = unique(as.character(phenodata$celltype)), iter.max = 200)$prop.est.mvw)
    fiterror <- data.frame(sample = colnames(results), RMSE = sapply(1:ncol(bulk), function(i) {
      u <- sweep(ref, MARGIN = 2, results[, i], "*")
      k <- apply(u, 1, sum)
      sqrt((mean((k - bulk[, i])^2)))
    }))
  } else if (decon_method == "scaden") {
    if(is.null(tmpdir)) {
      warning("tmpdir not supplied. Creating tmpdir in current working directory.")
      tmpdir <- paste0(getwd(), "/tmpdir/")
      dir.create(tmpdir)
    } else if(!dir.exists(tmpdir)) stop("tmpdir does not exist")
    else if(length(list.files(tmpdir)) > 0) stop("tmpdir exists, but is not empty.")
    cwd <- getwd()
    on.exit(setwd(cwd))
    setwd(tmpdir)
    data.table::fwrite(as.matrix(Matrix::t(ref)), "./decon_counts.txt", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
    data.table::fwrite(data.frame(Celltype = phenodata$celltype), "./decon_celltypes.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    data.table::fwrite(as.data.frame(bulk), "./decon_bulk_data.txt", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
    if(verbose) message("running scaden, and it might take a while. Set verbose = TRUE to track progress. ")
    system(paste0("scaden simulate --data ./ --pattern '*_counts.txt'", " -n ", nsamples), ignore.stdout = !verbose, ignore.stderr = !verbose)
    system(paste0("scaden process data.h5ad decon_bulk_data.txt"), ignore.stdout = !verbose, ignore.stderr = !verbose)
    system(paste0("scaden train processed.h5ad --steps 5000 --model_dir model"), ignore.stdout = !verbose, ignore.stderr = !verbose)
    system(paste0("scaden predict --model_dir model decon_bulk_data.txt"), ignore.stdout = !verbose, ignore.stderr = !verbose)
    results <- t(read.delim("./scaden_predictions.txt", header = TRUE, row.names = 1))
    fiterror <- NA
    if(remove_tmpdir) system(paste0("rm -rf ", tmpdir))
  }
  return(list(results, fiterror))
}
