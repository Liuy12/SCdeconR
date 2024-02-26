filter_cells <- function(filter_param) {
  cellstoremove <- which(filter_param > median(filter_param) + 3 * mad(filter_param) | filter_param < median(filter_param) - 3 * mad(filter_param))
  cellstoremove
}

Normalization <- function(data) {
  data <- DGEList(data)
  CtrlGenes <- grep("ERCC-", rownames(data))

  if (length(CtrlGenes) > 1) {
    spikes <- data[CtrlGenes, ]
    spikes <- calcNormFactors(spikes, method = "TMM")
    data$samples$norm.factors <- spikes$samples$norm.factors
  } else {
    data <- calcNormFactors(data, method = "TMM")
  }

  return(data)
}

markerfc <- function(fit2, log2_threshold = log2(1.5), output_name = "markers") {
  topTable_RESULTS <- limma::topTable(fit2, coef = 1:ncol(fit2$contrasts), number = Inf, adjust.method = "BH", p.value = 0.05, lfc = log2_threshold)
  AveExpr_pval <- topTable_RESULTS[, (ncol(topTable_RESULTS) - 3):ncol(topTable_RESULTS)]
  topTable_RESULTS <- topTable_RESULTS[, 1:(ncol(topTable_RESULTS) - 4)]

  if (length(grep("ERCC-", topTable_RESULTS$gene)) > 0) {
    topTable_RESULTS <- topTable_RESULTS[-grep("ERCC-", topTable_RESULTS$gene), ]
  }

  markers <- apply(topTable_RESULTS, 1, function(x) {
    temp <- sort(x)
    ((temp[ncol(topTable_RESULTS)] - temp[ncol(topTable_RESULTS) - 1]) >= log2_threshold) | (abs(temp[1] - temp[2]) >= log2_threshold)
  })

  topTable_RESULTS <- topTable_RESULTS[markers, ]

  markers <- cbind.data.frame(
    rownames(topTable_RESULTS),
    t(apply(topTable_RESULTS, 1, function(x) {
      temp <- max(x)
      if (temp < log2_threshold) {
        temp <- c(min(x), colnames(topTable_RESULTS)[which.min(x)])
      } else {
        temp <- c(max(x), colnames(topTable_RESULTS)[which.max(x)])
      }
      temp
    }))
  )

  colnames(markers) <- c("gene", "log2FC", "CT")
  markers$log2FC <- as.numeric(as.character(markers$log2FC))
  markers <- markers %>% arrange(CT, desc(log2FC))

  markers$AveExpr <- AveExpr_pval$AveExpr[match(markers$gene, rownames(AveExpr_pval))]
  markers$gene <- as.character(markers$gene)
  markers$CT <- as.character(markers$CT)

  return(markers)
}



marker_strategies <- function(marker_distrib, marker_strategy) {
  if (marker_strategy == "all") {

    # using all markers that were found
    markers <- marker_distrib
  } else if (marker_strategy == "pos_fc") {

    # using only markers with positive FC (=over-expressed in cell type of interest)
    markers <- marker_distrib %>%
      filter(log2FC > 0) %>%
      as.data.frame()
  } else if (marker_strategy == "top_50p_logFC") {

    # top 50% of markers (per CT) based on logFC
    markers <- marker_distrib %>%
      filter(log2FC > 0) %>%
      arrange(CT, desc(log2FC)) %>%
      group_by(CT) %>%
      top_n(ceiling(n() * 0.5), wt = log2FC) %>%
      as.data.frame()
  } else if (marker_strategy == "top_50p_AveExpr") {

    # top 50% of markers based on average gene expression (baseline expression)
    markers <- marker_distrib %>%
      filter(log2FC > 0) %>%
      arrange(CT, desc(AveExpr)) %>%
      group_by(CT) %>%
      top_n(ceiling(n() * 0.5), wt = log2FC) %>%
      as.data.frame()
  }

  return(markers)
}

#' Transformation of gene expression data
#'
#' Methods to use for data transformation.
#'
#' @param matrix a matrix-like objector of gene expression values with rows representing genes, columns representing samples or cells
#' @param option character value specifying the transformation method to use. Has to be one of "none", "log", "sqrt", "vst".
#' @return a matrix-like object with the same dimension of input object after data transformation.
#'
#' @details
#' refer to \code{\link{scdecon}} for more details.
#' @export


transformation <- function(matrix, option) {

  #############################################################
  ##########    DATA TRANSFORMATION (on full data)   ##########
  if (option == "none") {
    matrix <- matrix
  }

  if (option == "log2") {
    matrix <- log2(matrix + 1)
  }

  if (option == "sqrt") {
    matrix <- sqrt(matrix)
  }

  if (option == "vst") {
    matrix <- DESeq2::varianceStabilizingTransformation(as.matrix(matrix))
  }

  return(matrix)
}


#' Normalization of gene expression data
#'
#' Methods to use for data normalization.
#'
#' @param matrix a matrix-like objector of gene expression values with rows representing genes, columns representing samples or cells
#' @param option character value specifying the normalization method to use. Has to be one of "none", "LogNormalize", "TMM", "median_ratios", "TPM",
#' "SCTransform", "scran", "scater", "Linnorm".
#' @param gene_length a data.frame with two columns. The first column represents gene names that match with provided bulk data. The second column
#' represents length of each gene. Only applicable when norm_method is selected as "TPM"
#' @param seed random seed used for simulating FFPE artifacts. Only applicable when ffpe_artifacts is set to TRUE.
#' @param ffpe_artifacts logical value indicating whether to add simulated ffpe artifacts in the bulk data. Only applicable to simulation experiments in
#' evaluating the effect of FFPE artifacts.
#' @param model pre-constructed ffpe model data. Can be downloaded from github: https://github.com/Liuy12/SCdeconR/blob/master/data/ffpemodel.rda
#'
#' @return a matrix-like object with the same dimension of input object after data normalization.
#'
#' @details
#' refer to \code{\link{scdecon}} for more details.
#' @export


scaling <- function(matrix, option, gene_length = NULL, seed = 1234, ffpe_artifacts = FALSE, model = NULL) {

  ##########    Remove rows & columns full of zeroes   ##########
  matrix <- matrix[rowSums(matrix) != 0, ]
  # OLS with error if all elements within a row are equal (e.g. all 0, or all a common value after log/sqrt/vst transformation)
  matrix <- matrix[!apply(matrix, 1, function(x) var(x) == 0), ]
  matrix <- matrix[, colSums(matrix) != 0]

  if (option == "LogNormalize") {
    matrix <- as.matrix(expm1(LogNormalize(matrix, verbose = FALSE))) # for v3
  } else if (option == "TMM") {
    matrix <- DGEList(counts = matrix)
    CtrlGenes <- grep("ERCC-", rownames(data))

    if (length(CtrlGenes) > 1) {
      spikes <- data[CtrlGenes, ]
      spikes <- calcNormFactors(spikes, method = "TMM")
      matrix$samples$norm.factors <- spikes$samples$norm.factors
    } else {
      matrix <- calcNormFactors(matrix, method = "TMM")
    }

    matrix <- cpm(matrix)
  } else if (option == "median_ratios") { # requires integer values
    metadata <- data.frame(sampleid = colnames(matrix))
    CtrlGenes <- grep("ERCC-", rownames(matrix))
    matrix <- DESeq2::DESeqDataSetFromMatrix(matrix, colData = metadata, design = ~1)
    if (length(CtrlGenes) > 1 & sum(rowSums(counts(matrix[CtrlGenes, ]) != 0) >= 0.5 * (ncol(matrix))) >= 2) {
      dds <- DESeq2::estimateSizeFactors(matrix, type = "ratio", controlGenes = CtrlGenes)
    } else {
      dds <- DESeq2::estimateSizeFactors(matrix, type = "ratio")
    }
    matrix <- DESeq2::counts(dds, normalized = TRUE)
  } else if (option == "TPM") {
    matrix <- matrix[rownames(matrix) %in% gene_length$GeneName, ]
    matrix <- matrix %>%
      DGEList() %>%
      calcNormFactors() %>%
      rpkm(gene.length = gene_length$Length[match(rownames(matrix), gene_length$GeneName)])
    matrix <- sweep(matrix, 2, apply(matrix, 2, sum), "/") * 10^6
    rownames(matrix) <- toupper(rownames(matrix))
    if (ffpe_artifacts) {
      set.seed(seed)
      ### convert to log2 scale
      matrix <- log2(matrix + 0.1)
      ### load GAM fitted model
      #model <- data("ffpemodel")
      diffmat <- matrix(
        sample(c(-1, 1), nrow(matrix) * ncol(matrix), replace = TRUE, prob = c(0.5, 0.5)) *
          rnorm(1:length(c(matrix)),
            mean = predict.gam(model, newdata = data.frame(avgtpm = c(matrix))), sd = 0.5
          ),
        nrow = nrow(matrix), ncol = ncol(matrix), byrow = FALSE
      )
      print(diffmat[1:5, 1:5])
      matrix <- 2^(matrix - diffmat)
    }
    ####################################################################################
    ## scRNA-seq specific
  } else if (option == "SCTransform") {
    matrix <- as(matrix, "dgCMatrix")
    matrix <- sctransform::vst(matrix, return_corrected_umi = TRUE, show_progress = FALSE)$umi_corrected
    matrix <- as(matrix, "matrix")
  } else if (option == "scran") {
    sf <- scran::computeSumFactors(as.matrix(matrix), clusters = NULL)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(matrix)))
    sizeFactors(sce) <- sf
    sce <- scater::normalize(sce, exprs_values = "counts", return_log = FALSE)
    matrix <- SingleCellExperiment::normcounts(sce)
  } else if (option == "scater") {
    size_factors <- scater::librarySizeFactors(matrix)
    matrix <- scater::normalizeCounts(as.matrix(matrix), size_factors = size_factors, return_log = FALSE)
  } else if (option == "Linnorm") { # It is not compatible with log transformed datasets.
    ## use Linnorm function from Linnorm package
    matrix <- expm1(Linnorm::Linnorm(as.matrix(matrix))) # Main function contains log1p(datamatrix)
  }
  return(matrix)
}
