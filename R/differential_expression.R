#' Differential expression analysis
#'
#' Performing differential expression analysis adjusting for cell proportion differences, and other covariates using a additive model.
#'
#' @param bulk a matrix-like object of gene expression values with rows representing genes, columns representing samples
#' @param prop a matrix-like object of cell proportion values with rows representing cell types, columns representing samples. Default to NULL,
#' not adjust for cell proportion.
#' @param sampleinfo a data.frame of metadata for the samples. Rows represents samples; columns represents covariates to adjust for. The first column of
#' sampleinfo should contains group information for differential analysis.
#' @param control a character value indicating the control group in sampleinfo. Set to NULL to perform ANOVA-like analysis.
#' @param case a character value indicating the case group in sampleinfo. Set to NULL to perform ANOVA-like analysis.
#' @param de_method a character value indicating the method to use for testing differential expression. Should be one of "edgeR", "DESeq2", "limma_voom", "limma"
#' @param padj_method method for adjusting multiple hypothesis testing. Default to "BH". See \code{\link{p.adjust}} for more details.
#' @param ... parameters pass to DE methods.
#'
#' @return a list of two elements:
#' \enumerate{
#'  \item a data.frame containing normalized gene expression data.
#'  \item a data.frame containing detailed differential expression statistics. Columns represent "Gene name", "log2 fold change", "log2 average expression",
#' "p value", "adjusted p value" respectively
#' }
#'
#' @details
#' To perform ANOVA like analysis (differences between any groups), set \code{control} & \code{case} options to \code{NULL} and choose one of the following methods:
#' edgeR, limma_voom or limma. DESeq2 does not provide direct support for this type of comparison.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ref_list <- c(paste0(system.file("extdata", package = "SCdeconR"), "/refdata/sample1"),
#'               paste0(system.file("extdata", package = "SCdeconR"), "/refdata/sample2"))
#' phenopath1 <- paste0(system.file("extdata", package = "SCdeconR"),
#' "/refdata/phenodata_sample1.txt")
#' phenopath2 <- paste0(system.file("extdata", package = "SCdeconR"),
#' "/refdata/phenodata_sample2.txt")
#' phenodata_list <- c(phenopath1,phenopath2)
#'
#' # construct integrated reference using harmony algorithm
#' refdata <- construct_ref(ref_list = ref_list,
#'                       phenodata_list = phenodata_list,
#'                       data_type = "cellranger",
#'                       method = "harmony",
#'                       group_var = "subjectid",
#'                       nfeature_rna = 50,
#'                       vars_to_regress = "percent_mt", verbose = FALSE)
#' phenodata <- data.frame(cellid = colnames(refdata),
#'                         celltypes = refdata$celltype,
#'                         subjectid = refdata$subjectid)
#' ## construct a vector with same proportions across cell types
#' prop1 <- data.frame(celltypes = unique(refdata$celltype),
#'                    proportion = rep(0.125, 8))
#' ## simulate 20 bulk samples based on specified cell type proportion
#' bulk_sim1 <- bulk_generator(ref = GetAssayData(refdata, slot = "data", assay = "SCT"),
#'                             phenodata = phenodata,
#'                             num_mixtures = 20,
#'                             prop = prop1, replace = TRUE)
#' ## generate another vector with high proportion for a certian cell type
#' prop2 <- data.frame(celltypes = unique(refdata$celltype),
#'                     proportion = c(0.8, 0.1, 0.1, rep(0, 5)))
#' bulk_sim2 <- bulk_generator(ref = GetAssayData(refdata, slot = "data", assay = "SCT"),
#'                             phenodata = phenodata,
#'                             num_mixtures = 20,
#'                             prop = prop2, replace = TRUE)
#' ## compare data for differential analysis
#' bulk_sim <- list(cbind(bulk_sim1[[1]], bulk_sim2[[1]]),
#'                  cbind(bulk_sim1[[2]], bulk_sim2[[2]]))
#' ## force to be integer for DE purposes
#' bulk <- round(bulk_sim[[1]], digits=0)
#' colnames(bulk) <- paste0("sample", 1:ncol(bulk))
#' ## predict cell type proportions using "OLS" algorithm
#' decon_res <- scdecon(bulk = bulk,
#'                      ref = GetAssayData(refdata, slot = "data", assay = "SCT"),
#'                      phenodata = phenodata,
#'                      filter_ref = TRUE,
#'                      decon_method = "OLS",
#'                      norm_method = "none",
#'                      trans_method = "none",
#'                      marker_strategy = "all")
#' ## create sampleinfo
#' sampleinfo <- data.frame(condition = rep(c("group1", "group2"), each =20))
#' row.names(sampleinfo) <- colnames(bulk)
#' library(DESeq2)
#' deres <- run_de(bulk = bulk,
#'                prop = decon_res[[1]],
#'                sampleinfo = sampleinfo,
#'                control = "group1",
#'                case = "group2",
#'                de_method = "edgeR")
#'
#'  ## run differential analysis without adjusting for cell proportion differences
#'  deres_notadjust <- run_de(bulk = bulk,
#'                            prop = NULL,
#'                            sampleinfo = sampleinfo,
#'                            control = "group1",
#'                            case = "group2",
#'                            de_method = "edgeR")
#'
#' ## scatter plot to compare the effect of adjusting cell proportion differences
#' comparedeg_scatter(results1 = deres[[2]],
#'                    results2 = deres_notadjust[[2]],
#'                    result_names = c("adjust for cell proportion", "not adjust for cell proportion"),
#'                    fc_cutoff = 1.5,
#'                    pval_cutoff = 0.05,
#'                    pvalflag = TRUE,
#'                    interactive = T)
#'
#' ## generate cell-type specific gene expression
#' ct_exprs_list <- celltype_expression(bulk = bulk,
#'                                      ref = GetAssayData(refdata, slot = "data", assay = "SCT"),
#'                                      phenodata = phenodata,
#'                                      prop = decon_res[[1]],
#'                                      UMI_min = 0,
#'                                      CELL_MIN_INSTANCE = 1)
#'
#'
#' }


run_de <- function(
    bulk,
    prop = NULL,
    sampleinfo,
    control = NULL,
    case = NULL,
    de_method = c("edgeR", "DESeq2", "limma_voom", "limma"),
    padj_method = "BH", ...) {
    if (length(intersect(colnames(bulk), rownames(sampleinfo))) != length(union(colnames(bulk), rownames(sampleinfo)))) stop("column names of bulk do not match with row names of sampleinfo.")
    colnames(sampleinfo)[1] <- "group"
    sampleinfo <- sampleinfo[match(colnames(bulk), rownames(sampleinfo)),,drop = FALSE]
    if (!is.null(prop)) {
        if (length(intersect(colnames(bulk), colnames(prop))) != length(union(colnames(bulk), colnames(prop)))) stop("column names of bulk do not match with column names of prop.")
        prop <- prop[, match(colnames(bulk), colnames(prop))]
    }
    if (is.null(control) && (!is.null(case)) || (!is.null(control)) && (is.null(case))) {
        stop("You need to either 1) set values for both control and case, or 2) leave both control and case as NULL for anova-like analysis")
    }
    if ((!is.null(control)) && (!is.null(case)) && ((!control %in% sampleinfo$group) | (!case %in% sampleinfo$group))) {
        stop("control and case defined, but not not found in the first column of sampleinfo")
    }
    if (length(unique(sampleinfo$group)) == 2 && is.null(control) && is.null(case)) {
        stop("Only two conditions found within first column of sampleinfo. You need to set values for both control and case")
    }
    if (de_method == "limma" && !is.null(prop)) stop("limma is recommended for testing cell-type specific gene expression, but prop is provided. Use celltype_expression first.")
    if (!de_method %in% c(c("edgeR", "DESeq2", "limma_voom", "limma"))) stop("de_method has to be one of edgeR, DESeq2, limma_voom or limma")
    if (de_method == "edgeR") {
        results <- edger_fun(bulk, prop, sampleinfo, control, case, padj_method, ...)
    } else if (de_method == "DESeq2") {
        if (length(unique(sampleinfo$group)) > 2) stop("DE analysis with group > 2 is not directly supported by DESeq2") else results <- deseq2_fun(bulk, prop, sampleinfo, control, case, padj_method, ...)
    } else if (de_method == "limma_voom"){
        results <- limma_voom_fun(bulk, prop, sampleinfo, control, case, padj_method, ...)
    } else {
        results <- limma_fun(bulk, sampleinfo, control, case, padj_method, ...)
    }
}


#' Generate a scatter plot comparing two differential expression results
#'
#' Generate a scatter plot of fold changes comparing two differential expression results, e.g. w/wo adjusting for cell proportion differences.
#'
#' @param results1 a data.frame containing differential expression results with five columns: "Gene name", "log2 fold change", "log2 average expression",
#' "p value", "adjusted p value". The second element of the output from function \code{\link{run_de}}.
#' @param results2 similar to \code{results1}.
#' @param result_names a vector of length 2 indicating the names of the two differential results. If NULL, names will be set to c("results1", "results2")
#' @param fc_cutoff fold change cutoff to identify differential expressed genes.
#' @param pval_cutoff p value cutoff to identify differential expressed genes.
#' @param pvalflag a logical value indicating whether to use adjusted p value in selecting differential expressed genes.
#' @param interactive a logical value indicating whether to generate an interactive plot.
#'
#' @details See examples from \code{\link{run_de}}.
#' @export
#'
#'

comparedeg_scatter <- function(
    results1,
    results2,
    result_names = NULL,
    fc_cutoff,
    pval_cutoff,
    pvalflag = TRUE,
    interactive = FALSE) {
    if (ncol(results1) != 5 || ncol(results2) != 5) stop("inputs need to have five columns: genename, log2foldchange, log2avgexp, pval & padj")
    colnames(results1) <- colnames(results2) <- c("genename", "log2foldchange", "log2avgexp", "pval", "padj")
    if ((!is.null(result_names)) && length(result_names) != 2) stop("The length of result_names has to equal to 2") else if (is.null(result_names)) result_names <- c("result1", "result2")
    results_combined <- inner_join(results1, results2, by = "genename")
    if (pvalflag) {
        idx1 <- which(abs(results_combined$log2foldchange.x) >= log2(fc_cutoff) & results_combined$padj.x <= pval_cutoff)
        idx2 <- which(abs(results_combined$log2foldchange.y) >= log2(fc_cutoff) & results_combined$padj.y <= pval_cutoff)
    } else {
        idx1 <- which(abs(results_combined$log2foldchange.x) >= log2(fc_cutoff) & results_combined$pval.x <= pval_cutoff)
        idx2 <- which(abs(results_combined$log2foldchange.y) >= log2(fc_cutoff) & results_combined$pval.y <= pval_cutoff)
    }
    if (length(idx1) == 0 && length(idx2) == 0) stop("No DEGs identified, try loosing the thresholds")
    results_combined$category <- "Not significant"
    results_combined$category[intersect(idx1, idx2)] <- "Both"
    results_combined$category[setdiff(idx1, idx2)] <- paste0(result_names[1], " only")
    results_combined$category[setdiff(idx2, idx1)] <- paste0(result_names[2], " only")
    color_df <- data.frame(
        group = c("Both", paste0(result_names[1], " only"), paste0(result_names[2], " only"), "Not significant"),
        colors = c("red", "seagreen1", "royalblue", "grey")
    )
    results_combined$category <- factor(results_combined$category, levels = color_df$group[color_df$group %in% unique(results_combined$category)])
    lim_values <- range(results_combined[, c("log2foldchange.x", "log2foldchange.y")], na.rm = TRUE)
    options(warn = -1)
    gp <- ggplot() +
        geom_point(aes(x = log2foldchange.x, y = log2foldchange.y, color = category, text = paste0("GeneSymbol: ", genename)), size = 5, data = results_combined) +
        theme_classic() + labs(x= result_names[1], y = result_names[2]) +
        scale_color_manual(values = color_df$colors[color_df$group %in% unique(results_combined$category)]) +
        lims(x = lim_values, y = lim_values) +
        geom_hline(yintercept = c(log2(fc_cutoff), -log2(fc_cutoff)), linetype = 2) +
        geom_vline(xintercept = c(log2(fc_cutoff), -log2(fc_cutoff)), linetype = 2) +
        geom_abline(slope = 1, intercept = 0)
    options(warn = 0)
    if (interactive) {
        return(plotly::ggplotly(gp))
    } else {
        return(gp)
    }
}

#' Compute cell type specific gene expression
#'
#' Compute cell type specific gene expression based on predicted cell proportions and reference data.
#'
#' @param bulk a matrix-like object of bulk RNA-seq data with rows representing genes, columns representing samples
#' @param ref a matrix-like object of scRNA-seq data with rows representing genes, columns representing cells.
#' @param phenodata a data.frame with rows representing cells, columns representing cell attributes. It should at least contain the first two
#' columns as:
#' \enumerate{
#'  \item cell barcodes
#'  \item cell types
#' }
#' @param prop a matrix-like object of cell proportion values with rows representing cell types, columns representing samples.
#' @param ... additional parameters passed to \code{create.RCTD} from \code{spacexr}.
#'
#' @return a list with length equal to number of unique cell types in phenodata. Each element in the list represents gene expression matrix for each unique cell type.
#'
#' @details this function is inspired by cell-type specific gene expression estimation for doublet mode in \code{spacexr}. See examples from \code{\link{run_de}}.
#'
#' @export
#'

celltype_expression <- function(bulk, ref, phenodata, prop, ...) {
    if (ncol(phenodata) < 2) stop("phenodata should contain at least the first two columns: cellid, celltype.")
    colnames(phenodata)[1:2] <- c("cellid", "celltype")
    if (length(unique(phenodata$cellid)) != nrow(phenodata)) stop("values of cellid in phenodata not unique.")
    if (length(intersect(colnames(ref), phenodata$cellid)) != length(union(colnames(ref), phenodata$cellid))) stop("column names of reference data do not match with cellid of the phenodata.")
    if (length(union(colnames(bulk), colnames(prop))) != length(intersect(colnames(bulk), colnames(prop)))) stop("observed inconsistency between column names of bulk and prop")
    prop <- prop[, match(colnames(bulk), colnames(prop))]
    phenodata <- phenodata[match(colnames(ref), phenodata$cellid), ]
    celltypes <- factor(phenodata[, "celltype"])
    names(celltypes) <- phenodata$cellid
    ### use functions from spacexr package
    reference <- eval(parse(text = 'spacexr::Reference(ref, celltypes, require_int = FALSE, min_UMI = 0)'))
    bulk_fake_spatial <- eval(parse(text = 'spacexr::SpatialRNA(counts = bulk, use_fake_coords = TRUE, require_int = FALSE)'))
    rctd_obj <- eval(parse(text = 'spacexr::create.RCTD(bulk_fake_spatial, reference, fc_cutoff = 0, ...)'))
    rctd_obj <- eval(parse(text = 'spacexr::fitBulk(rctd_obj)'))
    celltypes <- rctd_obj@cell_type_info$info[[2]]
    ct_renorm <- rctd_obj@cell_type_info$renorm[[1]]
    ct_mat <- lapply(1:length(celltypes), function(i) {
        mat1 <- matrix(ct_renorm[, colnames(ct_renorm) == celltypes[i]], nrow = nrow(ct_renorm), ncol = 1)
        mat2 <- matrix(prop[rownames(prop) == celltypes[i], ], nrow = 1, ncol = ncol(prop))
        mat1 %*% mat2
    })
    ct_sum <- Reduce("+", ct_mat)
    ct_exp <- vector("list", length(celltypes))
    names(ct_exp) <- celltypes
    for (i in 1:length(celltypes)) ct_exp[[i]] <- rctd_obj@spatialRNA@counts * (ct_mat[[i]] / ct_sum)
    return(ct_exp)
}


deseq2_fun <- function(counts, prop, sampleinfo, control, case, padj_method, ...) {
    message("Carrying out differential expression analysis using DESeq2!")
    if (!is.null(prop)) {
        ## remove one cell-type to make design matrix full rank
        prop <- prop[-which.min(apply(prop, 1, median)), ]
        sampleinfo <- cbind(sampleinfo, t(prop))
    }
    colnames(sampleinfo) <- make.names(colnames(sampleinfo))
    formula_de <- as.formula(paste0("~ ", paste0(colnames(sampleinfo), collapse = "+")))
    dse <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = sampleinfo, design = formula_de)
    dse <- DESeq2::estimateSizeFactors(dse)
    normdata <- DESeq2::counts(dse, normalized = TRUE)
    dse <- DESeq2::estimateDispersions(dse, ...)
    dse <- DESeq2::nbinomWaldTest(dse)
    res <- DESeq2::results(dse, cooksCutoff = FALSE, contrast = c("group", case, control))
    res <- as.data.frame(res)
    teststats <- res %>%
        select(-c(lfcSE, stat, padj)) %>%
        dplyr::rename(log2foldchange = log2FoldChange, pval = pvalue) %>%
        mutate(log2avgexp = log2(baseMean + 0.1), .after = log2foldchange) %>%
        mutate(padj = p.adjust(pval, method = padj_method)) %>%
        mutate(genename = rownames(res), .before = log2foldchange) %>%
        select(-baseMean)
    list(normdata, teststats)
}


edger_fun <- function(counts, prop, sampleinfo, control, case, padj_method, ...) {
    message("Carrying out differential expression analysis using edgeR!")
    cds <- DGEList(counts = counts, group = sampleinfo$group)
    cds <- calcNormFactors(cds)
    if (!is.null(prop)) {
        ## remove one cell-type to make design matrix full rank
        prop <- prop[-which.min(apply(prop, 1, median)), ]
        sampleinfo <- cbind(sampleinfo, t(prop))
    }
    colnames(sampleinfo) <- make.names(colnames(sampleinfo))
    if (is.null(control) && is.null(case)) {
        formula_de <- as.formula(paste0("~ ", paste0(colnames(sampleinfo), collapse = "+")))
        design <- model.matrix(formula_de, data = sampleinfo)
        colnames(design)[2:length(unique(sampleinfo$group))] <- gsub("^group", "", colnames(design)[2:length(unique(sampleinfo$group))])
        cds <- estimateGLMRobustDisp(cds, design, ...)
        fit <- glmQLFit(cds, design)
        lrt <- glmQLFTest(fit, coef = 2:length(unique(sampleinfo$group)))
    } else {
        formula_de <- as.formula(paste0("~ 0 + ", paste0(colnames(sampleinfo), collapse = "+")))
        design <- model.matrix(formula_de, data = sampleinfo)
        colnames(design)[1:length(unique(sampleinfo$group))] <- gsub("^group", "", colnames(design)[1:length(unique(sampleinfo$group))])
        cds <- estimateGLMRobustDisp(cds, design, ...)
        fit <- glmQLFit(cds, design)
        lrt <- eval(parse(text = paste0("glmQLFTest(fit,contrast=makeContrasts(",case,"-",control,",levels=design))")))
    }
    normdata <- cpm(cds)
    teststats <- lrt$table %>%
        select(-F) %>%
        dplyr::rename(logavgexp = logCPM, pval = PValue) %>%
        rename_with(~ gsub("logFC", "log2foldchange", .x)) %>%
        mutate(padj = p.adjust(pval, method = padj_method)) %>%
        mutate(genename = rownames(lrt$table), .before = log2foldchange)
    return(list(normdata, teststats))
}


limma_voom_fun <- function(counts, prop, sampleinfo, control, case, padj_method, ...) {
    message(paste0("Carrying out differential expression analysis using limma_voom!", "\n"))
    if (!is.null(prop)) {
        ## remove one cell-type to make design matrix full rank
        prop <- prop[-which.min(apply(prop, 1, median)), ]
        sampleinfo <- cbind(sampleinfo, t(prop))
    }
    colnames(sampleinfo) <- make.names(colnames(sampleinfo))
    nf <- calcNormFactors(counts)
    y <- limma::voom(counts, plot = FALSE, lib.size = colSums(counts) * nf)
    normdata <- as.matrix(2^(y$E))
    if (is.null(control) && is.null(case)) {
        formula_de <- as.formula(paste0("~ ", paste0(colnames(sampleinfo), collapse = "+")))
        design <- model.matrix(formula_de, data = sampleinfo)
        colnames(design)[2:length(unique(sampleinfo$group))] <- gsub("^group", "", colnames(design)[2:length(unique(sampleinfo$group))])
        fit <- limma::lmFit(y, design, ...)
        fit <- limma::eBayes(fit)
        lmres <- limma::topTable(fit, coef = 2:length(unique(sampleinfo$group)), n = nrow(counts), sort.by = "none")
        ncoef <- length(unique(sampleinfo$group)) - 1
    } else {
        formula_de <- as.formula(paste0("~ 0 + ", paste0(colnames(sampleinfo), collapse = "+")))
        design <- model.matrix(formula_de, data = sampleinfo)
        colnames(design)[1:length(unique(sampleinfo$group))] <- gsub("^group", "", colnames(design)[1:length(unique(sampleinfo$group))])
        fit <- limma::lmFit(y, design, ...)
        fit2 <- eval(parse(text = paste0("limma::contrasts.fit(fit, contrasts = limma::makeContrasts(", case,"-",control,",levels=design))")))
        fit2 <- limma::eBayes(fit2)
        lmres <- limma::topTable(fit2, coef = 1, n = nrow(counts), sort.by = "none")
        ncoef <- 1
    }
    teststats <- lmres %>%
        select(-c(F, adj.P.Val)) %>%
        dplyr::rename(pval = P.Value, log2avgexp = AveExpr) %>%
        rename_with(~ gsub("^", "log2foldchange.", .x), .cols = 1:all_of(ncoef)) %>%
        mutate(padj = p.adjust(pval, method = padj_method)) %>%
        mutate(genename = rownames(lmres), .before = log2foldchange)
    return(list(normdata, teststats))
}


limma_fun <- function(counts, prop, sampleinfo, control, case, p_adj_method, ...) {
    message(paste0("Carrying out differential expression analysis using limma!", "\n"))
    normdata <- preprocessCore::normalize.quantiles(as.matrix(counts))
    rownames(normdata) <- rownames(counts)
    colnames(normdata) <- colnames(counts)
    normdata_log2 <- log2(normdata + 0.01)
    colnames(sampleinfo) <- make.names(colnames(sampleinfo))
    if (is.null(control) && is.null(case)) {
        formula_de <- as.formula(paste0("~ ", paste0(colnames(sampleinfo), collapse = "+")))
        design <- model.matrix(formula_de, data = sampleinfo)
        colnames(design)[2:length(unique(sampleinfo$group))] <- gsub("^group", "", colnames(design)[2:length(unique(sampleinfo$group))])
        fit <- limma::lmFit(normdata_log2, design, ...)
        fit <- limma::eBayes(fit)
        lmres <- limma::topTable(fit, coef = 2:length(unique(sampleinfo$group)), n = nrow(counts), sort.by = "none")
        ncoef <- length(unique(sampleinfo$group)) - 1
    } else {
        formula_de <- as.formula(paste0("~ 0 + ", paste0(colnames(sampleinfo), collapse = "+")))
        design <- model.matrix(formula_de, data = sampleinfo)
        colnames(design)[1:length(unique(sampleinfo$group))] <- gsub("^group", "", colnames(design)[1:length(unique(sampleinfo$group))])
        fit <- limma::lmFit(normdata_log2, design, ...)
        fit2 <- eval(parse(text = paste0("limma::contrasts.fit(fit, contrasts = limma::makeContrasts(", case,"-",control,",levels=design))")))
        fit2 <- limma::eBayes(fit2)
        lmres <- limma::topTable(fit2, coef = 1, n = nrow(counts), sort.by = "none")
        ncoef <- 1
    }
    teststats <- lmres %>%
        select(-c(F, adj.P.Val)) %>%
        dplyr::rename(pval = P.Value, log2avgexp = AveExpr) %>%
        rename_with(~ gsub("^", "log2foldchange.", .x), .cols = 1:all_of(ncoef)) %>%
        mutate(padj = p.adjust(pval, method = padj_method)) %>%
        mutate(genename = rownames(lmres), .before = log2foldchange)
    return(list(normdata, teststats))
}
