
#' Generate artificial bulk RNA-seq samples based on simulation
#'
#' Generate artificial bulk RNA-seq samples with random or pre-defined cell-type proportions for benchmarking deconvolution algorithms
#'
#' @param ref a matrix-like object of gene expression values with rows representing genes, columns representing cells.
#' @param phenodata a data.frame with rows representing cells, columns representing cell attributes. It should at least contain the first two columns as:
#' \enumerate{
#'  \item cell barcodes
#'  \item cell types
#' }
#' @param num_mixtures total number of simulated bulk samples. Have to be multiple of \code{num_mixtures_sprop}. Default to 500.
#' @param num_mixtures_sprop number of simulated bulk samples with the same simulated cell type proportions. Only applicable when \code{prop} is not specified.
#' Those samples will be used to estimate bias & variance. Default to 10.
#' @param pool_size number of cells to use to construct each artificial bulk sample. Default to 100.
#' @param seed seed to use for simulation. Default to 1234.
#' @param prop a data.frame with two columns. The first column includes unique cell types in phenodata; the second column includes cell type proportions.
#' If specified, bulk samples will be simulated based on the specified cell proportions.
#' @param replace logical value indicating whether to sample cells with replacement. Default to FALSE, to sample cells without replacement.
#'
#' @return a list of two objects:
#' \enumerate{
#'  \item simulated bulk RNA-seq data, with rows representing genes, columns representing samples
#'  \item cell type proportions used to simulate the bulk RNA-seq data, with rows representing cell types, columns representing samples
#' }
#'
#' @details
#' If \code{prop} is not specified, cell type proportions will be firstly randomly generated with at least two cell types present. Then, for each cell proportion
#' vector, \code{num_mixtures_sprop} number of samples is simulated. Eventually, a total of \code{num_mixtures} number of samples is simulated. If prop is
#' specified, then a total of \code{num_mixtures} number of samples will be simulated based on the same cell proportion vector specified.
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
#' bulk_sim <- bulk_generator(ref = GetAssayData(refdata, slot = "data", assay = "SCT"),
#'                            phenodata = phenodata,
#'                            num_mixtures = 20,
#'                            num_mixtures_sprop = 1)
#' }


bulk_generator <- function(
    ref,
    phenodata,
    num_mixtures = 500,
    num_mixtures_sprop = 10,
    pool_size = 100,
    seed = 1234,
    prop = NULL,
    replace = FALSE){
  colnames(phenodata)[1:2] <- c("cellid", "celltype")
  if (length(unique(phenodata$cellid)) != nrow(phenodata)) stop("values of cellid in phenodata not unique.")
  if (length(intersect(colnames(ref), phenodata$cellid)) != length(union(colnames(ref), phenodata$cellid))) stop("column names of reference data do not match with cellid of the phenodata.")
  if (length(unique(phenodata$celltype)) < 2) stop("At least two unique cell types need to present in phenodata.")
  set.seed(seed)
  celltypes <- unique(phenodata$celltype)

  if (!is.null(prop)) {
    colnames(prop) <- c("ct", "expected")
    bulk <- lapply(1:num_mixtures, function(i) {
      chosen_cells <- sapply(which(prop$expected != 0), function(x) {
        ncells <- prop$expected[x] * pool_size
        chosen <- sample(
          phenodata$cellid[phenodata$celltype == prop$ct[x]],
          ncells, replace = replace
        )
        chosen
      }) %>% unlist()
      rowSums(ref[, colnames(ref) %in% chosen_cells]) %>% as.data.frame()
    })
    bulk <- do.call(cbind, bulk)
    colnames(bulk) <- paste("mix1_", 1:num_mixtures, sep = "")
    ## set CT as rownames
    rownames(prop) <- prop$ct
    prop <- prop[, rep(2, num_mixtures)]
    colnames(prop) <- colnames(bulk)
  } else {
    cell_distribution <- data.frame(table(phenodata$celltype), stringsAsFactors = FALSE)
    colnames(cell_distribution) <- c("ct", "maxn")

    tissues <- list()
    proportions <- list()
    rm(list = "prop")

    for (y in 1:(num_mixtures / num_mixtures_sprop)) {
      # Only allow feasible mixtures based on cell distribution
      while (!exists("prop")) {
        num_ct_mixture <- sample(x = 2:length(celltypes), 1)
        selectedct <- sample(celltypes, num_ct_mixture, replace = FALSE)

        prop <- runif(num_ct_mixture, 1, 99)
        prop <- round(prop / sum(prop), digits = log10(pool_size)) # sum to 1
        prop <- data.frame(ct = selectedct, expected = prop, stringsAsFactors = FALSE)

        missingct <- celltypes[!celltypes %in% selectedct]
        missingct <- data.frame(ct = missingct, expected = rep(0, length(missingct)), stringsAsFactors = FALSE)

        prop <- rbind.data.frame(prop, missingct)
        potential_mix <- merge(prop, cell_distribution)
        potential_mix$size <- potential_mix$expected * pool_size

        if (!all(potential_mix$maxn >= potential_mix$size) || sum(prop$expected) != 1) {
          rm(list = "prop")
        }
      }
      # Using info in prop to build bulk simultaneously
      # randomly sample 10 times to estimate bias and variance for each vector of prop
      bulk <- lapply(1:num_mixtures_sprop, function(i) {
        chosen_cells <- sapply(which(prop$expected != 0), function(x) {
          ncells <- prop$expected[x] * pool_size
          chosen <- sample(
            phenodata$cellid[phenodata$celltype == prop$ct[x]],
            ncells, replace = replace
          )
          chosen
        }) %>% unlist()
        rowSums(ref[, colnames(ref) %in% chosen_cells]) %>% as.data.frame()
      })
      bulk <- do.call(cbind, bulk)
      colnames(bulk) <- paste("mix", y, "_", 1:num_mixtures_sprop, sep = "")
      prop <- prop[, c("ct", "expected")]
      propn <- nrow(prop)
      prop <- prop[rep(seq_len(propn), times = num_mixtures_sprop), ]
      prop$mix <- paste("mix", y, "_", rep(1:num_mixtures_sprop, each = propn), sep = "")
      tissues[[y]] <- bulk
      proportions[[y]] <- prop

      rm(list = c("bulk", "prop", "missingct"))
    }

    prop <- do.call(rbind.data.frame, proportions)
    bulk <- do.call(cbind.data.frame, tissues)

    prop <- reshape2::dcast(prop, ct ~ mix,
      value.var = "expected",
      fun.aggregate = sum
    )
    prop <- as.data.frame(prop, row.names = 1)
    prop <- prop[, mixedsort(colnames(prop))]
  }
  return(list(bulk = bulk, prop = prop))
}

#' Statistical evaluations of predicted cell proportions
#'
#' Compute RMSE, bias & variance metrics for predicted cell proportions by comparing with expected cell proportions.
#'
#' @param prop_pred a matrix-like object of predicted cell proportion values with rows representing cell types, columns representing samples.
#' @param prop_sim a matrix-like object of simulated/expected cell proportion values with rows representing cell types, columns representing samples.
#'
#' @return a list of two objects:
#' \enumerate{
#'  \item a data.fame of summary metrics containing RMSE, bias & variance grouped by cell types and mixture ids (simulated samples with the same expected cell proportions).
#'  \item a data.frame of aggregated RMSE values across all cell types within each sample.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## generate artificial bulk samples
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
#'
#' ## compute metrics
#' metrics_res <- compute_metrics(decon_res[[1]], bulk_sim[[2]])
#' }



compute_metrics <- function(prop_pred, prop_sim) {
  if (length(union(colnames(prop_pred), colnames(prop_sim))) != length(intersect(colnames(prop_pred), colnames(prop_sim)))) stop("different column names between prop_pred and prop_sim")
  prop_sim <- prop_sim[, match(colnames(prop_pred), colnames(prop_sim))]
  prop_pred <- prop_pred[mixedsort(rownames(prop_pred)), ]
  prop_pred <- melt(as.matrix(prop_pred))
  colnames(prop_pred) <- c("celltype", "mixture", "observed_values")
  prop_sim <- prop_sim[mixedsort(rownames(prop_sim)), ]
  prop_sim <- melt(as.matrix(prop_sim))
  colnames(prop_sim) <- c("celltype", "mixture", "expected_values")
  prop_pred$expected_values <- prop_sim$expected_values
  prop_pred$mixture_samep <- gsub("_.*", "", prop_pred$mixture)
  prop_pred_sum <- prop_pred %>%
    group_by(celltype, mixture_samep) %>%
    summarise(
      RMSE = sqrt(mean((observed_values - expected_values)^2)),
      bias = mean(observed_values - expected_values),
      variance = mean((observed_values - mean(observed_values))^2)
    )
  prop_pred_sum1 <- prop_pred %>%
    group_by(mixture) %>%
    summarise(RMSE = sqrt(mean((observed_values - expected_values)^2)))
  return(list(prop_pred_sum, prop_pred_sum1))
}
