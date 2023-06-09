% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_bulk.R
\name{compute_metrics}
\alias{compute_metrics}
\title{Statistical evaluations of predicted cell proportions}
\usage{
compute_metrics(prop_pred, prop_sim)
}
\arguments{
\item{prop_pred}{a matrix-like object of predicted cell proportion values with rows representing cell types, columns representing samples.}

\item{prop_sim}{a matrix-like object of simulated/expected cell proportion values with rows representing cell types, columns representing samples.}
}
\value{
a list of two objects:
\enumerate{
 \item a data.fame of summary metrics containing RMSE, bias & variance grouped by cell types and mixture ids (simulated samples with the same expected cell proportions).
 \item a data.frame of aggregated RMSE values across all cell types within each sample.
}
}
\description{
Compute RMSE, bias & variance metrics for predicted cell proportions by comparing with expected cell proportions.
}
\examples{
\dontrun{
## generate artificial bulk samples
data(refdata)
phenodata <- data.frame(cellid = colnames(refdata),
                        celltypes = refdata$celltype,
                        subjectid = refdata$subjectid)
bulk_sim <- bulk_generator(ref = GetAssayData(refdata, slot = "data", assay = "SCT"), 
                           phenodata = phenodata, 
                           num_mixtures = 20, 
                           num_mixtures_sprop = 1)

## perform deconvolution based on "OLS" algorithm
decon_res <- scdecon(bulk = bulk_sim[[1]],
                     ref = GetAssayData(refdata, slot = "data", assay = "SCT"),
                     phenodata = phenodata,
                     filter_ref = TRUE,
                     decon_method = "OLS",
                     norm_method = "none",
                     trans_method = "none",
                     marker_strategy = "all")

## compute metrics 
metrics_res <- compute_metrics(decon_res[[1]], bulk_sim[[2]])
}
}
