library(digitalDLSorteR)
library(SingleCellExperiment)

sce <- SingleCellExperiment(assays = ref, 
                            rowData = data.frame(Gene_ID = rownames(ref)),
                            colData = data.frame(Cell_ID = colnames(ref),
                                                 Cell_Type = phenodata$celltype))

bulk_data <- SummarizedExperiment(assays = bulkdata,
                                  rowData = data.frame(Gene_ID = rownames(bulkdata)),
                                  colData = data.frame(Sample_ID = colnames(bulkdata)))

ddls_obj <- createDDLSobject(sc.data = sce,sc.cell.ID.column = "Cell_ID", sc.gene.ID.column = "Gene_ID", sc.cell.type.column = "Cell_Type", bulk.data = bulk_data, bulk.sample.ID.column = "Sample_ID", bulk.gene.ID.column = "Gene_ID")

ddls_obj <- estimateZinbwaveParams(
  object = ddls_obj,
  cell.ID.column = "Cell_ID",
  gene.ID.column = "Gene_ID",
  cell.type.column = "Cell_Type",
  cell.cov.columns = NULL,
  gene.cov.columns = NULL,
  threads = max(1, parallel::detectCores() -1),
  subset.cells = 1000,
  verbose = TRUE
)

#### simulate new sc data
ddls_obj <- simSCProfiles(
  object = ddls_obj,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  n.cells = 100,
  suffix.names = "_Simul",
  verbose = TRUE
)

## for reproducibility
set.seed(123)

## prior knowledge for prob.design argument
## need to modify based on your experiment
probMatrix <- data.frame(
  Cell_Type = unique(ddls_obj@single.cell.real$Cell_Type),
  from = c(15,1,15,1,1,15,1,15),
  to = c(70,15,30,15,15,30,15,70)
)

### generate proportion matrix
ddls_obj <- generateBulkCellMatrix(
  object = ddls_obj,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  prob.design = probMatrix,
  num.bulk.samples = 1000,
  n.cells = 100,
  verbose = TRUE
)

ddls_obj <- simBulkProfiles(
  object = ddls_obj, type.data = "both", pseudobulk.function = "AddRawCount"
)

reticulate::use_python("/path/to/conda/env/python/")
ddls_obj <- trainDDLSModel(object = ddls_obj, scaling = "standardize")
ddls_obj <- calculateEvalMetrics(object = ddls_obj)

### test on bulk data

ddls_obj <- deconvDDLSObj(
  object = ddls_obj, 
  name.data = "Bulk.DT",
  normalize = TRUE,
  scaling = "standardize",
  verbose = FALSE
)
