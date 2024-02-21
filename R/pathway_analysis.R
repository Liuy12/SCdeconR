
#' Methods to manipulate .gmt files
#'
#' Reformat, read & write .gmt file.
#'
#' @param gmtfile path to a gene set definition file in .gmt format.
#' @param gmt a gmt object returned by \code{read_gmt}.
#' @param outputfile full path including file name to export reformatted .gmt file.
#' @param replace a logical value indicating whether to replace the output file if it already exists. Default to FALSE.
#'
#' @return for read_gmt, returns a list object with length equal to the total number of gene sets within the .gmt file. Each list contains three elements:
#' "id", "name", "genes". read_gmt & write_gmt are reimplemented based on functions from \code{ActivePathways} package.
#'
#' @details \code{reformat_gmt} replaces blank spaces within the gene-set names to help string-matching methods in downstream plot functions \code{\link{gsea_sumplot}}
#' ,\code{\link{gsea_rwplot}}, \code{\link{gsea_heatmap}}.
#' @export
#'

reformat_gmt <- function(gmtfile, outputfile, replace = FALSE) {
  gmtinput <- read_gmt(gmtfile)
  for(i in 1:length(gmtinput)) gmtinput[[i]]$id <- gsub(" ", "_", gmtinput[[i]]$id)
  write_gmt(gmtinput, outputfile, replace = replace)
}

#' @rdname reformat_gmt
#' @export

read_gmt <- function(gmtfile) {
        gmt <- strsplit(readLines(gmtfile), "\t")
        names(gmt) <- sapply(gmt, function(i) i[1])
        gmt <- lapply(gmt, function(i) { list(id = i[1], name = i[2], genes = i[-c(1, 2)]) })
        gmt
}

#' @rdname reformat_gmt
#' @export

write_gmt <- function(gmt, outputfile, replace = FALSE) {
  if(file.exists(outputfile) && (!replace)) stop("outputfile exists!")
  sink(outputfile)
  for (term in gmt) {
      cat(term$id, term$name, paste(term$genes, collapse = "\t"), sep = "\t")
      cat("\n")
  }
  sink()
}


#' Prepare .rnk file for GSEA preranked analysis
#'
#' Prepare .rnk file for GSEA preranked analysis
#'
#' @param teststats a data.frame containing differential expression results with five columns: "Gene name", "log2 fold change", "log2 average expression",
#' "p value", "adjusted p value". The second element of the output from function \code{\link{run_de}}.
#' @param outputfile full path including file name to .rnk file.
#' @param replace a logical value indicating whether to replace the output file if it already exists. Default to FALSE.
#' @export


prepare_rnk <- function(teststats, outputfile, replace = FALSE){
  if(file.exists(outputfile) && (!replace)) stop("outputfile exists!")
  if (ncol(teststats) != 5) stop("inputs need to have five columns: genename, log2foldchange, log2avgexp, pval & padj")
  colnames(teststats) <- c("genename", "log2foldchange", "log2avgexp", "pval", "padj")
  teststats <- teststats[!duplicated(teststats$genename), ]
  write.table(teststats[,c("genename", "log2foldchange")], outputfile, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}


#' Generate a scatter plot comparing two gene set enrichment analysis results
#'
#' Generate a scatter plot of normalized enrichment scores comparing two results for gene set enrichment analysis, e.g. w/wo adjusting for cell proportion differences.
#'
#' @param gseares_path1 path to GSEA output.
#' @param gseares_path2 path to a second GSEA output.
#' @param result_names a vector of length 2 indicating the names of the two GSEA results. If NULL, names will be set to c("results1", "results2")
#' @param nes_cutoff normalized enrichment score cutoff to identify enriched gene-sets.
#' @param pval_cutoff p value cutoff to identify enriched gene-sets.
#' @param pvalflag a logical value indicating whether to use adjusted p value in selecting enriched gene-sets. Default to TRUE.
#' @param interactive a logical value indicating whether to generate an interactive plot. Default to FALSE.
#'
#' @export
#'
#' @details this function does not support output from GSEA \href{https://github.com/GSEA-MSigDB/GSEA_R}{R implementation}
#'


comparegsea_scatter <- function(
  gseares_path1,
  gseares_path2,
  result_names = NULL,
  nes_cutoff = 2,
  pval_cutoff = 0.1,
  pvalflag = TRUE,
  interactive = FALSE){
    if ((!is.null(result_names)) && length(result_names) != 2) stop("The length of result_names has to equal to 2") else if (is.null(result_names)) result_names <- c("result1", "result2")
    gseares1 <- do.call(rbind, lapply(list.files(gseares_path1, pattern = "gsea_report_for_.*.tsv", full.names = TRUE), function(i) read.delim(i, header = TRUE)))
    gseares2 <- do.call(rbind, lapply(list.files(gseares_path2, pattern = "gsea_report_for_.*.tsv", full.names = TRUE), function(i) read.delim(i, header = TRUE)))
    gseares3 <- inner_join(gseares1, gseares2, by = "NAME")
    if (pvalflag) {
        idx1 <- which(abs(gseares3$NES.x) >= log2(nes_cutoff) & gseares3$FDR.q.val.x <= pval_cutoff)
        idx2 <- which(abs(gseares3$NES.y) >= log2(nes_cutoff) & gseares3$FDR.q.val.y <= pval_cutoff)
    } else {
        idx1 <- which(abs(gseares3$NES.x) >= log2(nes_cutoff) & gseares3$NOM.p.val.x <= pval_cutoff)
        idx2 <- which(abs(gseares3$NES.y) >= log2(nes_cutoff) & gseares3$NOM.p.val.y <= pval_cutoff)
    }
   if (length(idx1) == 0 && length(idx2) == 0) stop("No enriched genesets identified, try loosing the thresholds")
    gseares3$category <- "Not significant"
    gseares3$category[intersect(idx1, idx2)] <- "Both"
    gseares3$category[setdiff(idx1, idx2)] <- paste0(result_names[1], " only")
    gseares3$category[setdiff(idx2, idx1)] <- paste0(result_names[2], " only")
    color_df <- data.frame(
        group = c("Both", paste0(result_names[1], " only"), paste0(result_names[2], " only"), "Not significant"),
        colors = c("red", "seagreen1", "royalblue", "grey")
    )
    gseares3$category <- factor(gseares3$category, levels = color_df$group[color_df$group %in% unique(gseares3$category)])
    lim_values <- range(gseares3[, c("NES.x", "NES.y")], na.rm = TRUE)
    options(warn = -1)
    gp <- ggplot() +
        geom_point(aes(x = NES.x, y = NES.y, color = category, text = paste0("Geneset: ", NAME)), size = 5, data = gseares3) +
        theme_classic() + labs(x = result_names[1], y = result_names[2]) +
        scale_color_manual(values = color_df$colors[color_df$group %in% unique(gseares3$category)]) +
        lims(x = lim_values, y = lim_values) +
        geom_hline(yintercept = c(log2(nes_cutoff), -log2(nes_cutoff)), linetype = 2) +
        geom_vline(xintercept = c(log2(nes_cutoff), -log2(nes_cutoff)), linetype = 2) +
        geom_abline(slope = 1, intercept = 0)
    options(warn = 0)
    if (interactive) {
        return(ggplotly(gp))
    } else {
        return(gp)
    }
}

#' Summary plot of gene set enrichment analysis
#'
#' Summary plot of selected up/down regulated gene-sets for gene set enrichment analysis.
#'
#' @param gseares_path path to GSEA output.
#' @param pos_sel a character vector of upregulated gene-set names.
#' @param neg_sel a character vector of downregulated gene-set names.
#' @param pvalflag a logical value indicating whether to use adjusted p value in the plot. Default to TRUE.
#' @param interactive a logical value indicating whether to generate an interactive plot. Default to FALSE.
#'
#' @export
#'
#' @details this function does not support output from GSEA \href{https://github.com/GSEA-MSigDB/GSEA_R}{R implementation}


gsea_sumplot <- function(
  gseares_path,
  pos_sel,
  neg_sel,
  pvalflag = TRUE,
  interactive = FALSE){
  gseares_up <- read.delim(list.files(gseares_path, pattern = "gsea_report_for_.*_pos_.*.tsv", full.names = TRUE), stringsAsFactors = FALSE)
  gseares_down <- read.delim(list.files(gseares_path, pattern = "gsea_report_for_.*_neg_.*.tsv", full.names = TRUE), stringsAsFactors = FALSE)
  idx_up <- na.omit(match(pos_sel, gseares_up$NAME))
  message(paste0("found ", length(idx_up), " # of up-regulated genesets"))
  idx_down <- na.omit(match(neg_sel, gseares_down$NAME))
  message(paste0("found ", length(idx_down), " # of down-regulated genesets"))
  if(length(idx_up) != length(pos_sel) || length(idx_down) != length(neg_sel))
    stop("Please check spelling of selected gene sets")
  gseares_up <- gseares_up[idx_up,]
  gseares_down <- gseares_down[idx_down,]
  tmp <- rbind(gseares_up, gseares_down)
  tmp$NAME <- gsub("_", " ", tmp$NAME)
  if(pvalflag) pval <- "FDR.q.val" else pval <- "NOM.p.val"
  options(warn = -1)
  gp <- ggplot(tmp, aes(NES, NAME)) +
    geom_point(aes(colour = get(pval), size = SIZE)) + labs(colour = pval) +
    scale_color_gradientn(colours=rainbow(4), limits=c(0, 1)) +
    geom_vline(xintercept = 0, size = 0.5, colour = "gray50") +
    theme(panel.background = element_rect(fill = "gray95", colour = "gray95"),
          panel.grid.major = element_line(size = 0.25, linetype = "solid", colour = "gray90"),
          panel.grid.minor = element_line(size = 0.25,linetype = "solid", colour = "gray90"),
          axis.title.y = element_blank()) +
    scale_y_discrete(limits = rev(tmp$NAME))
  options(warn = 0)
  if(interactive) return(ggplotly(gp)) else return(gp)
}


#' Heatmap to demonstrate enrichment of selected gene-sets
#'
#' Heatmap to demonstrate enrichment of selected gene-sets.
#'
#' @param normdata a matrix-like object of normalized & untransformed bulk RNA-seq data, with rows representing genes and columns representing samples.
#' The first element of the output from function \code{\link{run_de}}.
#' @param teststats a data.frame containing differential expression results with five columns: "Gene name", "log2 fold change", "log2 average expression",
#' "p value", "adjusted p value". The second element of the output from function \code{\link{run_de}}.
#' @param gmtfile path to gmt file used for GSEA analysis.
#' @param numgenes Number of genes to include in the heatmap. Will choose \code{numgenes} # of top up-regulated genes, as well as \code{numgenes} # of down-regulated genes
#' @param gsname_up a character value indicating selected up-regulated gene-set.
#' @param gsname_down a character value indicating selected down-regulated gene-set.
#' @param anncol a data.frame of sample meta information to include as column annotation bars. See option annCol from \code{\link[NMF]{aheatmap}} for more details.
#' @param color color used for heatmap. See option color option from \code{\link[NMF]{aheatmap}} for details.
#' @param anncolors optional data.frame to define colors for column annotations in \code{anncol}.
#' @param rankcol a logical value indicating whether to sort samples based on correlation between fold change & gene expression for better visualization. Default to TRUE.
#' @param zscore_range a vector of length two indicating the desired range of z-score transformed data. Default to c(-3, 3).
#' @export
#'
#' @details this function does not support output from GSEA \href{https://github.com/GSEA-MSigDB/GSEA_R}{R implementation}


gsea_heatmap <- function(
  normdata,
  teststats,
  gmtfile,
  numgenes,
  gsname_up,
  gsname_down,
  anncol,
  color,
  anncolors = NULL,
  rankcol = TRUE,
  zscore_range = c(-3, 3)){
  if (ncol(teststats) != 5) stop("inputs need to have five columns: genename, log2foldchange, log2avgexp, pval & padj")
  colnames(teststats) <- c("genename", "log2foldchange", "log2avgexp", "pval", "padj")
  if (length(intersect(rownames(normdata), teststats$genename)) != length(union(rownames(normdata), teststats$genename))) stop("row names of normdata not consistent with genename of teststats")
  teststats <- teststats[match(rownames(normdata), teststats$genename),]
  gmt <- read_gmt(gmtfile)
  or <- order(teststats$log2foldchange, decreasing = TRUE)
  normdata_sel <- normdata[c(or[1:numgenes], tail(or, numgenes)),]
  idx <- which(toupper(names(gmt)) == toupper(gsname_up))
  if(!length(idx)) stop(paste0(gsname_up, " not found in provided gmtfile"))
  ann_gs1 <- sapply(rownames(normdata_sel), function(i) {
    idx1 <- grep(paste0("^", i, "$"), gmt[[idx]]$genes)
    if(length(idx1)) "In GS" else "Not in GS"
  })
  idx <- which(toupper(names(gmt)) == toupper(gsname_down))
  if(!length(idx)) stop(paste0(gsname_down, " not found in provided gmtfile"))
  ann_gs2 <- sapply(rownames(normdata_sel), function(i) {
    idx1 <- grep(paste0("^", i, "$"), gmt[[idx]]$genes)
    if(length(idx1)) "In GS" else "Not in GS"
  })
  rowann <- data.frame(gsname_up = ann_gs1,
                       gsname_down = ann_gs2)
  if(is.null(anncolors)) anncolors <- list(gsname_up = c('black', 'white'), gsname_down = c('green','white'))
  tmp <- log2(normdata_sel + 0.1)
  if(rankcol) {
    colrk <- order(sapply(1:(ncol(tmp)), function(i) cor(tmp[,i], teststats$log2foldchange[c(or[1:numgenes], tail(or, numgenes))])), decreasing = FALSE)
  } else colrk <- 1:ncol(tmp)
  tmp <- t(scale(t(tmp)))
  tmp[tmp < zscore_range[1]] <- zscore_range[1]
  tmp[tmp > zscore_range[2]] <- zscore_range[2]
  NMF::aheatmap(tmp, annRow = rowann, annCol = anncol, scale = "none", labCol = NA, labRow = NA, annColors = anncolors, color = color, Rowv = NA, Colv = NA)
  }


#' GSEA random-walk plot
#'
#' Generate high-quality GSEA random-walk figures.
#'
#' @param gseares_path path to GSEA output.
#' @param gsname a character value indicating the gene-set name to generate random-walk plot.
#' @param class_name a character value indicating the class of the gene-set, e.g. "GO".
#' @param metric_range optional range of the ranking metric.
#'
#' @export
#'
#' @details this function does not support output from GSEA \href{https://github.com/GSEA-MSigDB/GSEA_R}{R implementation}. Scripts initially implemented by Thomas Kuilman.
#'
#'


gsea_rwplot <- function(gseares_path, gsname, class_name, metric_range = NULL) {
  if (!dir.exists(gseares_path)) {
    stop(paste0(gsearespath, " not found"))
  }
  if (missing(gsname)) {
    stop("gsname argument is required")
  }

  ## Load .rnk data
  path_rnk <- list.files(path = file.path(gseares_path, "edb"), pattern = ".rnk$", full.names = TRUE)
  gsea_rnk <- read.delim(file = path_rnk, header = FALSE)
  colnames(gsea_rnk) <- c("hgnc.symbol", "metric")
  if (is.null(metric_range)) metric_range <- range(gsea_rnk$metric)

  ## Load .edb data
  path_edb <- list.files(path = file.path(gseares_path, "edb"), pattern = ".edb$", full.names = TRUE)
  gsea_edb <- read.delim(file = path_edb, header = FALSE, stringsAsFactors = FALSE)
  gsea_edb <- unlist(gsea_edb)
  gsea_metric <- gsea_edb[grep("METRIC=", gsea_edb)]
  gsea_metric <- unlist(strsplit(gsea_metric, " "))
  gsea_metric <- gsea_metric[grep("METRIC=", gsea_metric)]
  gsea_metric <- gsub("METRIC=", "", gsea_metric)
  gsea_edb <- gsea_edb[grep("<DTG", gsea_edb)]

  # Select the right gene set
  if (length(gsea_edb) == 0) stop(paste0("The gene set name:", gsname, " was not found, please provide a correct name"))
  if (length(grep(paste0(gsub(".\\$(.*$)", "\\1", gsname), " "), gsea_edb)) > 1) warning(paste("More than 1 gene set matched the gene.set, argument; the first match is plotted"))
  gsea_edb <- gsea_edb[grep(gsname, gsea_edb)[1]]
  # Get template name
  gsea_edb <- gsub(".*TEMPLATE=(.*)", "\\1", gsea_edb)
  gsea_edb <- unlist(strsplit(gsea_edb, " "))
  gsea_template <- gsea_edb[1]

  # Get gene set name
  gsea_gene_set <- gsea_edb[2]
  gsea_gene_set <- gsub("GENESET=gene_sets.gmt#", "", gsea_gene_set)

  # Get enrichment score
  gsea_enrichment_score <- gsea_edb[3]
  gsea_enrichment_score <- gsub("ES=", "", gsea_enrichment_score)

  # Get normalized enrichment score
  gsea_normalized_enrichment_score <- gsea_edb[4]
  gsea_normalized_enrichment_score <- gsub("NES=", "", gsea_normalized_enrichment_score)

  # Get nominal p-value
  gsea_pvalue <- gsea_edb[5]
  gsea_pvalue <- gsub("NP=", "", gsea_pvalue)
  gsea_pvalue <- as.numeric(gsea_pvalue)

  # Get FDR
  gsea_fdr <- gsea_edb[6]
  gsea_fdr <- gsub("FDR=", "", gsea_fdr)
  gsea_fdr <- as.numeric(gsea_fdr)

  # Get hit indices
  gsea_edb <- gsea_edb[grep("HIT_INDICES=", gsea_edb):length(gsea_edb)]
  gsea_hit_indices <- gsea_edb[seq_len(grep("ES_PROFILE=", gsea_edb) - 1)]
  gsea_hit_indices <- gsub("HIT_INDICES=", "", gsea_hit_indices)
  gsea_hit_indices <- as.integer(gsea_hit_indices)

  # Get ES profile
  gsea_edb <- gsea_edb[grep("ES_PROFILE=", gsea_edb):length(gsea_edb)]
  gsea_es_profile <- gsea_edb[seq_len(grep("RANK_AT_ES=", gsea_edb) - 1)]
  gsea_es_profile <- gsub("ES_PROFILE=", "", gsea_es_profile)
  gsea_es_profile <- as.numeric(gsea_es_profile)

  ## Create GSEA plot
  # Save default for resetting
  def_par <- par(no.readonly = TRUE)

  # Create a division of the device
  gsea_layout <- layout(matrix(c(1, 2, 3, 4)), heights = c(1.7, 0.5, 0.2, 2))

  # Create plots
  par(mar = c(0, 5, 2, 2))
  plot(c(1, gsea_hit_indices, length(gsea_rnk$metric)),
       c(0, gsea_es_profile, 0), type = "l", col = "red", lwd = 1.5, xaxt = "n",
       xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
       main = list(gsea_gene_set, font = 1, cex = 1),
       panel.first = {
         abline(h = seq(round(min(gsea_es_profile), digits = 1),
                        max(gsea_es_profile), 0.1),
                col = "gray95", lty = 2)
         abline(h = 0, col = "gray50", lty = 2)
       })
  plot_coordinates <- par("usr")
  if(gsea_enrichment_score < 0) {
    text(length(gsea_rnk$metric) * 0.01, plot_coordinates[3] * 0.98,
         paste("Nominal p-value:", gsea_pvalue, "\nFDR:", gsea_fdr, "\nES:",
               gsea_enrichment_score, "\nNormalized ES:",
               gsea_normalized_enrichment_score), adj = c(0, 0))
  } else {
    text(length(gsea_rnk$metric) * 0.99, plot_coordinates[4] - ((plot_coordinates[4] - plot_coordinates[3]) * 0.03),
         paste("Nominal p-value:", gsea_pvalue, "\nFDR:", gsea_fdr, "\nES:",
               gsea_enrichment_score, "\nNormalized ES:",
               gsea_normalized_enrichment_score, "\n"), adj = c(1, 1))
  }

  par(mar = c(0, 5, 0, 2))
  plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
       ylab = "", xlim = c(1, length(gsea_rnk$metric)))
  abline(v = gsea_hit_indices, lwd = 0.75)

  par(mar = c(0, 5, 0, 2))
  rank_colors <- gsea_rnk$metric - metric_range[1]
  rank_colors <- rank_colors / (metric_range[2] - metric_range[1])
  rank_colors <- ceiling(rank_colors * 255 + 1)
  tryCatch({
    rank_colors <- colorRampPalette(c("blue", "white", "red"))(256)[rank_colors]
  }, error = function(e) {
    stop("Please use the metric_range argument to provide a metric range that",
         "includes all metric values")
  })
  # Use rle to prevent too many objects
  rank_colors <- rle(rank_colors)
  barplot(matrix(rank_colors$lengths), col = rank_colors$values, border = NA, horiz = TRUE, xaxt = "n", xlim = c(1, length(gsea_rnk$metric)))
  box()
  text(length(gsea_rnk$metric) / 2, 0.7,
       labels = ifelse(!missing(class_name), class_name, gsea_template))
  text(length(gsea_rnk$metric) * 0.01, 0.7, "Positive", adj = c(0, NA))
  text(length(gsea_rnk$metric) * 0.99, 0.7, "Negative", adj = c(1, NA))

  par(mar = c(5, 5, 0, 2))
  rank_metric <- rle(round(gsea_rnk$metric, digits = 2))
  plot(gsea_rnk$metric, type = "n", xaxs = "i",
       xlab = "Rank in ordered gene list", xlim = c(0, length(gsea_rnk$metric)),
       ylim = metric_range, yaxs = "i",
       ylab = if(gsea_metric == "None") {"Ranking metric"} else {gsea_metric},
       panel.first = abline(h = seq(metric_range[1] / 2,
                                    metric_range[2] - metric_range[1] / 4,
                                    metric_range[2] / 2), col = "gray95", lty = 2))

  barplot(rank_metric$values, col = "lightgrey", lwd = 0.1, xaxs = "i",
          xlab = "Rank in ordered gene list", xlim = c(0, length(gsea_rnk$metric)),
          ylim = c(-1, 1), yaxs = "i", width = rank_metric$lengths, border = NA,
          ylab = ifelse(gsea_metric == "None", "Ranking metric", gsea_metric), space = 0, add = TRUE)
  box()
  # Reset to default
  par(def_par)
}
