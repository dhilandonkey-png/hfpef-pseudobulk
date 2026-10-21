#!/usr/bin/env Rscript
# Pseudobulk differential expression with limma-voom, one model per cell type,
# as in Hahn, Chaffin et al. (Circ Res 2026): ~ group + sex + pool, BH FDR 0.05.
#
# Usage:
#   Rscript limma_voom.R <pseudobulk_out_dir> <group_col> <case_level> <control_level> [covariate ...]
#   Rscript limma_voom.R results/ disease HFpEF control sex pool
#
# Needs: install.packages("BiocManager"); BiocManager::install(c("limma", "edgeR"))

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: limma_voom.R <dir> <group_col> <case> <control> [covariates...]")
in_dir <- args[1]
group_col <- args[2]
case_level <- args[3]
control_level <- args[4]
covariates <- if (length(args) > 4) args[5:length(args)] else character(0)

run_de <- function(counts_file, meta_file, label) {
  counts <- as.matrix(read.csv(counts_file, row.names = 1, check.names = FALSE))
  meta <- read.csv(meta_file, row.names = 1, check.names = FALSE)
  meta <- meta[colnames(counts), , drop = FALSE]
  meta$group <- relevel(factor(meta[[group_col]]), ref = control_level)

  # drop covariates with a single level in this cell type (e.g. one pool)
  covs <- covariates[sapply(covariates, function(c) length(unique(meta[[c]])) > 1)]
  for (c in covs) if (!is.numeric(meta[[c]])) meta[[c]] <- factor(meta[[c]])
  design <- model.matrix(as.formula(paste("~", paste(c("group", covs), collapse = " + "))), data = meta)

  dge <- DGEList(counts)
  keep <- filterByExpr(dge, design = design)
  dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE], method = "TMM")

  v <- voom(dge, design)
  fit <- eBayes(lmFit(v, design))
  coef <- paste0("group", case_level)
  res <- topTable(fit, coef = coef, number = Inf, sort.by = "P")  # adj.P.Val is BH
  res$gene <- rownames(res)
  res$celltype <- label
  cat(sprintf("%-25s samples=%3d genes=%6d DE(FDR<0.05)=%d\n",
              label, ncol(dge), nrow(dge), sum(res$adj.P.Val < 0.05)))
  res
}

out_dir <- file.path(in_dir, "de_limma_voom")
dir.create(out_dir, showWarnings = FALSE)

ct_dir <- file.path(in_dir, "by_celltype")
for (f in list.files(ct_dir, pattern = "_counts\\.csv\\.gz$", full.names = TRUE)) {
  label <- sub("_counts\\.csv\\.gz$", "", basename(f))
  res <- run_de(f, file.path(ct_dir, paste0(label, "_metadata.csv")), label)
  write.csv(res, file.path(out_dir, paste0(label, "_de.csv")), row.names = FALSE)
}

res <- run_de(file.path(in_dir, "all_celltypes_counts.csv.gz"),
              file.path(in_dir, "all_celltypes_metadata.csv"), "all_celltypes")
write.csv(res, file.path(out_dir, "all_celltypes_de.csv"), row.names = FALSE)
