suppressPackageStartupMessages(library(SeuratObject))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2L)
obj <- readRDS(args[1])
counts <- obj@assays$RNA@counts
genes <- strsplit(paste(
  "ACTB ACTG1 GAPDH UBC B2M HLA-A HLA-B HLA-C HLA-E CD3D CD3E CD3G CD247 LCK TRAC CD2 CD7 CD52 CD27 CD74",
  "PTPRC CCR7 IL7R NKG7 KLRD1 GZMA GNLY GZMK GZMB CXCR4 LTB IL32 MALAT1 RPLP0 RPLP1 RPLP2 RPS3 RPS6 RPS8 RPS12",
  "RPS18 RPL3 RPL4 RPL5 RPL7 RPL11 RPL13A RPL18 RPL23 RPL27 EEF1A1 EEF1G TPT1 TMSB4X TMSB10 PPIA BTF3 CST3 LST1 TYROBP"
), " ")[[1]]
stopifnot(length(genes) == 60L)
keep <- intersect(genes,rownames(counts))
out <- data.frame(cell=colnames(counts),check.names=FALSE)
for (g in keep) out[[g]] <- as.integer(counts[g,] > 0)
write.table(out,gzfile(args[2]),sep="\t",row.names=FALSE,quote=FALSE)
