# Export the released data column for each model gene plus alternate labels.
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))
a <- commandArgs(TRUE)
x <- readRDS(a[1])
source <- x@assays$RNA@counts
requested <- readLines(a[2], encoding="UTF-8")
barcodes <- read.delim(a[3], sep="\t", stringsAsFactors=FALSE)
stopifnot(identical(colnames(source), rownames(x@meta.data)))
raw_bc <- sub("_.*$", "", colnames(source))
ix <- match(raw_bc, barcodes$barcode)
cohort <- as.character(x@meta.data$orig.ident)
valid <- (cohort == "Massier et al. #1" & barcodes$source_number[ix] == 1L) |
         (cohort == "Massier et al. #4" & barcodes$source_number[ix] == 4L)
idx <- which(!is.na(ix) & valid)
keys <- intersect(rownames(source), requested)
stopifnot(length(idx)>0, !anyDuplicated(raw_bc[idx]), !anyDuplicated(keys))
selected <- source[keys, idx, drop=FALSE]
writeMM(selected, paste0(a[4], ".mtx"))
writeLines(keys, paste0(a[4], "_keys.txt"))
df <- data.frame(barcode=raw_bc[idx], source=barcodes$source_number[ix[idx]],
                 total=Matrix::colSums(source[,idx,drop=FALSE]), cluster=as.character(Idents(x)[idx]))
write.table(df,paste0(a[4], "_obs.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
