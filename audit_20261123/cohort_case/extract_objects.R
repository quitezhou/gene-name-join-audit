# Export the deposited RNA counts by row name; preserve stored numeric values.
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly=TRUE)
inv <- read.delim("deposit_inventory.tsv",check.names=FALSE)
wanted <- unique(c(readLines("candidate_symbols.txt"),readLines("control_symbols.txt")))
cache_dir <- Sys.getenv("HWAT_CACHE",unset=Sys.getenv("TEMP",unset=tempdir()))
for (assay in args) {
    filename <- paste0(assay,".rds")
    rec <- inv[basename(inv$path)==filename,]
    stopifnot(nrow(rec)==1L)
    f <- file.path(cache_dir,paste0("hwat_",rec$file_id,"_",filename))
    x <- readRDS(f)
    a <- x@assays$RNA@counts
    meta <- x@meta.data
    stopifnot(identical(colnames(a),rownames(meta)))
    keep <- which(grepl("^Massier",meta$orig.ident))
    genes <- intersect(wanted,rownames(a))
    prefix <- file.path(cache_dir,paste0("hwat_extract_",assay))
    writeLines(genes,paste0(prefix,"_features.txt"))
    write.table(data.frame(cell=colnames(a)[keep],cohort=as.character(meta$orig.ident[keep])),
                paste0(prefix,"_cells.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
    writeMM(a[genes,keep,drop=FALSE],paste0(prefix,"_counts.mtx"))
    cat(assay,dim(a),length(keep),length(genes),"\n")
}
