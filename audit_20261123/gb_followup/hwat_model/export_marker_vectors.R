suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
obj <- readRDS(args[1])
out <- args[2]
a <- obj@assays$RNA
markers <- c("ATXN7","POLR2J3","PINX1","AHRR","ARMCX5-GPRASP2")
features <- intersect(as.vector(outer(markers,c("",".1"),paste0)),rownames(a))
cols <- data.frame(cell=colnames(a), cluster=as.character(Idents(obj)),source=as.character(obj@meta.data$orig.ident))
for (gene in features) {
    cols[[paste0("counts_",gene)]] <- as.numeric(a@counts[gene,])
    cols[[paste0("data_",gene)]] <- as.numeric(a@data[gene,])
}
write.table(cols,gzfile(out),row.names=FALSE,sep="\t",quote=FALSE)
cat("object",basename(args[1]),"active assay",obj@active.assay,"assays",names(obj@assays),"cells",nrow(cols),"cluster levels",paste(levels(Idents(obj)),collapse=","),"\n")
