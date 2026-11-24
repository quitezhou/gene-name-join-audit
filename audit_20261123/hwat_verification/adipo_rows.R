# Read just the row names of the published adipocyte count slot.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SeuratObject))
o <- readRDS(args[1])
a <- o@assays$RNA@counts
k <- c("AHRR","ARMCX5-GPRASP2","ATXN7","PINX1","POLR2J3")
tab <- do.call(rbind, lapply(k,function(label){
    variants <- unique(c(label,paste0(label,".1"),paste0(label,"-1")))
    do.call(rbind,lapply(variants,function(v){
        i <- match(v,rownames(a))
        data.frame(label=label,row_name=v,one_based_row=i,atlas_nrows=nrow(a),atlas_ncols=ncol(a))
    }))
}))
write.table(tab,args[2],sep="\t",row.names=FALSE,quote=FALSE)
# Keep all cells as stored and export only five duplicate families for a direct
# comparison of base and suffixed rows against the deposit's stable-ID matrices.
sel <- which(rownames(a) %in% c(k,paste0(k,".1"),paste0(k,"-1")))
root <- dirname(args[2])
cache <- file.path(root,"source")
dir.create(cache,showWarnings=FALSE)
writeLines(rownames(a)[sel],file.path(cache,"adipo_names.txt"))
write.table(data.frame(cell=colnames(a),cohort=o@meta.data$orig.ident),file.path(cache,"adipo_cells.tsv"),row.names=FALSE,sep="\t",quote=FALSE)
writeMM(a[sel,,drop=FALSE],file.path(cache,"adipo_selected.mtx"))
