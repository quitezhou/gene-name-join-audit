# Run after check_source_axes.py. Cache path and output directory passed as arguments.
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(limma))
args <- commandArgs(trailingOnly=TRUE)
cache <- args[1]
out <- args[2]
load(file.path(cache,"sce_layer.Rdata"))
load(file.path(cache,"eb_contrasts.Rdata"))
x <- read.delim(file.path(out,"family_layer_counts.tsv"), stringsAsFactors=FALSE)
target <- unique(x$anchor_id)
protection <- x[x$arm == "protected",]
observed <- assay(sce_layer,"counts")
loc <- cbind(match(protection$anchor_id,rownames(observed)), match(protection$section_layer,colnames(observed)))
stopifnot(nrow(loc)==length(target)*76L, !anyNA(loc), all(observed[loc]==protection$umi))
source_check <- data.frame(section_layer=protection$section_layer,id=protection$anchor_id,
                           H5_sum=protection$umi,archive_count=observed[loc])
write.table(source_check,file.path(out,"archive_count_check.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
expr <- assay(sce_layer,"logcounts")
stopifnot(identical(rownames(expr),rownames(eb_contrasts)), nrow(expr)==22331L)
models <- list()
for (arm in c("protected","filtered")) {
    working <- expr
    if (arm=="filtered") {
        rows <- x[x$arm=="filtered",]
        for (i in seq_len(nrow(rows))) {
            key <- rows$source_id[i]
            cc <- rows$section_layer[i]
            stopifnot(rows$umi[i] == (if (nzchar(key)) observed[key,cc] else 0))
            working[rows$anchor_id[i],cc] <- if (nzchar(key)) expr[key,cc] else 0
        }
    }
    fit <- lmFit(working,design=eb_contrasts$design,block=eb_contrasts$block,correlation=eb_contrasts$correlation)
    models[[arm]] <- eBayes(contrasts.fit(fit,eb_contrasts$contrasts))
}
cn <- "WM-Layer3"
# Use the authors' contrast orientation and retain the full testing family for BH.
fit_table <- function(z,arm) {
    pp <- z$p.value[,cn]
    data.frame(arm=arm, id=target, contrast=cn, coefficient=z$coefficients[target,cn],
               p=pp[target], p_BH=p.adjust(pp,"BH")[target], row.names=NULL)
}
tbl <- do.call(rbind, list(fit_table(eb_contrasts,"published"),fit_table(models$protected,"protected_replay"),
                          fit_table(models$filtered,"filter_then_name")))
write.table(tbl,file.path(out,"layer_model.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
evidence <- data.frame(id=rownames(expr),paper_beta=eb_contrasts$coefficients[,cn],
                       protected_beta=models$protected$coefficients[,cn],
                       counterfactual_beta=models$filtered$coefficients[,cn],
                       paper_p=eb_contrasts$p.value[,cn],protected_p=models$protected$p.value[,cn],
                       counterfactual_p=models$filtered$p.value[,cn],row.names=NULL)
handle <- gzfile(file.path(out,"full_layer_fit.tsv.gz"),"wt")
write.table(evidence,handle,sep="\t",row.names=FALSE,quote=FALSE)
close(handle)
resid <- data.frame(genes=nrow(expr),groups=ncol(expr),target_IDs=length(target),
                    checked_source_counts=nrow(loc),max_error_counts=max(abs(observed[loc]-protection$umi)),
                    contrast=cn,consensus_correlation=eb_contrasts$correlation,
                    max_coef_residual=max(abs(models$protected$coefficients[,cn]-eb_contrasts$coefficients[,cn])),
                    max_raw_p_residual=max(abs(models$protected$p.value[,cn]-eb_contrasts$p.value[,cn])),
                    protected_calls=sum(p.adjust(models$protected$p.value[,cn],"BH")<.05),
                    paper_calls=sum(p.adjust(eb_contrasts$p.value[,cn],"BH")<.05),
                    filtered_calls=sum(p.adjust(models$filtered$p.value[,cn],"BH")<.05),
                    changed_calls=sum((p.adjust(models$protected$p.value[,cn],"BH")<.05) !=
                                      (p.adjust(models$filtered$p.value[,cn],"BH")<.05)),
                    R_version=as.character(getRversion()),limma=as.character(packageVersion("limma")))
write.table(resid,file.path(out,"layer_model_check.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
print(resid);print(tbl)
