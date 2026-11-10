# Invoked only by the frozen final evaluator after it validates the lock.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L, file.exists(args[3]))
x <- readRDS(args[1])
required <- c('name','keep','label_class','padj','pval','sign','zscore_diff')
stopifnot(all(required %in% names(x)),!anyDuplicated(x$name))
optional <- c('tile_id','label','type','class','BC_count','zscore_neg',
              'zscore_ref','log2.value','var_id','chr','pos','start','end',
              'ref','alt','ref_dbsnp','alt_dbsnp','rsid','spaced_mut',
              'mutation_reference','mutation_valentina','mutation_panel',
              'mutation_gwas','snp_perimeter')
columns <- c(required,intersect(optional,names(x)))
# Select atomic columns only; retain all original rows and missing values.
columns <- columns[vapply(x[columns],function(z) is.atomic(z)&&is.null(dim(z)),logical(1))]
stopifnot(all(required %in% columns))
options(digits=17)
write.csv(as.data.frame(x[columns]),args[2],row.names=FALSE,na='NA')
cat('Frozen-source outcome export:',nrow(x),'rows;',length(columns),'columns\n')
