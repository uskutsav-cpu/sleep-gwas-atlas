# Independent base-R check of the exact published paired-model contrast.
# This checks interpretation of the design, not a DESeq2 coefficient refit.
root <- '/Users/swethasunilkumar/Documents/Codex/2026-10-27/loo'
out <- file.path(root, 'outputs/sleep_frailty_ml_direction')
path <- file.path(out, 'sources/MPRAmodel_e853506.R')
code <- readLines(path, warn=FALSE)
stopifnot(any(grepl('levels = c("ref","alt")',code,fixed=TRUE)))
stopifnot(any(grepl('design(dds) <- ~material + material:sample.n + material:allele',code,fixed=TRUE)))
stopifnot(any(grepl('contrast=list("materialRNA.allelealt","materialDNA.allelealt")',code,fixed=TRUE)))
sam <- data.frame(material=factor(rep(c('DNA','RNA'),each=10)),
                  allele=factor(rep(c('ref','alt'),10),levels=c('ref','alt')),
                  sample.n=factor(rep(rep(LETTERS[1:5],each=2),2)))
design <- model.matrix(~material + material:sample.n + material:allele,sam)
stopifnot(qr(design)$rank==ncol(design))
stopifnot(all(c('materialRNA:allelealt','materialDNA:allelealt') %in% colnames(design)))
errors <- numeric()
for (dna in c(-.8,0,.35)) for (rna in c(-.4,0,.95)) {
  # Deterministic log2 means with replicate shifts exactly represented by design.
  y <- 12 + ifelse(sam$material=='RNA',-1,0) +
       (as.numeric(sam$sample.n)-1)*.2 +
       ifelse(sam$allele=='alt',ifelse(sam$material=='RNA',rna,dna),0)
  coefficients <- qr.solve(design,y)
  observed <- coefficients['materialRNA:allelealt']-coefficients['materialDNA:allelealt']
  errors <- c(errors,abs(observed-(rna-dna)))
}
stopifnot(max(errors)<1e-12)
write.table(data.frame(passed=TRUE,deterministic_design_cases=length(errors),
                      max_contrast_error=max(errors),
                      contrast='log2(ALT RNA / REF RNA) - log2(ALT DNA / REF DNA)',
                      original_DESeq2_coefficients_refitted=FALSE),
            file.path(out,'paired_contrast_independent_R_receipt.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat('PAIRED CONTRAST ALGEBRA INDEPENDENT R CHECK PASS\n')
