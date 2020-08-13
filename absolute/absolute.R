# https://software.broadinstitute.org/cancer/cga/absolute_run
# https://github.com/broadinstitute/PhylogicNDT/issues/4
library(ABSOLUTE)

seg.dat.fn <- "cldoxdmso32.cr.igv.seg"
sigma.p <- 0
max.sigma.h <- 0.02
min.ploidy <- 0.95
max.ploidy <- 10
primary.disease <- "brc"
platform <- "Illumina_WES"
sample.name <- "cldoxdmso32"
results.dir <- "."
max.as.seg.count <- 1500
max.non.clonal <- 0
max.neg.genome <- 0
copy_num_type <- "allelic"

RunAbsolute(seg.dat.fn,
	    sigma.p,
	    max.sigma.h,
	    min.ploidy,
	    max.ploidy,
	    primary.disease,
	    platform,
	    sample.name,
	    results.dir,
	    max.as.seg.count,
	    max.non.clonal,
	    max.neg.genome,
	    copy_num_type,
	    verbose=TRUE)
