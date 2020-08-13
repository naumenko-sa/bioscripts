# https://software.broadinstitute.org/cancer/cga/absolute_run
# https://github.com/broadinstitute/PhylogicNDT/issues/4
library(ABSOLUTE)

DoAbsolute <- function(scan, sif) {
    registerDoSEQ()
    library(ABSOLUTE)
    plate.name <- "DRAWS"
    genome <- "hg38"
    platform <- "Illumina_WES"

    primary.disease <- sif[scan, "PRIMARY_DISEASE"]
    sample.name <- sif[scan, "SAMPLE_NAME"]
    sigma.p <- 0
    max.sigma.h <- 0.02
    min.ploidy <- 0.95
    max.ploidy <- 10
    max.as.seg.count <- 1500
    max.non.clonal <- 0
    max.neg.genome <- 0
    copy_num_type <- "allelic"
    seg.dat.fn <- file.path("output", scan, "hapseg",
                            paste(plate.name, "_", scan, "_segdat.RData", sep=""))
                            results.dir <- file.path(".", "output", scan, "absolute")
                            print(paste("Starting scan", scan, "at", results.dir))
                            log.dir <- file.path(".", "output", "abs_logs")
                            if (!file.exists(log.dir)) {
                                   dir.create(log.dir, recursive=TRUE)
                            }
                            if (!file.exists(results.dir)) {
                                  dir.create(results.dir, recursive=TRUE)
                            }
    sink(file=file.path(log.dir, paste(scan, ".abs.out.txt", sep="")))
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
        sink()
}

arrays.txt <- "./paper_example/mix250K_arrays.txt"
sif.txt <- "./paper_example/mix_250K_SIF.txt"
## read in array names
scans <- readLines(arrays.txt)[-1]
sif <- read.delim(sif.txt, as.is=TRUE)
library(foreach)
## library(doMC)
## registerDoMC(20)
foreach (scan=scans, .combine=c) %dopar% {
    DoAbsolute(scan, sif)
}

obj.name <- "DRAWS_summary"
results.dir <- file.path(".", "output", "abs_summary")
absolute.files <- file.path(".", 
			    "output", 
			    scans, 
			    "absolute",
			    paste(scans, ".ABSOLUTE.RData", sep=""))

library(ABSOLUTE)
CreateReviewObject(obj.name, 
		   absolute.files, 
		   results.dir, 
		   "allelic", 
		   verbose=TRUE)
## At this point you'd perform your manual review and mark up the file 
## output/abs_summary/DRAWS_summary.PP-calls_tab.txt by prepending a column with
## your desired solution calls. After that (or w/o doing that if you choose to accept
## the defaults, which is what running this code will do) run the following command:
calls.path = file.path("output", "abs_summary", "DRAWS_summary.PP-calls_tab.txt")
modes.path = file.path("output", "abs_summary", "DRAWS_summary.PP-modes.data.RData")
output.path = file.path("output", "abs_extract")
ExtractReviewedResults(calls.path, "test", modes.path, output.path, "absolute", "allelic")
