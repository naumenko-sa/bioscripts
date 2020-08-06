library(rtracklayer)
mySession <- browserSession("UCSC")
genome(mySession) <- "hg38"
simpleRepeats <- track(ucscTableQuery(mySession,
                                      track = "Simple Repeats", 
                                      table="simpleRepeat"))
export(simpleRepeats, "hg38_simpleRepeats.bed")
