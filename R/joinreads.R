# Function to run the application FLASH to overlap paired-end reads
#
# r1 Path to read pair 1 
# r2 Path to read pair 2
# d Path to directory for output files
# t Set number of working threads
# o Output prefix
# O Allow "outies" (read pairs extend in opposite directions joined by overlap) (default:true)
# x Maximum allowed ratio between mismatched base pairs and read length (default:0.25)
# e Minumum overlap length for "outie" oriented reads (default:35bp) 

joinreads <- function(r1,r2,d,O,t,overlap=700){
  paste("flash -O -x 0.25 -e 35 --max-overlap=",overlap," --output-prefix=",o," ",r1," ",r2, " >> ", file.path(d,"preprocessing_output.txt"),
        sep="")
}