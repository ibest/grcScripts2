

###### WORKFLOW
# 1. Screen reads for contaminants (minimally phiX)
# 2. deduplicate the reads
# 3. Run sickle (Trimming by quality)
# 4. Flash to join reads
# 5. Kmer normalization
# 6  Merge and report
######
opt=list()
opt$samplesFile = "samples.txt"
opt$samplesColumn = "SAMPLE_ID"

#samples <- loadSamplesFile(opt$samplesFiles,opt$samplesColumn)


