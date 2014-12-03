#!/usr/bin/env Rscript

#*******************************************************************************
# Preprocess Experiment
# preprocesses both fastq files (PE) and SFF files
# performs
# 1. Merging and de-deduplcation of reads files
# 2. 'cleaning' using seqyclean
# 2.a checks for phiX in the sequence and removes if found
# 2.b checks for possible other contaminants or vector sequence if supplied
# 2.c trims sequence ends by quality and plolyA/T (lucy scheme)
# 3. overlaps reads using FLASH (custom)
# 4. produces final read files (SE,PE1, PE2) and reports results
#*******************************************************************************

suppressPackageStartupMessages(library("optparse"))
# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")
option_list <- list(
    make_option(c("-f", "--file"), type="character", default="samples.txt",
                help="The filename of the sample file [default %default]",
                dest="samplesFile"),
    make_option(c("-d", "--directory"), type="character", default="00-RawData",
                help="Directory where the raw sequence data is stored [default %default]",
                dest="Raw_Folder"),
    make_option(c("-q", "--quality"), type="integer", default=24,
                help="Quality score to use during lucy trimming (seqyclean) [default %default]",
                dest="qual"),
    make_option(c("-m", "--minimumLength"), type="integer", default=150,
                help="Discard reads less then minimum length (seqyclean) [default %default]",
                dest="minL"),
    make_option(c("-o", "--overlap"), type="integer", default=700,
                help="maximum overlap (flash) [default %default]",
                dest="overlap"),
    make_option(c("-O", "--skip-overlap"), action="store_true", default=FALSE,
                help="do not perform the overlapping (flash) [default %default]",
                dest="noOverlap"),
    make_option(c("-p", "--processors"), type="integer", default=1,
                help="number of processors to use [default %default]",
                dest="procs"),
    make_option(c("-s", "--skip-duduplicates"), action="store_true", default=FALSE,
                help="do not perform the de-duplication step [default %default] ",
                dest="skip_dedup"),
    make_option(c("-c", "--contaminants-folder"), type="character", default=NULL,
                help="folder name with contaminant sequences in fasta format (seqyclean) [default %default]",
                dest="contaminants"),
    make_option(c("-v", "--vector-folder"), type="character", default=NULL,
                help="folder name with vector sequences in fasta format (seqyclean) [default %default]",
                dest="vector"),
    make_option(c("-a", "--polyA"), action="store_true", default=FALSE,
                help="perform polyA trimming (seqyclean) [default %default]",
                dest="polyA"),
    make_option(c("--i64"), action="store_true",default=FALSE,
                help="input read Q scores are offset by 64 (seqyclean) [default %default]",
                dest="i64")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library("Biostrings"))
suppressPackageStartupMessages(library("ShortRead"))
suppressPackageStartupMessages(library("parallel"))

