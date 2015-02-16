###### R functions involved in Preprocessing of Illumina Reads

###### Required Software Versions and Paths
contaminant_path = ""
contaminant_version = "1.0.0"
duplicate_path = ""
duplicate_version = "2.0.0"
sickle_path = "sickle"
sickle_version = "1.33"
flash_path = "flash2"
flash_version = "2.2.00"
######


## Using 'which' make sure appName exists on the path (is executable)
does_app_exist <- function(appName){
  is.null(suppressWarnings(attr(system(paste("which",appName,sep=" "),intern=T,ignore.stdout=T,ignore.stderr=T),"status")))
}

## Compare version numbers
version_check <- function(appVersion, requiredVersion){
  ## assume version number is '.' delimeted
  aV <- strsplit(appVersion,split=".",fixed=T)[[1]]
  rV <- strsplit(requiredVersion,split=".",fixed=T)[[1]]
  for (i in seq.int(1,length(rV))){
    if (is.na(aV[i])) return(FALSE)
    if (aV[i] < rV[i]) return(FALSE)
    if (aV[i] > rV[i]) return(TRUE)
  }
  return(TRUE)
}

#***********************************************************
# Wrappers for preprocessing apps
#***********************************************************
######
# Each app has 4 functions
# app_check - version checker
# app_call - builds the command line string
# app_output - parses the output
# app_run - pulls everything together and runs the app
######

#***********************************************************
# Contaminant Screening (min PhiX)
#***********************************************************

contaminant_check <- function(app=contaminant_path,v=contaminant_version){
  if (!does_app_exist(contaminant_path))
    stop(paste(contaminant_path,"was not found on the path"))
  res <- suppressWarnings(system(paste(contaminant_path,"-v"),intern = TRUE,ignore.stderr=T))
  if (length(res) == 0)
    stop(paste("could not identify version information for",contaminant_path))
  version <- sapply(strsplit(res,split=" +"),"[[",2L)[1]
  return(version_check(version, duplicate_version))
}

#***********************************************************
# Screening of duplicate sequence
#***********************************************************

## first check to make sure the app exists and is the right version
## v minimum version number
screen_duplicates_PE_check <- function(app=duplicate_path,v=duplicate_version){
    if (!does_app_exist(duplicate_path))
      stop(paste(duplicate_path,"was not found on the path"))
    res <- suppressWarnings(system(path(duplicate_path,"--version"),intern = TRUE,ignore.stderr=T))
    if (length(res) == 0)
      stop(paste("could not identify version information for", duplicate_path))
    version <- sapply(strsplit(res,split=" +"),"[[",2L)[1]
    return(version_check(version, duplicate_version))
}

## run the python application screen_duplicates
## d directory containing paired reads
## o output directory
## b start duplicate comparison position
## l length of duplicate comparison
## s skip duplicate detection entirey (effectively just combine reads into a single fastq file)
## a reads were downloaded from the SRA (requires new IDS)
## q quite version mode
screen_duplicates_PE_call <- function(duplicate_path,d,o,b=10,l=25,s=FALSE,a=FALSE,q=FALSE){
  return(paste(duplicate_path,
               "-d", d,
               "-o", o,
               "-b", b,
               "-l", l,
               paste0(
                 ifelse(s,"-s ",""),
                 ifelse(a,"-a ",""),
                 ifelse(q,"--quite ","")),
               sep=" "))
}

## parse the stdout output from screen_duplicates_pe
## output_file file to parse
screen_duplicates_PE_output <- function(output_lines){
  dedup_res <- output_lines[substring(output_lines,1,6) == "Final:"]
  dedup_res <- rev(dedup_res)[1]
  dedup_data <- as.numeric(strsplit(dedup_res,split=" \\| | ")[[1]][seq(3,9,by=2)])
  names(dedup_data) <- c("Pairs","Duplicates","Percent","Reads_Sec")
  return(dedup_data)
}

#screen_duplicates_PE_check()
#dup_out = system(screen_duplicates_PE_call(parameters),intern=TRUE)
#dup_result = screen_duplicates_PE_output(dup_out)

#***********************************************************
# sickle, quality trimmer and min length check
#***********************************************************
sickle_check <- function(app=sickle_path,v=sickle_version){
  if (!does_app_exist(sickle_path))
    stop(paste(sickle_path,"was not found on the path"))
  res <- suppressWarnings(system(paste(sickle_path,"--version"),intern = TRUE,ignore.stderr=T))
  if (length(res) == 0)
    stop(paste("could not identify version information for",sickle_path))
  version <- sapply(strsplit(res[1],split=" +"),"[[",3L)[1]
  return(version_check(version, sickle_version))
}

## run sickle on illumina data
# r1, path to read1
# r2, path to read2
# o, output location
# minL, minimum length
# q, quality trim param
# folder, folder containing contaminant and vector fasta files (contaminants.fa and/or vector.fa)
# sample, sample name
# i64, is data in illumina 64 format (CASAVA 1.3-1.7), otherwise sanger (CASAVA 1.8)
sickle_pe_call <- function(sickle_path,r1,r2,o,minL=150,q=24,i64=FALSE) {
    return(paste(sickle_path,"pe",
                 "-g",
                 "-t",ifelse(i64,"illumina","sanger"),
                 "-q", q,
                 "-l", minL,
                 "-f", r1,
                 "-r", r2,
                 "-o", paste0(o,"_PE1.fastq.gz"),
                 "-p", paste0(o,"_PE2.fastq.gz"),
                 "-s", paste0(o,"_SE.fastq.gz"),
                 sep=" "))
}

sickle_output <- function(output_lines){
  analyzed <- output_lines[substring(output_lines,1,26) == "Total input FastQ records:"]
  analyzed <- as.numeric(sapply(strsplit(gsub(",|\\(|\\)","",analyzed),split=" "),"[[",6L))
  discard_pair <- output_lines[substring(output_lines,1,31) == "FastQ paired records discarded:"]
  discard_pair <- as.numeric(sapply(strsplit(gsub(",|\\(|\\)","",discard_pair),split=" "),"[[",6L))
  discard_single <- output_lines[substring(output_lines,1,31) == "FastQ single records discarded:"]
  discard_single <- as.numeric(sapply(strsplit(gsub(",|\\(|\\)","",discard_single),split=" "),"[",c(8L,11L)))

  sickle_data <- c(analyzed,discard_single,discard_pair,(analyzed-sum(discard_single,discard_pair)),sum(discard_single),(analyzed-sum(discard_single,discard_pair))/analyzed,(2*analyzed-sum(discard_single,2*discard_pair))/(2*analyzed))
  names(sickle_data) <- c("Pairs","PE1_discarded","PE2_discarded","Both_discarded","Pairs_kept","Single_kept","Pairs_kept_percent","Reads_kept_percent")
  return(sickle_data)

}

#sickle_check()
#sickle_out = system(sickle_pe_call(params),intern=TRUE)
#sickle_result = sickle_output(sickle_out)

#sickle_output = system(sickle_pe_call("Wireworm_nodup_noContaminate_PE1.fastq.gz","Wireworm_nodup_noContaminate_PE2.fastq.gz",o="Wireworm_nodub_noContaminant_sickle"),inter=T)


#***********************************************************
# Flash2
#***********************************************************
flash2_check <- function(app=flash_path,v=flash_version){
  if (!does_app_exist(flash_path))
    stop(paste(flash_path,"was not found on the path"))
  res <- suppressWarnings(system(paste(flash_path,"--version"),intern = TRUE,ignore.stderr=T))
  if (length(res) == 0)
    stop(paste("could not identify version information for",flash_path))
  version <- sub("v","",sapply(strsplit(res[1],split=" +"),"[[",2L)[1])
  return(version_check(version, flash_version))
}

flash2_call <- function(flash_path,r1,r2,o,discard=TRUE,allow_outies=TRUE,max_overlap=700,min_overlap_outie=35,min_overlap_innie=10,mismatch_density=0.25,i64=FALSE,threads=1){
  return(paste(flash_path,
               "-z",
               "-t", threads,
               "-p",ifelse(i64,"64","33"),
               paste0(
                 ifelse(discard,"","-D"),
                 ifelse(allow_outies,"-O","")),
               "-x",mismatch_density,
               "-e",min_overlap_outie,
               "-m",min_overlap_innie,
               "-M",max_overlap,
               "-o", o,
               r1,
               r2,
               sep=" "))
}

flash2_output <- function(output_lines){
  oflash <- tail(which(output_lines=="[FLASH] Read combination statistics:"),1)
  flash_res <- output_lines[c((oflash+1),(oflash+3),(oflash+2),(oflash+4),(oflash+5))]
  flash_data <- as.integer(sapply(strsplit(flash_res,split=" +|%"),"[[",4L))
  names(flash_data) <- c("Pairs","Reads_Discarded","Reads_Combined","Innie","Outie")
  return(flash_data)
}

#flash2_check()
#flash2_out = system(flash2_call(params),intern=TRUE)
#flash2_result = flash2_output(flash2_out)

#flash2_out <- system(flash2_call("Wireworm_nodup_noContaminant_sickle_PE1.fastq.gz","Wireworm_nodup_noContaminant_sickle_PE2.fastq.gz","Wireworm_nodup_noContaminant_sickle_flash",threads= 40),intern=TRUE)


#***********************************************************
# Seqyclean
#***********************************************************
## run seqyclean on illumina data
# r1, path to read1
# r2, path to read2
# o, output location
# minL, minimum length
# q, quality trim param
# folder, folder containing contaminant and vector fasta files (contaminants.fa and/or vector.fa)
# sample, sample name
# i64, is data in illumina 64 format
seqyclean_illumina <- function(r1,r2,o,minL=150, q=24,folder, sample, i64) {
  i64_param = ""
  if (i64){
    i64_param="-i64"
  }
  vc_param = ""
  if(file.exists(file.path(folder,"contaminants.fa"))){
    vc_param=paste(vc_param,"-c",file.path(getwd(),folder,"contaminants.fa"),sep=" ")
  }
  if(file.exists(file.path(folder,"vector.fa"))){
    vc_param=paste(vc_param,"-v",file.path(getwd(),folder,"vector.fa"),sep=" ")
  }
  paste("seqyclean --ow -qual", q, q, i64_param, vc_param,"-minimum_read_length",minL,"--new2old_illumina -1",r1,"-2",r2,"-o",o, ">>", file.path(folder,sample,"preprocessing_output.txt"),sep=" ")
}
## run the application FLASH to join overlapping paired end reads

join_reads <- function(r1,r2,o,overlap=275,d){
  paste("flash --max-overlap=",overlap," --output-prefix=",o," ",r1," ",r2, " >> ", file.path(d,"preprocessing_output.txt"),sep="")
}

## generate final linked files in cleaned folder with report
link_illumina <- function(se1,se2,r1,r2,o){
  require("ShortRead")
  ## first merge the 2 SE files
  se <- file.path(dirname(se1),"merged_SE_files.fastq")
  if (!file.exists(se)){
    fq <- readFastq(c(se1,se2))
    writeFastq(fq,se)
  }
  paste("ln -sf",file.path("../..",se),paste(o,"merged_SE.fastq",sep="_"),";ln -sf",file.path("../..",r1),paste(o,"notcombined_PE1.fastq",sep="_"),";ln -sf",file.path("../..",r2),paste(o,"notcombined_PE2.fastq",sep="_"),";",sep=" ")
}

## run seqyclean on 454 data
# path to sff file
# o, output location
# minL, minimum length
# q, quality trim param
# folder, folder containing contaminant and vector fasta files (contaminants.fa and/or vector.fa)
# sample, sample name
seqyclean_454 <- function(sff,o,minL=225,q=24,folder,sample){
  vc_param = ""
  if(file.exists(file.path(folder,"contaminants.fa"))){
    vc_param=paste(vc_param,"-c",file.path(getwd(),folder,"contaminants.fa"),sep=" ")
  }
  if(file.exists(file.path(folder,"vector.fa"))){
    vc_param=paste(vc_param,"-v",file.path(getwd(),folder,"vector.fa"),sep=" ")
  }
  paste("seqyclean -qual",q,q,vc_param,"-minimum_read_length",minL,"-454",sff,"-o",o, ">>", file.path(folder,sample,"preprocessing_output_454.txt"),sep=" ")
}
link_454 <- function(sff,o){
  paste("ln -sf",sff,paste(o,"sff",sep="."),sep=" ")
}

final_report_fun <- function(f,o){
  paste("read_info.py", "-d",f,">>",file.path(o,"preprocessing_output_final_report.txt"))
}

## download phiX from NCBI for screen Illumina data, always perform
get_phiX <- function(genbank_id = "NC_001422"){
  URL <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
               genbank_id, "&rettype=fasta&retmode=text",
               sep = "")
  res <- scan(file = URL, what = "", sep = "\n", quiet = TRUE)
  phix <- DNAStringSet(paste(res[-1],collapse=""))
  names(phix) <- "contaminant_PhiX"
  return(phix)
}

##### Actual Preprocessing Workflow


