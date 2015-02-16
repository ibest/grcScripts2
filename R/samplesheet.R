
file_extension <- function (x, compression=TRUE)
{
  if (compression)
    x <- sub('.gz$|.bz2$|.zip$|.xz$', '',x)
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

multisub <- function (pattern, replacement, text.var)
{  
  if (length(replacement) == 1){
    replacement <- rep(replacement, length(pattern))
  }else if (length(replacement) != length(pattern)){
    write(paste("Lengths of input strings differ (and replacement is not 1). Please check your input.\n"), stderr())
    stop("Quiting")
  }
  for (i in seq_along(pattern)){
    text.var <- sub(pattern[i], replacement[i], text.var, fixed = FALSE)
  }
  text.var
}

pair_reads <- function(files,pair.patterns=c("_R1|_R2","READ1|READ2","_1.fastq|_2.fastq","_PE1|_PE2"))
{
  pairs <- cbind(pair.patterns,paste0("@",seq_along(pair.patterns),"@"))
  red.pairs <- multisub(pairs[,1],pairs[,2],files)
  tb_red.pairs <- table(red.pairs)
  if (any(!(tb_red.pairs %in% c(1,2)))){
    write(paste("Something other than single-end or paired-end files found\n"), stderr())
    stop("Quiting")
  }
  SE.files <- names(tb_red.pairs)[which(tb_red.pairs == 1)]
  PE.files <- lapply(names(tb_red.pairs)[which(tb_red.pairs == 2)], function(file){
    ptype <- which(sapply(pairs[,2],grep,file)==1)
    unname(sapply(unname(unlist(strsplit(pairs[ptype,1],split="|",fixed=TRUE))),sub,pattern=pairs[ptype,2],x=file,fixed=TRUE))
  })
  return(list(SE=SE.files,PE=PE.files))
}

######################################################################
## loadSampleFile
## reads in the sample sheet, check for minimum expected format (minimum scolumn or 'SAMPLE_ID' in header)
## Parameters
##  file: sample sheet filename, scolumn column name for the sample ids, tab delimited
"loadSamplesFile" <- function(file="samples.txt",scolumn="SAMPLE_ID"){
  ##
  if ( !file.exists(file) ) {
    write(paste("Sample file",file,"does not exist\n"), stderr())
    stop("Quiting")
  }  
  ### rows can be commented out with #
  targets <- read.table(file,sep="\t",header=TRUE,as.is=TRUE)
  if( !(scolumn %in% colnames(targets)) ){
    write(paste("Expecting", scolumn,"in the header of samples file\n"), stderr())
    stop("Quiting")
  }
  samples <- apply(targets,1,function(x) { list(name=x[scolumn],metadata=x)})
  names(samples) <- targets$SAMPLE_ID
  write(paste("Samples sheet contains", length(samples), "samples to process"),stdout())
  if (length(samples) == 0)
    stop("Quiting")
  return(samples)  
}

######################################################################
## loadDataFiles
## provided the sample sheet,
## then check to make sure data (sequence reads, bam files, etc.) are available using the specified column
## Parameters
##  samples: Result of loadSamplesFile
##	data_folder: path to folder containing reads
##  type: if processing raw data, folderids use a differnt column identifier 'raw'
##  rcolumn: sample sheet column to use that specified raw data folders
"loadDataFiles" <- function(samples, data_folder="00-RawData",type=c("raw","fastq","fasta","bam"),rcolumn="SEQUENCE_ID"){
  ###
  type = match.arg(type)
  if(!(data_folder %in% dir())){
    write(paste("Data folder:", data_folder,", not found\n"), stderr())
    stop("Quiting")    
  }
  ### RAW READS, USES INFO FROM A DIFFERENT COLUMN IN THE SAMPLESHEET
  if (type=="raw"){
    if(!(rcolumn %in% names(samples[[1]][["metadata"]]))){
      write(paste("Expecting", rcolumn,"in the header of samples file\n"), stderr())
      stop("Quiting")
    }
    if (sum(sapply(samples,"[[","metadata")[rcolumn,] %in% dir(path=data_folder),na.rm=TRUE) != length(samples)){
      write(paste(rcolumn,"does not match the",data_folder,"data folder structure\n\n"), stderr())
      write(paste("SAMPLE",rcolumn,"FOUND\n",sep="\t"),stderr())
      write(paste(apply(data.frame(names(samples),sapply(samples,"[[","metadata")[rcolumn,],sapply(samples,"[[","metadata")[rcolumn,] %in% dir(path=data_folder)),1,paste,collapse="\t"),collapse="\n"),stderr())
      write("\n",stderr())
      stop("Quitting")
    }
    samples <- lapply(samples,function(sample){
      if (!file.info(file.path(data_folder,sample$metadata[rcolumn]))$isdir){
        write(paste("Data should be stored in subfolders of",data_folder,"\n"), stderr())
        stop("Quiting")
      }
      files <- dir(file.path(data_folder,sample$metadata[rcolumn]),full.names=TRUE)
      ### only process fastq files
      files <- files[file_extension(files) == "fastq"]
      sample$reads <- pair_reads(files)
      sample
    })
    return(samples)
  }
  ### FASTQ Reads in SAMPLE FOLDERS
  if(type=="fastq"){
    if (sum(sapply(samples,"[[","name") %in% dir(path=data_folder),na.rm=TRUE) != length(samples)){
      write(paste("samples file does not match the",data_folder,"data folder structure\n\n"), stderr())
      write(paste("Sample","FOUND\n",sep="\t"),stderr())
      write(paste(apply(data.frame(sapply(samples,"[[","name"),sapply(samples,"[[","name") %in% dir(path=data_folder)),1,paste,collapse="\t"),collapse="\n"),stderr())
      stop("Quitting")
    }
    samples <- lapply(samples,function(sample){
      if (!file.info(file.path(data_folder,sample$name))$isdir){
        write(paste("Data should be stored in subfolders of",data_folder,"\n"), stderr())
        stop("Quiting")
      }
      files <- dir(file.path(data_folder,sample$name),full.names=TRUE)
      ### only process fastq files
      files <- files[file_extension(files) == "fastq"]
      sample$reads <- pair_reads(files)
      sample
    })
    return(samples)    
  }
  ### FASTA Reads in SAMPLE FOLDERS
  if(type=="fasta"){
    if (sum(sapply(samples,"[[","name") %in% dir(path=data_folder),na.rm=TRUE) != length(samples)){
      write(paste("samples file does not match the",data_folder,"data folder structure\n\n"), stderr())
      write(paste("Sample","FOUND\n",sep="\t"),stderr())
      write(paste(apply(data.frame(sapply(samples,"[[","name"),sapply(samples,"[[","name") %in% dir(path=data_folder)),1,paste,collapse="\t"),collapse="\n"),stderr())
      stop("Quitting")
    }
    samples <- lapply(samples,function(sample){
      if (!file.info(file.path(data_folder,sample$name))$isdir){
        write(paste("Data should be stored in subfolders of",data_folder,"\n"), stderr())
        stop("Quiting")
      }
      files <- dir(file.path(data_folder,sample$name),full.names=TRUE)
      ### only process fasta files
      files <- files[file_extension(files) == "fasta"]
      sample$reads <- pair_reads(files)
      sample
    })
    return(samples)        
  }  
  ### BAM files in SAMPLE FOLDERS
  if(type=="bam"){
    write(paste("Still need to write bam detector\n"), stderr())
    stop("Quiting")
  }
}
