#!/bin/bash

# A POSIX variable
OPTIND=1  # Reset in case getopts has been used previously in the shell.

# Initialize variables:
threads=4
pe1=NA
pe2=NA
se="-"


while getopts "h?t:c:e:o:1:2:U:" opt; do
    case "$opt" in
    h|\?)
        echo "there is no help"
        exit 0
        ;;
    t)  threads=$OPTARG
        ;;
    c)  contaminant=$OPTARG
        ;;
    e)  extract_path=$OPTARG
        ;;
    o)  output=$OPTARG
        ;;
    1)  pe1=$OPTARG
        ;;
    2)  pe2=$OPTARG
        ;;
    U)  se=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

btargs=$@

wget -O phiX.tmp.fasta "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=9626372"
if [ $? -eq 0 ] ; then
 cat $contaminant phiX.tmp.tmp > contaminant_screen.fasta
 bowtie2-build contaminant_screen.fasta contaminant_screen_index
 bowtie_index=contaminant_screen_index
else
 echo "ERROR: No copy of phiX sequence available, and download from NCBI failed. Goodbye"
 exit 1
fi

#build the command string
bt_command="bowtie2 -I 0 -X 1500 --very-sensitive-local -p $threads -x $bowtie_index -q "
if [ $se = "-" ] ; then # expect paired end reads
      suffix=" -1 $pe1 -2 $pe2"
else
      suffix=" -U $se"
fi
source /usr/modules/init/bash
module load bowtie2 grc/2.0

#echo "The command:"
#echo "$bt_command $suffix | samtools view -bS "
#$bt_command $suffix | samtools view -bS - > btout.bam

$bt_command $suffix | $extract_path -o $output

# clean up
if [ -f contaminant_screen_index.fasta ] ; then rm phiX.tmp.fasta contaminant_screen_index*; fi

exit 0
