set -ueo pipefail

#Table of contents
INPUTS={{reads.toc}}
LIBRARY={{library.value}}
QUALITY_SCORE={{quality.value}}
WINDOW_SIZE={{window.value}}
ADAPTER={{adapter.value}}


# The directory that contains the preliminary reports.
mkdir -p fastqc_reports
# Add extra flags if the group parameter was selected.
if [ "{{group.value}}" == "nogroup" ]; then
    FLAGS="--nogroup"
else
    FLAGS=""
fi
# Run fastqc on selected files.
cat ${INPUTS} | egrep 'fastq|fq' | parallel fastqc -q ${FLAGS} {} -o fastqc_reports



#Output for trimmed files
mkdir -p output
#Adapter sequence fasta file for illumina clip
printf ">adapter sequence \n"%s"" ${ADAPTER} >adapter.fa


if [ "$ADAPTER" == "No_adapter" ] && [ "$LIBRARY" == "SE" ]; then
    cat ${INPUTS} | egrep 'fastq|fq' | parallel trimmomatic SE {} ./output/{1/.}_trimmed.fq SLIDINGWINDOW:${WINDOW_SIZE}:${QUALITY_SCORE}
fi

if [ "$ADAPTER" == "No_adapter" ] && [ "$LIBRARY" == "PE" ]; then
    cat ${INPUTS} | sort | egrep 'fastq|fq' | parallel -N 2 "trimmomatic PE {1} {2} ./output/{1/.}paired.fq ./output/{1/.}unpaired.fq ./output/{2/.}paired.fq  ./output/{2/.}unpaired.fq SLIDINGWINDOW:${WINDOW_SIZE}:${QUALITY_SCORE}"
fi


if [ "$ADAPTER" != "No_adapter" ] &&  [ "$LIBRARY" == "SE" ]; then
    cat ${INPUTS} | egrep 'fastq|fq' | parallel trimmomatic SE {} ./output/{1/.}_trimmed.fq ILLUMINACLIP:adapter.fa:2:30:10
fi


if [ "$ADAPTER" != "No_adapter" ] &&  [ "$LIBRARY" == "PE" ]; then
    cat ${INPUTS} | sort | egrep 'fastq|fq' | parallel -N 2 "trimmomatic PE {1} {2} ./output/{1/.}_paired.fq ./output/{1/.}_unpaired.fq ./output/{2/.}paired.fq  ./output/{2/.}unpaired.fq ILLUMINACLIP:adapter.fa:2:30:10"
fi


#Fastqc on the output
cd output

#file list of the trimmed files
ls > file_list.txt

# Report directory.
mkdir -p fastqc_reports_after_trimming

#Running fastqc on reported files
cat file_list.txt | egrep 'fastq|fq' | parallel fastqc -q ${FLAGS} {} -o fastqc_reports_after_trimming



