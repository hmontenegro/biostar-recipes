set -ueo pipefail

#Table of contents
INPUTS={{reads.toc}}
LIBRARY={{library.value}}
QUALITY_SCORE={{quality.value}}
WINDOW_SIZE={{window.value}}
ADAPTER={{adapter.value}}
#Output for trimmed files
mkdir -p output

echo ${ADAPTER} >adapter.fa

#if [ "$LIBRARY" == "SE" ] && ["$ADAPTER"==""]; then
#    cat ${INPUTS} | egrep 'fastq|fq' | parallel trimmomatic SE {} ./output/{1/.}_trimmed.fq SLIDINGWINDOW:${WINDOW_SIZE}:${QUALITY_SCORE}
#if [ "$LIBRARY" == "PE" ] && ["$ADAPTER"==""]; then
#    cat ${INPUTS} | sort | egrep 'fastq|fq' | parallel -N 2 "trimmomatic PE {1} {2} ./output/{1/.}paired.fq ./output/{1/.}unpaired.fq ./output/{2/.}paired.fq  ./output/{2/.}unpaired.fq SLIDINGWINDOW:${WINDOW_SIZE}:${QUALITY_SCORE}"
#else
#
#fi

if [ "$ADAPTER" == "DEFAULT" ]; then
    if [ "$LIBRARY" == "SE" ]; then
        cat ${INPUTS} | egrep 'fastq|fq' | parallel trimmomatic SE {} ./output/{1/.}_trimmed.fq SLIDINGWINDOW:${WINDOW_SIZE}:${QUALITY_SCORE}
    fi
    if [ "$LIBRARY" == "PE" ]; then
        cat ${INPUTS} | sort | egrep 'fastq|fq' | parallel -N 2 "trimmomatic PE {1} {2} ./output/{1/.}paired.fq ./output/{1/.}unpaired.fq ./output/{2/.}paired.fq  ./output/{2/.}unpaired.fq SLIDINGWINDOW:${WINDOW_SIZE}:${QUALITY_SCORE}"
    fi
fi
if [ "$ADAPTER" != "DEFAULT" ]; then
    if [ "$LIBRARY" == "SE" ]; then
        cat ${INPUTS} | egrep 'fastq|fq' | parallel trimmomatic SE {} ./output/{1/.}_trimmed.fq ILLUMINACLIP:adapter.fa:2:30:10
    fi
    if [ "$LIBRARY" == "PE" ]; then
        cat ${INPUTS} | sort | egrep 'fastq|fq' | parallel -N 2 "trimmomatic PE {1} {2} ./output/{1/.}paired.fq ./output/{1/.}unpaired.fq ./output/{2/.}paired.fq  ./output/{2/.}unpaired.fq ILLUMINACLIP:adapter.fa:2:30:10"
    fi
fi
