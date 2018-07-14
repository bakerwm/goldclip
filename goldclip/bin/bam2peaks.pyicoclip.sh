#!/usr/bin/env bash
##
## call peaks using Pyicoclip
## Date: 2017-08-04 integrate to mapping_clip
## Date: 2017-03-20 initiate 
##
## check regions: TSS, rRNA, pseudo, Promoter, ncRNA, Intron, Intergenic, CDS, 5UTR, 3UTR

##----------------------------------------------------------------------------##
## functions
function echo2() {
    date_fmt='+%Y-%m-%d %H:%M:%S'
    case $2 in 
        error)     echo -e "[$(date "${date_fmt}")] error: $1" && exit 1 ;;
        warning)   echo -e "[$(date "${date_fmt}")] warning: $1" && exit 1 ;;
        *)         echo -e "[$(date "${date_fmt}")] $1" ;;
    esac
}


function homer_annotation_parser() {
    genome=$1 # genome
    ## arguments
    # annoDir="${BINDIR}/../data/${genome}/annotation/homer"
    [[ -z ${annoDir} ]] && echo2 "${genome} - genome not supported." "error"
    tts="${annoDir}/tts.ann.bed"
    rRNA="${annoDir}/rRNA.ann.bed"
    pseudo="${annoDir}/pseudo.ann.bed"
    prom="${annoDir}/promoters.ann.bed"
    ncRNA="${annoDir}/ncRNA.ann.bed"
    utr3="${annoDir}/utr3.ann.bed"
    utr5="${annoDir}/utr5.ann.bed"
    cds="${annoDir}/coding.ann.bed"
    intron="${annoDir}/introns.ann.bed"
    igr="${annoDir}/intergenic.ann.bed"
    ann=(${tts} ${rRNA} ${pseudo} ${prom} ${ncRNA} ${utr3} ${utr5} ${cds} ${intron} ${igr})
    ## check anno files
    check_anno=0
    for a in ${ann[@]} ; do
        [[ -a ${a} ]] && check_anno=$((check_anno+1))
    done
    ## require at least 8 annotation BED files
    [[ ${check_anno} -lt 8 ]] && \
        >&2 echo "[error] - annotation BED files not correct!" && \
        exit 1
    echo ${ann[@]}
}

function run_pyicoclip(){
    bed_in=$1
    region_bed=$2
    peak_out=$3
    pyicoclip_log="${peak_out}.pyicoclip.log"
    pyicoclip --stranded --p-value 0.001 --region ${region_bed} -f bed \
        ${bed_in} ${peak_out} &> ${pyicoclip_log}
    if [[ $(grep -P "Output at:" ${pyicoclip_log} | wc -l) == 1 ]] ; then
        echo 0
    else
        echo 1
    fi
}

## main
[[ $# -lt 3 ]] && echo "Usage: $(basename $0) <genome> <bam> <outdir> " && exit 0
genome=$1
bam_in=$(readlink -f $2)
path_out=$(readlink -f $3)
[[ -f ${bam_in} ]] || echo2 "[${bam_in}] file not exists" "error"
[[ -z $(samtools view -c ${bam_in}) ]] && echo2 "[${bam_in}] file not BAM" "error" # subset of BAM, out-of-MEM, danger
mkdir -p ${path_out} || echo2 "${path_out} - cannot create directory" "error"
prefix=$(basename ${2%.bam}) # original bam_name
path_sub="${path_out}/${prefix}"
path_sub2="${path_sub}/sub_peaks"
path_unfilt="${path_sub2}/unfilt_bedpk"
mkdir -p ${path_unfilt}  || echo2 "${path_unfilt} - cannot create directory" "error" 

BINDIR=$(readlink -fs $(dirname $0))
bed_fix=$(readlink -fs "${BINDIR}/bedFix.py")

## specify GOLDCLIP_PATH
if [ -z ${GOLDCLIP_PATH+x} ] ; then
    echo "\${GOLDCLIP_PATH} not set, please check it throughly"
    exit 1
else
    annoDir="${GOLDCLIP_PATH}/data/${genome}/annotation/homer"
fi

# annotation files
anno=($(homer_annotation_parser ${genome}))
#order: tts rRNA pseudo prom ncRNA intron igr cds utr3 utr5

pushd ${path_sub} &>/dev/null
echo2 "start:Pyicoclip ${prefix}"

# output files
peak_all="${prefix}.all.bed"
peak_uniq="${prefix}.unique.bed"
peak_fixed="${prefix}.fixed.bed"

if [[ -f ${peak_fixed} ]] ; then
    echo2 "file exists, pyicoclip skipped ..."
else
    # bam to bed
    bed_in="${prefix}.bed"
    bedtools bamtobed -i ${bam_in} > ${bed_in}

    # call peaks
    PID=()
    for a in ${anno[@]} ; do
        anno_id=$(basename ${a/.ann.bed})
        peak_out="${path_sub2}/${prefix}.to_${anno_id}.bedpk"
        if [[ -f ${peak_out} ]] ; then
            echo2 "file exists, pyicoclip skipped ... : ${peak_out}"
        else
            run_pyicoclip ${bed_in} ${a} ${peak_out} &>/dev/null &
        fi
        PID+=("$!") # sub process
    done
    wait ${PID[@]} # run in parallel

    # move files to other folder
    for before_bedpk in sub_peaks/before_modfdr*bedpk ; do
        [[ -f ${before_bedpk} ]] && mv ${before_bedpk} sub_peaks/unfilt_bedpk
    done
    for unfilt in sub_peaks/unfilter*bedpk ; do
        [[ -f ${unfilt} ]] && mv ${unfilt} sub_peaks/unfilt_bedpk
    done
    for pylog in sub_peaks/*.log ; do
        [[ -f ${pylog} ]] && mv ${pylog} sub_peaks/unfilt_bedpk
    done

    # # move empty bedpk files
    # for bedpk in sub_peaks/${prefix}*.bedpk ; do
    #     [[ $(file -b ${bedpk}) = empty ]] && mv ${bedpk} sub_peaks/unfilt_bedpk
    # done

    cat sub_peaks/*.bedpk > ${peak_all}
    pyicos remduplicates ${peak_all} ${peak_uniq} &> ${prefix}.remduplicates.log
    python ${bed_fix} ${peak_uniq} > ${peak_fixed}

    # remove temp files
    [[ -f ${peak_all} ]] && rm ${peak_all}
    [[ -f ${bed_in} ]] && rm ${bed_in}

fi
popd >/dev/null 2>&1 #
echo2 "finish:Pyicoclip ${prefix}"

## EOF
