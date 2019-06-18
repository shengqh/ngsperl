#!/bin/bash
set -e
# set -o pipefail
version=0.1
############################################################
#  Program: sv_pipeline
#  Version: 0.1
#  Author: Youhuang Bai
############################################################
SV_DIR=`dirname $0`
LC_ALL=C
## global usage
function usage() {
    echo "
Program: sv_pipeline
Version: 0.1
Author: Youhuang Bai (youhuang.bai.1@vumc.org)
usage:   sv_pipeline <command> [options]
command: mapping     align FASTQ files with BWA-MEM, and add MC MQ Tags with samblaster
         refine      name sort the bam and apply gatk recalibration to it 
         call        SV call with lumpy in smoove.
         merge       merge SVs with smoove
         genotype    genotype SVs with svtyper in smoove
         paste       paste all samples into one file
         annotate    annotate with genome annotation
         filter      filter deletion and duplication
options: -h       show this message
"
}

function mapping() {
    function mapping_usage() {
	echo "
usage:   SV_pipeline mapping -g <reference.fa> -f1 <in1.fq> -f2 [in2.fq]
"
    }
    
    # Check options passed in.
    if test -z "$2"
    then
        mapping_usage
    exit 1
    fi
    
    # set defaults
    
    ##parameters for SLURM
    # nodes=1
    # nodepertasks=8
    # ntasks=$[ $nodes*$nodepertasks ]
    # ntime="336:00:00"
    # mem=$[ 6*($nodepertasks+1) ]
    # singularity=0; # 1: ues singularity; 0 do not use singularity
    # HG="hg38"
    
    while getopts ":hg:T:t:R:O:o:pV:" OPTION
    do
        case "${OPTION}" in
            h)
                mapping_usage
                exit 1
                ;;
            g)
                GENOME="$OPTARG"
                ;;
            t)
                THREADS="$OPTARG"
                ;;
            R)
                READ_GROUP="$OPTARG"
                ;;
            O)
                WORK_DIR="$OPTARG"
                ;;
            o)
                OUTPUT="$OPTARG"
                ;;
            p)
                INTERLEAVED=1
                ;;
            V)
                HG="$OPTARG"
                ;;
        esac
    done
    
    if [ $HG == "hg19" ]; then
        outdir="1_mapping_hg19"
    elif [ $HG == "hg38" ]; then
        outdir="1_mapping_hg38"
    else
        echo "Unknown genome version！"
        exit 1
    fi
    
    if [[ "$INTERLEAVED" -eq 1 ]]; then
        REF="${@:${OPTIND}:1}"
        FQ="${@:$((${OPTIND}+1)):1}"
        if [[ -z "$WORK_DIR" ]]; then
            OUTPUT=`basename "$FQ"`
        fi
        if [[ -z "$FQ" ]]; then
            mapping_usage
            echo -e "Error: Fastq file $FQ not found.\n"
            exit 1
        fi
    else
        REF="${@:${OPTIND}:1}"
        FQ1="${@:$((${OPTIND}+1)):1}"
        FQ2="${@:$((${OPTIND}+2)):1}"
        if [[ -z "$OUTPUT" ]]; then
        OUTPUT=`basename "$FQ1"`
        fi
        if [[ -z "$FQ1" ]]; then
            mapping_usage
            echo -e "Error: Fastq file $FQ1 not found.\n"
            exit 1
        elif [[ -z "$FQ2" ]]; then
            mapping_usage
            echo -e "Error: Fastq file $FQ2 not found. (single-end reads not supported, use -p for interleaved FASTQ)\n"
            exit 1
        fi
    fi

    if [[ -z "$GENOME" ]] || [[ ! -f "$GENOME" ]]; then
    mapping_usage
    echo -e "Error: Reference file $GENOME not found.\n"
    exit 1
    fi

    # Check for readgroup flag
    if [[ -z $READ_GROUP ]]
    then
    mapping_usage
    echo -e "Error: no readgroup found. Please set a readgroup with the -R flag.\n"
    exit 1
    fi
    
    if [[ "$INTERLEAVED" -eq 1 ]]; then
        bwa mem -K 100000000 -Y -t $THREADS -R $READ_GROUP $GENOME $FQ | samblaster -a --addMateTags | samtools view -Sb - >${OUTPUT}
    else
        bwa mem -K 100000000 -Y -t $THREADS -R $READ_GROUP $GENOME $FQ1 $FQ2 | samblaster -a --addMateTags | samtools view -Sb - >${OUTPUT}
    fi
}

function refine() {
    function refine_usage() {
	echo "
usage:   SV_pipeline refine -g <reference.fa> -f1 <in1.fq> -f2 [in2.fq]
"
    }
    
    # Check options passed in.
    if test -z "$2"
    then
        refine_usage
    exit 1
    fi

    # set defaults
    
    ##parameters for SLURM
    # nodes=1
    # nodepertasks=8
    ntasks=4
    # ntime="336:00:00"
    mem=170
    # singularity=0; # 1: ues singularity; 0 do not use singularity
    # HG="hg38"
    
    while getopts ":hg:T:t:R:O:o:pV:M:" OPTION
    do
        case "${OPTION}" in
            h)
                refine_usage
                exit 1
                ;;
            g)
                GENOME="$OPTARG"
                ;;
            t)
                THREADS="$OPTARG"
                ;;
            O)
                WORK_DIR="$OPTARG"
                ;;
            M)
                TMP_DIR="$OPTARG"
                ;;
            o)
                OUTPUT="$OPTARG"
                ;;
            P)
                PATTERN="$OPTARG"
                ;;
            V)
                HG="$OPTARG"
                ;;
        esac
    done
    if test -z "$PATTERN"; then
        echo "target bam pattern needed to be provied for merge！"
        exit 1
    elif [[ -d $WORK_DIR/1_mapping_${HG}/]];then
        line=${PATTERN%?}
    else
        echo "could not find previous mapping directory!"
        exit 1
    fi
    
    if [ $HG == "hg19" ]; then
        outdir="2_refine_hg19"
    elif [ $HG == "hg38" ]; then
        outdir="2_refine_hg38"
    else
        echo "Unknown genome version！"
        exit 1
    fi
    
    
    if [[ -z "$GENOME" ]] || [[ ! -f "$GENOME" ]]; then
    mapping_usage
    echo -e "Error: Reference file $GENOME not found.\n"
    exit 1
    fi

    if [[ -z $TMP_DIR]]
    then
    TMP_DIR=`mktemp -d /tmp/SV.XXXXXXXXXXXX`
    fi

mkdir $TMP_DIR
date
ls $WORK_DIR/1_mapping_${HG}/${pattern}|xargs samtools merge -u -r $TMP_DIR/${line}.bam 
sambamba sort -u --sort-picard -t $ntasks -m ${mem}G -o $TMP_DIR/${line}.sorted.queryName.bam $TMP_DIR/${line}.bam
rm $TMP_DIR/${line}.bam
date
java -Xmx$((${mem}-4))g -jar /scratch/cqs/shengq2/local/bin/picard/picard.jar MarkDuplicates ASSUME_SORT_ORDER="queryname" COMPRESSION_LEVEL=0 I=$TMP_DIR/${line}.sorted.queryName.bam O=$TMP_DIR/${line}.mkDup.bam M=$TMP_DIR/${line}.mkDup.metrics.txt 
rm $TMP_DIR/${line}.sorted.queryName.bam
date
sambamba sort -t $ntasks -m ${mem}G -o $WORK_DIR/$outdir/result/$line/${line}.mkDup.sorted.bam $TMP_DIR/${line}.mkDup.bam

date
java -Xmx$(((${mem}-4)/4))g -jar /scratch/cqs/softwares/gatk3.jar -T BaseRecalibrator \
-nct 4 -R $GENOME -I $WORK_DIR/$outdir/result/$line/${line}.mkDup.sorted.bam \
-o $WORK_DIR/$outdir/result/$line/${line}.table \
-knownSites /scratch/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
-knownSites /scratch/cqs/references/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
-knownSites /scratch/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz \
--preserve_qscores_less_than 6 -dfrac .1 \
-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22
date
java -Xmx40g -jar /scratch/cqs/softwares/gatk3.jar -T PrintReads -nct 4 \\
-I $WORK_DIR/$outdir/result/$line/${line}.mkDup.sorted.bam \\
-o $WORK_DIR/$outdir/result/$line/${line}.mkDup.sorted.recal.bam \\
-BQSR $WORK_DIR/$outdir/result/$line/${line}.table \\
-R $GENOME \\
-preserveQ 6 \\
-SQQ 10 -SQQ 20 -SQQ 30 \\
--disable_indel_quals

rm -rf $TMP_DIR
date

}

function merge() {
    function refine_usage() {
	echo "
usage:   SV_pipeline merge -g <reference.fa> -f1 <in1.fq> -f2 [in2.fq]
"
    }
    
    # Check options passed in.
    if test -z "$2"
    then
        refine_usage
    exit 1
    fi

    # set defaults
    
    ##parameters for SLURM
    # nodes=1
    # nodepertasks=8
    ntasks=4
    # ntime="336:00:00"
    mem=170
    # singularity=0; # 1: ues singularity; 0 do not use singularity
    # HG="hg38"
    
    while getopts ":hg:T:t:R:O:o:pV:M:" OPTION
    do
        case "${OPTION}" in
            h)
                refine_usage
                exit 1
                ;;
            g)
                GENOME="$OPTARG"
                ;;
            t)
                THREADS="$OPTARG"
                ;;
            O)
                WORK_DIR="$OPTARG"
                ;;
            M)
                TMP_DIR="$OPTARG"
                ;;
            o)
                OUTPUT="$OPTARG"
                ;;
            P)
                PATTERN="$OPTARG"
                ;;
            V)
                HG="$OPTARG"
                ;;
        esac
    done
    if test -z "$PATTERN"; then
        echo "target bam pattern needed to be provied for merge！"
        exit 1
    elif [[ -d $WORK_DIR/1_mapping_${HG}/]];then
        line=${PATTERN%?}
    else
        echo "could not find previous mapping directory!"
        exit 1
    fi
    
    if [ $HG == "hg19" ]; then
        outdir="2_refine_hg19"
    elif [ $HG == "hg38" ]; then
        outdir="2_refine_hg38"
    else
        echo "Unknown genome version！"
        exit 1
    fi
    
    
    if [[ -z "$GENOME" ]] || [[ ! -f "$GENOME" ]]; then
    mapping_usage
    echo -e "Error: Reference file $GENOME not found.\n"
    exit 1
    fi



}