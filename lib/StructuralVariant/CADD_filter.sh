#/usr/bin/bash

set -e
unset PERL5LIB
# set -o pipefail


    function usage() {
    echo "
usage:  CADD_filter.sh -V <vcf.file> -T <./test> 
version: 0.1
        -h    help information
        -V    vcf file
        -t    cadd threshold value
        -T    TEMP dir
"
    }
    
    # Check options passed in.
    if test -z "$1"
    then
        usage
    exit 1
    fi
    
    
    while getopts ":hV:g:t:T:S:" OPTION
    do
        case "${OPTION}" in
            h)
                usage
                exit 1
                ;;
            V)
                VCF="$OPTARG"
                ;;
            g)
                hg="$OPTARG"
                ;;
            v)
                version="$OPTARG"
                ;;
            t)
                THRESHOLD="$OPTARG"
                ;;
            T)
                TMP_DIR="$OPTARG"
                ;;
            S)
                split=true
                ;;
        esac
    done
    if [[ -z $VCF ]]
    then
    echo "VCF file was not provided!"
    fi
    if [[ -z $hg ]]
    then
    hg="GRCh38"
    fi
    if [[ -z $version ]]
    then
    version="v1.5"
    fi
    if [[ -z $split ]]
    then
    split=true
    fi
    if [[ -z $THRESHOLD ]]
    then
    THRESHOLD=10
    fi
    if [[ -z $TMP_DIR ]]
    then
        TMP_DIR=`mktemp -d /tmp/CADD_filter.XXXXXXXXXXXX`
    else
        if [[ ! -d $TMP_DIR ]]; then 
            mkdir $TMP_DIR
        fi
    fi

echo $VCF
dir=$(dirname $VCF)
if [[ $VCF == *.vcf ]]; then 
    decomp=false
    file=$(basename $VCF .vcf)
elif [[ $VCF == *.vcf.gz ]]; then 
    decomp=true
    file=$(basename $VCF .vcf.gz)
else
    echo "Does your file end with .vcf or .vcf.gz? If not, please change it!"
    exit 1
fi

date

if [[ ! "$split" = true ]];then
    if [[ $hg == "GRCh38" ]]; then
        if [[ $decomp = true ]];then
            zcat $VCF | sed 's/^chr//' >$TMP_DIR/${file}_nochr.vcf
        else
            sed 's/^chr//' $VCF >$TMP_DIR/${file}_nochr.vcf
        fi
    else 
        ln -s $VCF $TMP_DIR/${file}_nochr.vcf
    fi
    /scratch/cqs/baiy7/tools/CADD-scripts/CADD.sh -g $hg -v $version -o ${dir}/${file}_nochr_CADD${version}.vcf ${TMP_DIR}/${file}_nochr.vcf
    LANG=c grep "^#" $VCF >${TMP_DIR}/header.vcf
    LANG=c grep -v "^#" $VCF >${TMP_DIR}/body.vcf
    zcat ${dir}/${file}_nochr_CADD${version}.vcf|LANG=c grep -v "^#" |awk -F '\t' '{ if($6 >$THRESHOLD) print NR}' >${TMP_DIR}/line_num.txt
    awk 'NR == FNR{a[$0]; next};FNR in a' ${TMP_DIR}/line_num.txt ${TMP_DIR}/body.vcf >${TMP_DIR}/body_fitler.vcf
    cat ${TMP_DIR}/header.vcf ${TMP_DIR}/body_fitler.vcf >${dir}/${file}_CADD${version}_filtered${THRESHOLD}.vcf
    #rm -rf $TMP_DIR
else
    if [[ $hg == "GRCh38" ]]; then
        if [[ $decomp = true ]]; then 
            if [[ ! -f $VCF.tbi ]]; then
                tabix -p vcf $VCF
            fi
            zcat $VCF |sed -n '/^[^#]/q;p' >$TMP_DIR/header.vcf
	    zcat $VCF |LANG=c grep -v "^#" >$TMP_DIR/body.vcf
	    zcat $VCF |sed 's/^chr//' >$TMP_DIR/${file}_nochr.vcf
            bgzip -c $TMP_DIR/${file}_nochr.vcf >$TMP_DIR/${file}_nochr.vcf.gz
            tabix -p vcf $TMP_DIR/${file}_nochr.vcf.gz
        else
            sed -n '/^[^#]/q;p' $VCF >$TMP_DIR/header.vcf
            LANG=c grep -v "^#" $VCF >$TMP_DIR/body.vcf
            if [ ! -s $TMP_DIR/${file}_nochr.vcf ]; then
                sed 's/^chr//' $VCF >$TMP_DIR/${file}_nochr.vcf
            fi
            if [ ! -s $TMP_DIR/${file}_nochr.vcf.gz ]; then
                bgzip -c $TMP_DIR/${file}_nochr.vcf >$TMP_DIR/${file}_nochr.vcf.gz
                tabix -p vcf $TMP_DIR/${file}_nochr.vcf.gz
            fi
        fi
      
        parallel -j8 -N1 tabix $TMP_DIR/${file}_nochr.vcf.gz {} \> $TMP_DIR/{}.vcf ::: $(tabix -l $TMP_DIR/${file}_nochr.vcf.gz)
        parallel -j8 -N1 /scratch/cqs/baiy7/tools/CADD-scripts/CADD.sh -o $TMP_DIR/{}_CADD${version}.vcf.gz $TMP_DIR/{}.vcf ::: $(tabix -l $TMP_DIR/${file}_nochr.vcf.gz)
        parallel -j8 -N1 gzip -d $TMP_DIR/{}_CADD${version}.vcf.gz ::: $(tabix -l $TMP_DIR/${file}_nochr.vcf.gz)
        ##get each line_number files
	if [ ! -s $TMP_DIR/${file}_CADD${version}.vcf ]; then
	        echo "1"
		tabix -l $TMP_DIR/${file}_nochr.vcf.gz |xargs  -I {} grep -v '^#' $TMP_DIR/{}_CADD${version}.vcf |cat > $TMP_DIR/${file}_CADD${version}.vcf
	fi
	if [ ! -s ${TMP_DIR}/line_num.txt ];then
		echo "2"	
		awk -F '\t' "{ if(\$6 >$THRESHOLD) print NR}" $TMP_DIR/${file}_CADD${version}.vcf > ${TMP_DIR}/line_num.txt
	fi
        ##get each body.vcf files
        awk 'NR == FNR{a[$0]; next};FNR in a' ${TMP_DIR}/line_num.txt $TMP_DIR/body.vcf >${TMP_DIR}/body_filter.vcf
        ##stop here, extract snp by large than threshold by each chr and merge into one final vcf
        cat ${TMP_DIR}/header.vcf $TMP_DIR/body_filter.vcf >${dir}/${file}_CADD${version}_filtered${THRESHOLD}.vcf
    else 
        echo "use hg38 please!"
        #ln -s $VCF $TMP_DIR/${file}_nochr.vcf
    fi
fi

date
