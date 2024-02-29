VAR="/panfs/accrepfs.vampire/data/cqs/ramirema/ken_lau/20240221_11033/11033-DA-1_S1_L005_R1_001.fastq.gz,/panfs/accrepfs.vampire/data/cqs/ramirema/ken_lau/20240221_11033/11033-DA-1_S1_L005_R2_001.fastq.gz"

IFS=, read -r var1 var2 <<< $VAR
echo $var1
echo $var2
DIR="$(dirname "${var1}")"
FILE="$(basename "${var1}")"
FILE_PREFIX="$(echo $FILE | sed -e "s/_S.*$//")"
echo $DIR
echo $FILE_PREFIX

