#/usr/bin/bash
nodes=1
nodepertasks=8
ntasks=$[ $nodes*$nodepertasks ]
ntime="336:00:00"
#mem=$[ 10*($nodepertasks+1) ]
mem=180
singularity=0; # 1: ues singularity; 0 do not use singularity
work_dir="/gpfs23/scratch/cqs/baiy7/Tim_proj/Family_WGS/FusorSV_ana"
hg="hg38"

if [ $hg == "hg19" ]; then
genome="/scratch/cqs/shengq2/references/gatk/b37/bwa_index_0.7.17/human_g1k_v37.fasta"
indir="1_mapping"
outdir="2_bamMerge"
else
genome="/scratch/cqs/baiy7/Tim_proj/Family_WGS/genome/hg38/Homo_sapiens_assembly38.fasta"
indir="1_mapping_hg38"
outdir="2_bamMerge_hg38"
fi

ls -d $work_dir/$indir/result/align/*|cut -d / -f12|cut -d _ -f1|sort |uniq |while read line;
do
	echo $line
	#fullName=${line##*/}
	#leftName=${fullName%%_*}
	#rightName=${fullName#*_}
#	echo $leftName $rightName
#	if [ ! -d $work_dir/$outdir/result/$line ]; then	
#		mkdir -p $work_dir/$outdir/result/$line
#		if [ ! -d $work_dir/$outdir/pbs ]; then
#			mkdir -p $work_dir/$outdir/pbs	
#		fi
		cmd=($(ls -d $work_dir/$indir/result/align/${line}*/*[0-9].bam))
		#cmd_discordants=($(ls -d ../align/${line}*/*[1-9].discordants.bam))
		#cmd_splitters=($(ls -d ../align/${line}*/*[1-9].splitters.bam))
#		mkdir -p /scratch/cqs/baiy7/Tim_proj/Family_WGS/FusorSV_ana/2_bamMerge_hg38/result/$line
		cat > $work_dir/$outdir/pbs/${line}_merge.pbs <<EOL
#!/bin/bash
#SBATCH --mail-user=test
#SBATCH --mail-type=ALL
#SBATCH --nodes=$nodes
#SBATCH --ntasks=$ntasks
#SBATCH --time=$ntime
#SBATCH --mem=${mem}G
#SBATCH -o $work_dir/$outdir/log/${line}_merge.log
source ~/path.txt	
localdir=/tmp/baiy7_\${SLURM_JOBID}
mkdir \${localdir}
date
samtools merge -u -r - ${cmd[0]} ${cmd[1]} ${cmd[2]} ${cmd[3]} | \\
sambamba sort -u --sort-picard -t $ntasks -m ${mem}G -o \${localdir}/${line}.sorted.queryName.bam /dev/stdin

date
java -Xmx$((${mem}-4))g -jar /scratch/cqs/shengq2/local/bin/picard/picard.jar MarkDuplicates ASSUME_SORT_ORDER="queryname" COMPRESSION_LEVEL=0 I=\${localdir}/${line}.sorted.queryName.bam O=\${localdir}/${line}.mkDup.bam M=\${localdir}/${line}.mkDup.metrics.txt 
rm \${localdir}/${line}.sorted.queryName.bam
date
sambamba sort -t $ntasks -m ${mem}G -o $work_dir/$outdir/result/$line/${line}.mkDup.sorted.bam \${localdir}/${line}.mkDup.bam
cp \${localdir}/${line}.mkDup.metrics.txt $work_dir/$outdir/result/$line/
date
java -Xmx40g -jar /scratch/cqs/softwares/gatk3.jar -T BaseRecalibrator -nct 4 \\
-R $genome \\
-I $work_dir/$outdir/result/$line/${line}.mkDup.sorted.bam \\
-o $work_dir/$outdir/result/$line/${line}.table	\\
-knownSites /scratch/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz \\
-knownSites /scratch/cqs/references/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \\
-knownSites /scratch/cqs/references/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz \\
--preserve_qscores_less_than 6 -dfrac .1 \\
-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22
date
java -Xmx40g -jar /scratch/cqs/softwares/gatk3.jar -T PrintReads -nct 4 \\
-I $work_dir/$outdir/result/$line/${line}.mkDup.sorted.bam \\
-o $work_dir/$outdir/result/$line/${line}.mkDup.sorted.recal.bam \\
-BQSR $work_dir/$outdir/result/$line/${line}.table \\
-R /scratch/cqs/baiy7/Tim_proj/Family_WGS/genome/hg38/Homo_sapiens_assembly38.fasta \\
-preserveQ 6 \\
-SQQ 10 -SQQ 20 -SQQ 30 \\
--disable_indel_quals

date
#rm -r \${localdir}
#rm $work_dir/$outdir/result/$line/${line}.bam $work_dir/$outdir/result/$line/${line}.sorted.queryName.bam $work_dir/$outdir/result/$line/${line}.mkDup.bam
EOL



#	fi
done
