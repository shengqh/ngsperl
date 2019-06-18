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
		#cmd_discordants=($(ls -d ../align/${line}*/*[1-9].discordants.bam))
		#cmd_splitters=($(ls -d ../align/${line}*/*[1-9].splitters.bam))
#		mkdir -p /scratch/cqs/baiy7/Tim_proj/Family_WGS/FusorSV_ana/2_bamMerge_hg38/result/$line
		cat > $work_dir/$outdir/pbs/${line}_bam2cram.pbs <<EOL
#!/bin/bash
#SBATCH --mail-user=test
#SBATCH --mail-type=ALL
#SBATCH --nodes=$nodes
#SBATCH --ntasks=$ntasks
#SBATCH --time=48:00:00
#SBATCH --mem=40G
#SBATCH -o $work_dir/$outdir/log/${line}_bam2cram.log
ml GCC SAMtools
date
samtools view -S -C -T $genome -@ $ntasks -o $work_dir/$outdir/result/${line}/${line}.mkDup.sorted.recal.cram $work_dir/$outdir/result/${line}/${line}.mkDup.sorted.recal.bam
date
EOL

done
