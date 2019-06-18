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
excludeBed="/scratch/cqs/baiy7/tools/ceph18.b37.lumpy.exclude.2014-01-15.bed"
indir="2_bamMerge"
outdir="3_smoove_call"
else
genome="/scratch/cqs/baiy7/Tim_proj/Family_WGS/genome/hg38/Homo_sapiens_assembly38.fasta"
excludeBed="/scratch/cqs/baiy7/tools/exclude.cnvnator_100bp.GRCh38.20170403.bed"
indir="2_bamMerge_hg38"
outdir="3_smoove_call_hg38"
fi

ls -d $work_dir/$indir/result/*/ |cut -d / -f11|while read line;
do
	echo $line
	#fullName=${line##*/}
	#leftName=${fullName%%_*}
	#rightName=${fullName#*_}
#	echo $leftName $rightName
	if [ ! -d $work_dir/$outdir/result ]; then	
		mkdir -p $work_dir/$outdir/result/
	fi
        if [ ! -d $work_dir/$outdir/pbs ]; then
                mkdir -p $work_dir/$outdir/pbs
        fi
        if [ ! -d $work_dir/$outdir/log ]; then
                mkdir -p $work_dir/$outdir/log
        fi
		#cmd_discordants=($(ls -d ../align/${line}*/*[1-9].discordants.bam))
		#cmd_splitters=($(ls -d ../align/${line}*/*[1-9].splitters.bam))
#		mkdir -p /scratch/cqs/baiy7/Tim_proj/Family_WGS/FusorSV_ana/2_bamMerge_hg38/result/$line
		cat > $work_dir/$outdir/pbs/${line}_smoove_call.pbs <<EOL
#!/bin/bash
#SBATCH --mail-user=test
#SBATCH --mail-type=ALL
#SBATCH --nodes=$nodes
#SBATCH --ntasks=$ntasks
#SBATCH --time=48:00:00
#SBATCH --mem=40G
#SBATCH -o $work_dir/$outdir/log/${line}_smoove_call.log
date
if [ ! -d $work_dir/$outdir/result/call ]; then
        mkdir -p $work_dir/$outdir/result/call
fi
singularity exec /scratch/cqs/baiy7/tools/smoove.simg smoove call --outdir $work_dir/$outdir/result/call/ --exclude $excludeBed --name $line --fasta $genome -p 3 --excludechroms "hs37d5,~:,~^GL,~decoy,~alt" --genotype $work_dir/$indir/result/${line}/${line}.mkDup.sorted.recal.bam
date
EOL

                cat > $work_dir/$outdir/pbs/${line}_smoove_genotype.pbs <<EOL
#!/bin/bash
#SBATCH --mail-user=test
#SBATCH --mail-type=ALL
#SBATCH --nodes=$nodes
#SBATCH --ntasks=$ntasks
#SBATCH --time=72:00:00
#SBATCH --mem=80G
#SBATCH -o $work_dir/$outdir/log/${line}_smoove_genotype.log
date
if [ ! -d $work_dir/$outdir/result/genotype ]; then
        mkdir -p $work_dir/$outdir/result/genotype
fi
singularity exec /scratch/cqs/baiy7/tools/smoove.simg smoove genotype --outdir $work_dir/$outdir/result/genotype/ --name ${line}-joint --fasta $genome -d -x -p 7 --vcf $work_dir/$outdir/result/merge/merged.sites.vcf.gz $work_dir/$indir/result/${line}/${line}.mkDup.sorted.recal.bam
date
EOL

done
## smoove_merge pbs
cat > $work_dir/$outdir/pbs/smoove_merge.pbs <<EOL
#!/bin/bash
#SBATCH --mail-user=test
#SBATCH --mail-type=ALL
#SBATCH --nodes=$nodes
#SBATCH --ntasks=$ntasks
#SBATCH --time=48:00:00
#SBATCH --mem=40G
#SBATCH -o $work_dir/$outdir/log/smoove_merge.log
date
if [ ! -d $work_dir/$outdir/result/merge ]; then
        mkdir -p $work_dir/$outdir/result/merge
fi
singularity exec /scratch/cqs/baiy7/tools/smoove.simg smoove merge --outdir $work_dir/$outdir/result/merge/ --name merged --fasta $genome $work_dir/$outdir/result/call/*.genotyped.vcf.gz
date
EOL
