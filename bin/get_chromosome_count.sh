#!/bin/sh

for f in *.bam;
do 
  if [[ ! -s ${f}.bai ]]; then
    echo "samtools index $f"
    singularity exec -e /scratch/cqs_share/softwares/singularity/cqs-exomeseq.simg samtools index $f 
  fi
  if [[ ! -s ${f}.chromosome.count ]]; then
    echo "samtools idxstats $f > ${f}.chromosome.count"
    singularity exec -e /scratch/cqs_share/softwares/singularity/cqs-exomeseq.simg samtools idxstats $f > ${f}.chromosome.count
  fi
done
