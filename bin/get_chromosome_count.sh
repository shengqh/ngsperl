#!/bin/sh

for f in *.bam ;
do 
  af="$(pwd)/$f"
  if [[ ! -e ${af}.bai ]]; then
    echo "samtools index $af"
    singularity exec -c -B /gpfs52/data:/data,/workspace -e /data/cqs/softwares/singularity/cqs-exomeseq.simg samtools index "$af" 
  fi
  if [[ ! -e ${af}.chromosome.count ]]; then
    echo "samtools idxstats $af > ${af}.chromosome.count"
    singularity exec -c -B /gpfs52/data:/data,/workspace -e /data/cqs/softwares/singularity/cqs-exomeseq.simg samtools idxstats "$af" > "${af}.chromosome.count"
  fi
done
