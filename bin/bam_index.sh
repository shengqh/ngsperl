#!/bin/sh

for f in *.bam;
do 
  if [[ ! -s ${f}.bai ]]; then
    echo "samtools index $f"
    echo copy to local folder
    cp $f /workspace/shengq2/temp/
    localf="$(basename -- $f)"
    singularity exec -B /workspace -e /scratch/cqs_share/softwares/singularity/cqs-exomeseq.simg samtools index /workspace/shengq2/temp/$localf
    if [[ -s /workspace/shengq2/temp/${localf}.bai ]]; then
      echo copy index back
      cp /workspace/shengq2/temp/${localf}.bai .
      rm /workspace/shengq2/temp/${localf}.bai
    fi
    
    rm /workspace/shengq2/temp/$localf 
  fi
done
