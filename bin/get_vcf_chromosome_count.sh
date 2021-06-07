#!/bin/sh

for f in *.indel.pass.vcf.gz;
do 
  if [[ ! -s ${f}.chromosome.count ]]; then
    echo "$f"
    zcat $f | grep -v "^#" | cut -f1 | uniq -c | awk '{print $2"\t"$1}' > ${f}.chromosome.count
  fi
done
