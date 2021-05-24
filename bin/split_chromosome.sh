#!/bin/bash

mkdir chromosomes
cd chromosomes
csplit -s -z ../$1 '/>/' '{*}'
for i in xx* ; do \
  n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ; \
  mv "$i" "$n.fa" ; \
done
