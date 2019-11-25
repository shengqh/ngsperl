#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use Test::More tests => 2;
use Pipeline::PipelineUtils;

my $chromosomes = readChromosomeFromDictFile(dirname(__FILE__) . "/../../data/GRCm38.p6.genome.dict");
is (scalar(@$chromosomes), 22);
is (join(",", @$chromosomes), 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chrX,chrY,chrM');

1;
