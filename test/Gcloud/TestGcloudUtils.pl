#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 1;
use Gcloud::GcloudUtils;

is( replaceSampleNameBWACommand("bwa mem -K 100000000 -p -v 3 -t 16 -R '\@RG\\tID:TP0097_Norecur_Base\\tPU:illumina\\tLB:TP0097_Norecur_Base\\tSM:TP0097_Norecur_Base\\tPL:illumina' -Y \$bash_ref_fasta", "TEST"), 
"bwa mem -K 100000000 -p -v 3 -t 16 -R '\@RG\\tID:TEST\\tPU:illumina\\tLB:TEST\\tSM:TEST\\tPL:illumina' -Y \$bash_ref_fasta");

1;
