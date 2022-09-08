#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 2;
use SRA::SRAUtils;

is_deeply(SrrToUrl('SRR5029211'), ['ftp.sra.ebi.ac.uk/vol1/fastq/SRR502/001/SRR5029211/SRR5029211.fastq.gz']);
is_deeply(SrrToUrl('SRR5127903'), ['ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/003/SRR5127903/SRR5127903_1.fastq.gz', 'ftp.sra.ebi.ac.uk/vol1/fastq/SRR512/003/SRR5127903/SRR5127903_2.fastq.gz']);

1;
