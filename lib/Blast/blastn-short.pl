use strict;
use warnings;

use File::Basename;

my $output_file = $ARGV[0];
my $input_file    = $ARGV[1];

$output_file = "/scratch/cqs/shengq1/vickers/20150709_smallRNA_3018-KCV-35-38_human/identical_sequence_count_table/result/3018-KCV-35-38_sequence.blastn.tsv";
$input_file = "/scratch/cqs/shengq1/vickers/20150709_smallRNA_3018-KCV-35-38_human/identical_sequence_count_table/result/3018-KCV-35-38_sequence.count.fasta";

`blastn -task blastn-short -db nt -perc_identity 100 -remote -query $input_file -outfmt '6 qlen nident qacc sallacc salltitles' | awk '\$1 == \$2 {print}' | cut -f3- | uniq > $output_file`
