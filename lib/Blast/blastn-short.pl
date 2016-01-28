use strict;
use warnings;

use File::Basename;
use Bio::SeqIO;

my $output_file = $ARGV[0];
my $input_file  = $ARGV[1];

$output_file = "/scratch/cqs/shengq1/vickers/20150709_smallRNA_3018-KCV-35-38_human/identical_sequence_count_table/result/3018-KCV-35-38_sequence.blastn.tsv";
$input_file  = "/scratch/cqs/shengq1/vickers/20150709_smallRNA_3018-KCV-35-38_human/identical_sequence_count_table/result/3018-KCV-35-38_sequence.count.fasta";

my $seqio = Bio::SeqIO->new( -file => $input_file, '-format' => 'Fasta' );
if ( -e $output_file ) { unlink $output_file; }

while ( my $seq = $seqio->next_seq ) {
  my $fa_name = $seq->id . '.fasta';
  my $seqio_obj = Bio::SeqIO->new( -file => '>$fa_name', -format => 'fasta' );
  $seqio_obj->write_seq($seq);
  print $seq->id . "\n";
  #`blastn -task blastn-short -db nt -perc_identity 100 -remote -query $fa_name -outfmt '6 qlen nident qacc sallacc salltitles' | awk '\$1 == \$2 {print}' | cut -f3- | sort | uniq >> $output_file`;
}

