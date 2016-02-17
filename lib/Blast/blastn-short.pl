use strict;
use warnings;

use File::Basename;
use Bio::SeqIO;

my $output_file = $ARGV[0];
my $input_file  = $ARGV[1];
my $from        = $ARGV[2];
my $to          = $ARGV[3];

#$output_file = "/scratch/cqs/shengq1/temp/3018-KCV-35-38_sequence.blastn.tsv";
#$input_file  = "/scratch/cqs/shengq1/vickers/20150709_smallRNA_3018-KCV-35-38_human/identical_sequence_count_table/result/3018-KCV-35-38_sequence.count.fasta";

my $seqio = Bio::SeqIO->new( -file => $input_file, '-format' => 'Fasta' );
if ( -e $output_file ) { unlink $output_file; }

my @sequences = ();
while ( my $seq = $seqio->next_seq ) {
  push( @sequences, $seq );
}

my $total_count = scalar(@sequences);
if ( !defined $from ) {
  $from = 0;
  $to   = $total_count - 1;
}
elsif ( !defined $to ) {
  $to = $total_count - 1;
}

$from = $from < 0 ? 0 : $from;
$to = $to >= $total_count ? $total_count - 1 : $to;

for ( my $current_index = $from ; $current_index <= $to ; $current_index++ ) {
  my $seq       = $sequences[$current_index];
  my $fa_name   = $seq->id . '.fasta';
  my $fa_output = $seq->id . '.fasta.output';
  my $seqio_obj = Bio::SeqIO->new( -file => ">$fa_name", -format => 'fasta' );
  $seqio_obj->write_seq($seq);
  my $datastring = localtime();
  print $datastring . " : " . $current_index . "/" . $to . " : " . $seq->id . "\n";
  `blastn -task blastn-short -db nt -perc_identity 100 -remote -query $fa_name -outfmt '6 qlen nident qacc sallacc salltitles' | awk '\$1 == \$2 {print}' | cut -f3- | sort | uniq >> $fa_output`;

  if ( !-e $fa_output ) {
    print STDERR "blastn failed for $fa_name ";
  }
  else {
    `cat $fa_output >> $output_file`;
    unlink($fa_output);
  }

  unlink($fa_name);
}

