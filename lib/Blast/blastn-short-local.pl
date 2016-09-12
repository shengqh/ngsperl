use strict;
use warnings;
use threads;
use POSIX;

use File::Basename;
use Bio::SeqIO;
use Getopt::Long;

my $usage = "

Synopsis:

perl blast-short-local.pl -i fasta_file -o output_file -l localdb_folder -t number_of_thread

Options:
  -i|--input {file}       Input sequence file in fasta format
  -o|--output {file}      Output file
  -l|--localdb            Local blast database (if not defined, using remote server)
  -t|--thread {integer}   Number of thread
  -h|--help               This page.
";

Getopt::Long::Configure('bundling');

my $input_file;
my $output_file;
my $num_of_threads;
my $local_db;
my $help;

GetOptions(
  'i|input=s'   => \$input_file,
  'o|output=s'  => \$output_file,
  'l|localdb=s' => \$local_db,
  't|thread=s'  => \$num_of_threads,
  'h|help'      => \$help,
);
if ( defined $help ) {
  print $usage;
  exit(1);
}

#$output_file = "/scratch/cqs/shengq1/temp/3018-KCV-35-38_sequence.blastn.tsv";
#$input_file  = "/scratch/cqs/shengq1/vickers/20150709_smallRNA_3018-KCV-35-38_human/identical_sequence_count_table/result/3018-KCV-35-38_sequence.count.fasta";

if ( !defined $input_file || !-e $input_file ) {
  die "Input valid input file!";
}

if ( !defined $output_file ) {
  die "Define output file!";
}

if ( !defined $local_db ) {
  die "Define folder contains local blast database!";
}

die "$output_file already exists" if ( -e $output_file );

my $thread_option = "";
if ( defined $num_of_threads ) {
  $thread_option = "-num_threads $num_of_threads";
}

my $seqio = Bio::SeqIO->new( -file => $input_file, '-format' => 'Fasta' );
my @sequences = ();
while ( my $seq = $seqio->next_seq ) {
  push( @sequences, $seq );
}

my $total_sequence = scalar(@sequences);
for ( my $current_index = 0 ; $current_index < $total_sequence ; $current_index++ ) {
  my $seq       = $sequences[$current_index];
  my $fa_name   = $seq->id . '.fasta';
  my $fa_output = $seq->id . '.fasta.output';
  my $seqio_obj = Bio::SeqIO->new( -file => ">$fa_name", -format => 'fasta' );
  $seqio_obj->write_seq($seq);
  my $datastring = localtime();

  print $datastring . " : " . ( $current_index + 1 ) . "/" . $total_sequence . " : " . $seq->id . "\n";
`blastn -task blastn-short -db $local_db/nt -perc_identity 100 -query $fa_name $thread_option -outfmt '6 qlen nident qacc sallacc salltitles' | awk '\$1 == \$2 {print}' | cut -f3- | sort | uniq >> $fa_output`;
  if ( !-e $fa_output ) {
    print STDERR "blastn failed for $fa_name ";
  }
  else {
    `cat $fa_output >> $output_file`;
    unlink($fa_output);
  }

  unlink($fa_name);
}

print "Program Done!\n";
