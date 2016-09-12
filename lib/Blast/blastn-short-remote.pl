use strict;
use warnings;
use threads;
use POSIX;

use File::Basename;
use Bio::SeqIO;
use Getopt::Long;

my $usage = "

Synopsis:

perl blast-short-thread.pl -i fasta_file -o output_file -t number_of_thread

Options:
  -i|--input {file}       Input sequence file in fasta format
  -o|--output {file}      Output file
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
  'i|input=s'  => \$input_file,
  'o|output=s' => \$output_file,
  't|thread=s' => \$num_of_threads,
  'h|help'     => \$help,
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

die "$output_file already exists" if ( -e $output_file );

if ( !defined $num_of_threads ) {
  $num_of_threads = 1;
}

my $seqio = Bio::SeqIO->new( -file => $input_file, '-format' => 'Fasta' );
my @sequences = ();
while ( my $seq = $seqio->next_seq ) {
  push( @sequences, $seq );
}

my $total_count = scalar(@sequences);
if ( $total_count < $num_of_threads ) {
  $num_of_threads = $total_count;
}
my $step = ceil( $total_count / $num_of_threads );

# use the initThreads subroutine to create an array of threads.
my @threads = initThreads();

# Loop through the array:
my $from  = 0;
my @files = ();
foreach (@threads) {

  # Tell each thread to perform our 'doOperation()' subroutine.
  my $curoutfile = $output_file . $from;
  $_ = threads->create( \&doOperation, $from, $curoutfile );
  $from = $from + $step;

  push( @files, $curoutfile );
}

# This tells the main program to keep running until all threads have finished.
foreach (@threads) {
  $_->join();
}

foreach my $result_file (@files) {
  `cat $result_file >> $output_file`;
}

print "Program Done!\n";

####################### SUBROUTINES ############################
sub initThreads {
  my @initThreads;
  for ( my $i = 1 ; $i <= $num_of_threads ; $i++ ) {
    push( @initThreads, $i );
  }
  return @initThreads;
}

sub doOperation {
  my ( $from, $cur_output_file ) = @_;
  my $file = shift;
  my $to   = $from + $step - 1;
  if ( $to >= $total_count ) {
    $to = $total_count - 1;
  }

  # Get the thread id. Allows each thread to be identified.
  my $id = threads->tid();
  print "Thread $id start!\n";

  for ( my $current_index = $from ; $current_index <= $to ; $current_index++ ) {
    my $seq       = $sequences[$current_index];
    my $fa_name   = $seq->id . '.fasta';
    my $fa_output = $seq->id . '.fasta.output';
    my $seqio_obj = Bio::SeqIO->new( -file => ">$fa_name", -format => 'fasta' );
    $seqio_obj->write_seq($seq);
    my $datastring = localtime();
    
    print $datastring . " : thread " . $id . " : " . ( $current_index + 1 ) . "/" . ( $to + 1 ) . " : " . $seq->id . "\n";
    `blastn -task blastn-short -db nt -perc_identity 100 -remote -query $fa_name -outfmt '6 qlen nident qacc sallacc salltitles' | awk '\$1 == \$2 {print}' | cut -f3- | sort | uniq >> $fa_output`;
    if ( !-e $fa_output ) {
      print STDERR "blastn failed for $fa_name ";
    }
    else {
      `cat $fa_output >> $cur_output_file`;
      unlink($fa_output);
    }

    unlink($fa_name);
  }

  print "Thread $id done!\n";

  # Exit the thread
  threads->exit();
}

