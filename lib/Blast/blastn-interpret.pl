use strict;
use warnings;

use Getopt::Long;

my $usage = "

Synopsis:

perl blast-interpret.pl -i blastn_result_file -o interpret_result_file

Options:
  -i|--input {file}       Input blastn result file
  -o|--output {file}      Output table file
  -h|--help               This page.
";

Getopt::Long::Configure('bundling');

my $input_file;
my $output_file;
my $help;

GetOptions(
  'i|input=s'  => \$input_file,
  'o|output=s' => \$output_file,
  'h|help'     => \$help,
);

if ( defined $help ) {
  print $usage;
  exit(1);
}

if ( !defined $input_file || !-e $input_file ) {
  die "Input valid input file!";
}

if ( !defined $output_file ) {
  die "Define output file!";
}

#die "$output_file already exists" if ( -e $output_file );

my %res;
foreach my $file ( split( ",", $input_file ) ) {
  open( my $input, "<$file" ) or die "Cannot open $file";
  while (<$input>) {
    chomp;
    my @parts  = split( "\t", $_ );
    my $seq    = $parts[0];
    my $genome = $parts[2];
    $genome =~ s/PREDICTED: //g;
    $genome =~ s/Homologies in //g;

    if ( $genome =~ /human[ \t]/i || $genome =~ /homo sapiens/i ) {
      $genome = "Human";
    }
    elsif ( $genome =~ /rat[ \t]/i || $genome =~ /rattus norvegicus/i ) {
      $genome = "Rat";
    }
    elsif ( $genome =~ /mouse[ \t]/i || $genome =~ /Mus musculus/i ) {
      $genome = "Mouse";
    }
    elsif ( $genome =~ /zebrafish[ \t]/i || $genome =~ /Danio rerio/i ) {
      $genome = "Zebrafish";
    }
    else {
      my @tokens = split( " ", $genome );
      if ( scalar(@tokens) == 1 ) {
        if ( $tokens[0] ne "N/A" ) {
          $genome = $tokens[0];
        }
        else {
          next;
        }
      }
      else {
        if ( $tokens[1] =~ /^\d/ ) {
          $genome = $tokens[0];
        }
        else {
          $genome = $tokens[0] . " " . $tokens[1];
        }
      }
    }

    $res{$genome}{$seq} = 1;
  }
  close($input);
}

#merge genome whose mapped reads are contained in another genome
my %merged;
my $index = 0;
while (1) {
  my @sorted_genomes = sort {
    my $counta = keys %{ $res{$a} };
    my $countb = keys %{ $res{$b} };
    $countb <=> $counta
  } keys %res;

  my $genome_count = scalar(@sorted_genomes);
  if ( $index >= $genome_count ) {
    last;
  }

  my $name        = $sorted_genomes[$index];
  my %sequenceMap = %{ $res{$name} };
  my @sequences   = sort keys %sequenceMap;
  my $count       = scalar(@sequences);

  for ( my $another = $genome_count - 1 ; $another > $index ; $another-- ) {
    my $another_name      = $sorted_genomes[$another];
    my @another_sequences = sort keys %{ $res{$another_name} };
    my $another_count     = scalar(@another_sequences);

    my $issubset = 1;
    foreach my $seq (@another_sequences) {
      if ( not exists $sequenceMap{$seq} ) {
        $issubset = 0;
        last;
      }
    }

    if ($issubset) {
      if ( $count == $another_count ) {
        $name = $name . ";" . $another_name;
      }
      delete $res{$another_name};
    }
  }

  $merged{$name} = \@sequences;
  $index++;
}

sub getSequenceGenomeMap {
  my $merged = shift;
  my $result = {};
  foreach my $name ( keys %$merged ) {
    my @sequences = @{ $merged->{$name} };
    foreach my $seq (@sequences) {
      if ( not exists $result->{$seq} ) {
        $result->{$seq} = {};
      }
      $result->{$seq}->{$name} = 1;
    }
  }

  return $result;
}

#remove genome without unique sequence
while (1) {

  #build sequence/genome map
  my $sequenceGenomeMap = getSequenceGenomeMap( \%merged );

  #sort genomes by increasing of count
  my @sorted_genomes = sort {
    my $counta = scalar( @{ $merged{$a} } );
    my $countb = scalar( @{ $merged{$b} } );
    $counta <=> $countb;
  } keys %merged;

  #remove first genome without unique sequence
  my $removed = 0;
  foreach my $genome (@sorted_genomes) {
    my $unique    = 0;
    my @sequences = @{$merged{$genome}};
    foreach my $sequence (@sequences) {
      if ( scalar( keys %{ $sequenceGenomeMap->{$sequence} } ) == 1 ) {
        $unique = 1;
        last;
      }
    }

    if ( not $unique ) {
      $removed = 1;
      delete $merged{$genome};
      last;
    }
  }

  if ( not $removed ) {
    last;
  }
}

my $sequenceGenomeMap = getSequenceGenomeMap( \%merged );
my @sorted_genomes    = sort {
  my $counta = scalar( @{ $merged{$a} } );
  my $countb = scalar( @{ $merged{$b} } );
  $countb <=> $counta;
} keys %merged;

my $sequenceFile = $output_file;
$sequenceFile =~ s{\.[^.]*$}{.sequences.tsv};
open( my $output, ">$sequenceFile" ) or die "Cannot open $sequenceFile";
printf $output "sequence\t" . join( "\t", @sorted_genomes ) . "\n";

my @sorted_sequences = sort keys %$sequenceGenomeMap;
if(index($sorted_sequences[0], "_")){
  @sorted_sequences = sort {
    my ($counta) = $a =~ /_(\d+)$/;
    my ($countb) = $b =~ /_(\d+)$/;
    $countb <=> $counta;
  } @sorted_sequences;
};

foreach my $seq (@sorted_sequences ) {
  printf $output $seq;
  my %genomes = %{ $sequenceGenomeMap->{$seq} };
  foreach my $genome (@sorted_genomes) {
    if ( exists $genomes{$genome} ) {
      printf $output "\t1";
    }
    else {
      printf $output "\t0";
    }
  }
  printf $output "\n";
}

close($output);

open( $output, ">$output_file" ) or die "Cannot open $output_file";
printf $output "genome\tsequence_count\tunique_sequence_count\n";

foreach my $name (@sorted_genomes) {
  my @sequences           = @{ $merged{$name} };
  my $uniqueSequenceCount = 0;
  foreach my $seq (@sequences) {
    my $genomes = $sequenceGenomeMap->{$seq};
    if ( scalar( keys %$genomes ) == 1 ) {
      $uniqueSequenceCount++;
    }
  }

  my $count = scalar(@sequences);
  my $sequence_str = join( ";", @sequences );
  printf $output "%s\t%d\t%d\n", $name, $count, $uniqueSequenceCount;
}

close($output);
