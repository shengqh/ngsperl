use strict;
use warnings;

my $output_file = $ARGV[0];
my $input_file  = $ARGV[1];

#$input_file = "../../data/3116_human_41_50.blastn.tsv";
#$output_file = "../../data/3116_human_41_50.blastn.table.tsv";

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

    if ( $genome =~ /human/i || $genome =~ /homo sapiens/i ) {
      $genome = "Human";
    }
    elsif ( $genome =~ /rat/i || $genome =~ /rattus norvegicus/i ) {
      $genome = "Rat";
    }
    elsif ( $genome =~ /mouse/i || $genome =~ /Mus musculus/i ) {
      $genome = "Mouse";
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

    #print $seq, "\t", $genome, "\n";
    $res{$genome}{$seq} = 1;
  }
  close($input);
}

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
      if ( !exists $sequenceMap{$seq} ) {
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

open( my $output, ">$output_file" ) or die "Cannot open $output_file";
printf $output "genome\tunique_sequence_count\tsequences\n";

my @sorted_genomes = sort {
  my $counta = scalar( @{ $merged{$a} } );
  my $countb = scalar( @{ $merged{$b} } );
  $countb <=> $counta
} keys %merged;

foreach my $name (@sorted_genomes) {
  my @sequences    = @{ $merged{$name} };
  my $count        = scalar(@sequences);
  my $sequence_str = join( ";", @sequences );
  printf $output "%s\t%d\t%s\n", $name, $count, $sequence_str;
}

close($output);
