use strict;
use Bio::SeqIO;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(extractFastaFile)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

sub extractFastaFile {
  my ( $inputFile, $pattern, $outputFile, $overwrite ) = @_;

  if ( !defined $overwrite ) {
    $overwrite = 1;
  }

  if ( -s $outputFile && !$overwrite ) {
    print STDOUT "$outputFile exists, ignored.\n";
    return;
  }

  print STDOUT "Extracting $outputFile ...\n";
  open( my $output, ">$outputFile" ) or die "Could not create file '$outputFile' $!";

  my $seqio = Bio::SeqIO->new( -file => $inputFile, -format => 'fasta' );
  while ( my $seq = $seqio->next_seq ) {
    my $id       = $seq->id;
    my $desc     = $seq->desc;
    my $sequence = $seq->seq;

    if ( $id =~ /$pattern/ ) {
      print $output ">$id $desc\n$sequence\n";
    }
  }

  close($output);
  print STDOUT "$outputFile extracted.\n";
}
