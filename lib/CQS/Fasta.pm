package CQS::Fasta;

use strict;
use warnings;

use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(get_sequence_names)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub get_sequence_names {
  my $fastafile = shift;

  my @names = ();
  open( IN, $fastafile ) or die "Couldn't optn $fastafile\n";
  while (<IN>) {
    if (/^>(\S+)/) {
      push( @names, $1 );
    }
  }
  close(IN);
  return @names;
}
