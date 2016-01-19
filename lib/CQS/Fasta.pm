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
  open( my $in, $fastafile ) or die "Couldn't open $fastafile\n";
  while (<$in>) {
    if (/^>(\S+)/) {
      push( @names, $1 );
    }
  }
  close($in);
  return @names;
}
