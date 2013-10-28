package CQS::StringUtils;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw( merge_string print_hash filter_array)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub merge_string {
  my ( $delimiter, @values ) = @_;
  my $first  = 1;
  my $result = "";
  foreach my $value (@values) {
    if ($first) {
      $result = $value;
      $first  = 0;
    }
    else {
      $result = $result . $delimiter . $value;
    }
  }
  return ($result);
}

sub print_hash {
  my $hash = shift;
  foreach my $k ( sort keys %{$hash} ) {
    print "$k => @{$hash->{$k}}\n";
  }
}

sub filter_array {
  my ($sourceFiles, $pattern ) = @_;

  if ( !defined $pattern ) {
    return $sourceFiles;
  }

  my @filteredFiles = ();
  for my $candidateFile ( @{$sourceFiles} ) {
    #print $candidateFile . "\n";
    if ( $candidateFile =~ m/$pattern/ ) {
      #print $candidateFile . " passed\n";
      push( @filteredFiles, $candidateFile );
    }
  }

  return \@filteredFiles;
}

1;
