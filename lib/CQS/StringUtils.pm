package CQS::StringUtils;

use strict;
use warnings;
use Algorithm::Loops qw( NestedLoops );

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw( capture_regex_groups merge_string print_hash filter_array string_combination string_repeat)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub capture_regex_groups {
  my ( $source, $pattern ) = @_;
  $source =~ /$pattern/; # matches
  my $result = "";
  foreach my $exp (1..$#-) {
    no strict 'refs';
    $result = $result . $$exp;
  }

  return($result);
}

sub string_repeat {
  my ( $array, $repeatTime ) = @_;
  my $result = [];
  for my $str (@$array) {
    push @$result, ($str) x $repeatTime;
  }
  return $result;
}

sub string_combination {
  my ( $arrayOfStrArray, $delimiter ) = @_;
  my $result = [];

  my @b;
  NestedLoops( $arrayOfStrArray, sub { push @b, [@_] } );

  for my $strArray (@b) {
    my @notEmpty = grep { length($_) > 0 } @$strArray;
    push( @$result, join( $delimiter, @notEmpty ) );
  }
  return $result;
}

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
  my ( $sourceFiles, $pattern, $canReturnEmpty ) = @_;

  if ( !defined $pattern || $pattern eq "" ) {
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

  if ( $pattern !~ /version\$/ ) {
    if ( scalar(@filteredFiles) == 0 && !$canReturnEmpty ) {
      my $str = join( ', ', @{$sourceFiles} );
      warn("No file in array [$str] matched with pattern $pattern, empty array returned.");
    }
  }
  return \@filteredFiles;
}

1;
