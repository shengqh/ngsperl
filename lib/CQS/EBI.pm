package CQS::EBI;

use strict;
use warnings;
use LWP::Simple;
use CQS::Download;
use CQS::HashUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(download_ebi_dataset)
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

my $filter = sub {
  my @url_array = @_;
  my @myresult  = ();

  foreach my $link (@url_array) {
    if ( $link =~ /(.+\.raw\.\d+\.zip)$/ ) {
      push( @myresult, "http://www.ebi.ac.uk" . $1 );
    }
    elsif ( $link =~ /(.+\.idf.txt)$/ ) {
      push( @myresult, "http://www.ebi.ac.uk" . $1 );
    }
    elsif ( $link =~ /(.+\.sdrf.txt)$/ ) {
      push( @myresult, "http://www.ebi.ac.uk" . $1 );
    }
  }

  return (@myresult);
};

my $makefilename = sub {
  my $oldurl = shift;
  $oldurl =~ /\/([^\/]+)$/;
  return ($1);
};

sub download_ebi_dataset {
  my ( $rootdir, $dataset ) = @_;
  my $url = "http://www.ebi.ac.uk/arrayexpress/experiments/$dataset";
  my $dir = "$rootdir/$dataset";

  if ( !( -d $dir ) ) {
    mkdir $dir;
  }

  download_files( $url, $dir, $filter, $makefilename );
}

1;

__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

CQS::GEO - Perl extension for download gene expression omnibus dataset from ftp server automatically

=head1 SYNOPSIS

  use CQS::EBI;

  download_ebi_dataset( "C:/Dataset", "E-MTAB-365" );

  1;

=head1 DESCRIPTION

Documentation for CQS::EBI

=head2 EXPORT

download_ebi_dataset(targetDirectory, dataset)
  
=head1 SEE ALSO



=head1 AUTHOR

Quanhu SHENG<lt>quanhu.sheng@vanderbilt.edu / shengqh@gmail.com<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2012 by Quanhu SHENG

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.14.2 or,
at your option, any later version of Perl 5 you may have available.


=cut
