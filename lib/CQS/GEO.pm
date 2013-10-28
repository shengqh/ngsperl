package CQS::GEO;

use strict;
use warnings;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Net::FTP;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(download_geo_all
      download_geo_supplementary
      download_geo_matrix
      download_geo_ids)
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

my $ftp_site = "ftp.ncbi.nlm.nih.gov";

sub uncompress_file {
  my ( $sourcefile, $targetfile ) = @_;

  if ( -e $sourcefile ) {
    if ( $sourcefile =~ /.gz$/i ) {
      if ( gunzip $sourcefile => $targetfile ) {
        unlink($sourcefile);
      }
      else {
        print "gunzip failed from $sourcefile => $targetfile\n: $GunzipError\n";
      }
    }
  }
}

sub download_geo_dir {
  my ( $dataset, $localdir, $ftpdir, $uncompress ) = @_;

  my $glob = '*.*';

  my $ftp = Net::FTP->new($ftp_site)
    or die "Could not connect to $ftp_site: $!";

  $ftp->login( "anonymous", '-anonymous@' )
    or die "Could not login into $ftp_site: $!";

  my $curftpdir = $ftpdir . $dataset . "/";

  if ( $ftp->cwd($curftpdir) ) {
    my @remote_files = $ftp->ls($glob);

    print "$curftpdir => $localdir\n";

    $ftp->binary();

    unless ( -d $localdir ) {
      mkdir $localdir or die "failed to create directory $localdir";
    }

    chdir $localdir or die "failed to change to directory $localdir";

    foreach my $file (@remote_files) {
      my $localfile = $file;
      my $uncompressedfile = substr $localfile, 0, length($localfile) - 3;

      if ( -e $localfile ) {
        if ($uncompress) {
          uncompress_file( $localfile, $uncompressedfile );
        }

        next;
      }

      if ($uncompress) {
        if ( -e $uncompressedfile ) {
          next;
        }
      }

      print "\tDownloading : $localdir/$file\n";
      $ftp->get( $file, $localfile . ".tmp" )
        or die "failed to download file $dataset/$file";

      rename($localfile . ".tmp", $localfile);
      
      if ($uncompress) {
        uncompress_file( $localfile, $uncompressedfile );
      }

      print "\tDownloaded : $localdir/$file\n";
    }
  }
  else {
    print "NO DATA : $curftpdir\n";
  }

  $ftp->quit();
}

sub download_geo_supplementary {
  my ( $dataset, $targetRootDir ) = @_;
  download_geo_dir( $dataset, $targetRootDir,
    "/pub/geo/DATA/supplementary/series/", 0 );
}

sub download_geo_matrix {
  my ( $dataset, $targetRootDir ) = @_;
  download_geo_dir( $dataset, $targetRootDir, "/pub/geo/DATA/SeriesMatrix/",
    1 );
}

sub download_geo_all {
  my ( $dataset, $targetRootDir ) = @_;
  download_geo_supplementary( @_, 0 );
  download_geo_matrix( @_, 1 );
}

sub download_geo_ids {
  my ($targetdir, $excludeddir, @datasets) = @_;

  foreach my $dataset (@datasets) {
    my $localdir       = $targetdir . $dataset;
    my $curexcludeddir = $excludeddir . $dataset;

    if ( !( -d $localdir ) && ( -d $curexcludeddir ) ) {
      next;
    }

    unless ( -d $localdir ) {
      mkdir $localdir or die "failed to create directory $localdir";
    }

    if (
      !has_file(
        $localdir,
        sub {
          my $file = shift;
          return $file =~ /\.cel$/i;
        }
      )
      )
    {
      download_geo_supplementary( $dataset, $localdir, 1 );
    }

    download_geo_matrix( $dataset, $localdir, 1 );
  }
}

1;

__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

CQS::GEO - Perl extension for download gene expression omnibus dataset from ftp server automatically

=head1 SYNOPSIS

  use CQS::GEO;
  
  my $dataset = 'GSE1561';

  my $targetdir = 'D:/temp/';

  download_supplementary($dataset, $targetdir);

  1;

=head1 DESCRIPTION

Stub documentation for CQS::Download, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

download_supplementary(dataset, targetDirectory) : 
  Download supplementary of geo dataset to target directory.

download_matrix(dataset, targetDirectory) : 
  Download matrix data of geo dataset to target directory.
  
=head1 SEE ALSO



=head1 AUTHOR

Quanhu SHENG<lt>quanhu.sheng@vanderbilt.edu / shengqh@gmail.com<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2012 by Quanhu SHENG

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.14.2 or,
at your option, any later version of Perl 5 you may have available.


=cut
