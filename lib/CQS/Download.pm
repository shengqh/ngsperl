package CQS::Download;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(download_url_and_parse_href
      download_files
      download_dirs)
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use LWP::UserAgent;
use LWP::Simple;
use HTML::Parser;
use File::Temp;

sub download_url_and_parse_href {
  my ( $url, $urlfilter ) = @_;

  my $content = get($url);

  if ( defined($content) ) {
    my @url_array;
    my $parser = HTML::Parser->new(
      api_version => 3,
      start_h     => [
        sub {
          my ( $tag, $attr ) = @_;
          if ( ( $tag =~ /^a$/ ) and ( defined $attr->{'href'} ) ) {
            push( @url_array, $attr->{'href'} );
          }
        },
        "tagname, attr"
      ],
    );

    $parser->parse($content);

    if ($urlfilter) {
      @url_array = $urlfilter->(@url_array);
    }

    return ( 1, \@url_array );
  }

  return ( 0, 0 );
}

#give url and directory, download files to directory
sub download_files {
  my ( $rooturl, $rootdir, $urlfilter, $makefilename ) = @_;

  #print "In download_files : root url = $rooturl \n";
  #print "In download_files : root dir = $rootdir \n";

  unless ( -d $rootdir ) {
    mkdir $rootdir or die "failed to create directory $rootdir";
  }

  if ( !( $rootdir =~ /\/$/ ) ) {
    $rootdir = $rootdir . "/";
  }

  my ( $download_ok, $url_arrayref ) =
    download_url_and_parse_href( $rooturl, $urlfilter );

  $download_ok or die "failed to download and parse tags from $rooturl";

  my @tmp_url_array = @$url_arrayref;

  foreach (@tmp_url_array) {
    my $tmp_link = $_;

    #start with ?, ignore
    if ( $tmp_link =~ /^\?/ ) {
      next;
    }

    #end with /, is Directory, ignore
    if ( $tmp_link =~ /\/$/ ) {
      next;
    }

    #construct local file name
    my $filename;
    if ($makefilename) {
      $filename = $makefilename->($tmp_link);
      $filename = $rootdir . $filename;
    }
    else {
      $filename = $rootdir . $tmp_link;
    }

    #if file exists, ignore
    if ( -e $filename ) {
      next;
    }

    #construct file url
    my $tmpurl;
    if ( $tmp_link =~ /^http/ ) {
      $tmpurl = $tmp_link;
    }
    else {
      $tmpurl = $rooturl . $tmp_link;
    }

    #download file
    print "\n DOWNLOADING: $tmpurl";

    my $tmpFile = File::Temp->new( DIR => $rootdir )->filename;
    if ( -s $tmpFile ) {
      unlink($tmpFile);
    }

    getstore( $tmpurl, $tmpFile );

    if ( -s $tmpFile ) {
      rename( $tmpFile, $filename ) || die ( "Error in renaming $tmpFile to $filename" );
    }

    if ( -s $filename ) {
      print "\n INFO: $filename is downloaded";
    }
    else {
      print "\n ERROR: $filename is NOT downloaded";
    }

    sleep 1;
  }
}

#give url and directory, download files in subdirectories to target directory
sub download_dirs {
  my ( $rooturl, $rootdir, $dirfilter, $filefilter ) = @_;

  #print "In download_dirs : root url = $rooturl \n";
  #print "In download_dirs : root dir = $rootdir \n";

  unless ( -d $rootdir ) {
    mkdir $rootdir or die "failed to create directory $rootdir";
  }

  my ( $download_ok, $url_arrayref ) =
    download_url_and_parse_href( $rooturl, $dirfilter );

  $download_ok or die "failed to download and parse tags from $rooturl";

  my @tmp_url_array = @$url_arrayref;

  foreach (@tmp_url_array) {
    my $tmp_link = $_;

    #start with ?, ignore
    if ( $tmp_link =~ /^\?/ ) {
      next;
    }

    #end with '/' and without internal '/', is subdirectory
    if ( $tmp_link =~ /^[^\/]*\/$/ ) {
      my $cururl = $rooturl . $tmp_link;
      my $curdir = $rootdir . $tmp_link;

      download_files( $cururl, $curdir, $filefilter );
    }
  }
}

1;

__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

CQS::Download - Perl extension for batch download file from directory or subdirectory through http or https

=head1 SYNOPSIS

  use CQS::Download;
  
  my $microarrayurl = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/brca/cgcc/unc.edu/agilentg4502a_07_3/transcriptome/';
  my $microarraydir = 'D:/data/microarray/';

  my $microarrayfilter = sub{
    my @url_array = @_;
  
    my @myresult=();
  
    foreach my $link(@url_array)
    {
      if (!($link =~ /\/$/)){
        next;
      }
    
      if ($link =~ /Level_[12]/) 
      {
        next;
      }
    
      push(@myresult, $link);
    }
  
    return (@myresult);  
  };

  download_dirs($microarrayurl, $microarraydir, $microarrayfilter);

  1;

=head1 DESCRIPTION

Stub documentation for CQS::Download, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

download_url(url) : 
  Download source code from assigned url.

download_url_and_parse_href(url, urlfilter) : 
  Download source code from assigned url and filter sub-url based on urlfilter if urlfilter is assigned.
  
download_files(rooturl, rootdir, urlfilter) :
  Download source code from rooturl, parse file url, filter file url based on urlfilter if urlfilter is assigned,
  download file and save to rootdir.
  
download_dirs(rooturl, rootdir, dirfilter, filefilter) :
  Download source code from rooturl, parse sub-dir url, filter sub-dir url based on dirfilter if dirfilter is assigned,
  download source code from sub-dir url, parse file url, filter file url based on urlfilter if urlfilter is assigned,
  download file and save to rootdir.
  
=head1 SEE ALSO



=head1 AUTHOR

Quanhu SHENG<lt>quanhu.sheng@vanderbilt.edu / shengqh@gmail.com<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2012 by Quanhu SHENG

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.14.2 or,
at your option, any later version of Perl 5 you may have available.


=cut
