#!/usr/bin/perl
package CQS::RmdUtils;

use strict;
use warnings;
use File::Basename;
use File::Copy;
use String::Util ':all';

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw( build_rmd_file )] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub build_rmd_file {
  my ($config, $section, $result_dir, $result_file) = @_;

  my $rtemplate = $config->{$section}{report_rmd_file};
  if ( ! -e $rtemplate ) {
    $rtemplate = dirname(__FILE__) . "/" . $rtemplate;
  }
  my $rfile     = $result_dir . "/" . $result_file;

  if ( defined $config->{$section}{additional_rmd_files} ) {
    my @additional_rmd_files = split( ',|;', $config->{$section}{additional_rmd_files} );
    for my $rmd_file (@additional_rmd_files) {
      if (! -e $rmd_file) {
        $rmd_file=dirname(__FILE__) . "/" . $rmd_file;
      }
      my $additional_rmd_files_destination     = $result_dir . "/" . basename($rmd_file);
      unlink $additional_rmd_files_destination if -e $additional_rmd_files_destination;
      #print("copy " . $rmd_file . " to " . $additional_rmd_files_destination . "\n");
      copy( $rmd_file, $additional_rmd_files_destination ) or die "Copy " . $rmd_file . " failed: $! ";
    }
  }

  copy( $rtemplate, $rfile ) or die "Copy failed $rtemplate: $!";

  return($rfile);
}

1;
