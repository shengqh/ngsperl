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

  my $rtemplate = dirname(__FILE__) . "/" . $config->{$section}{report_rmd_file};
  my $rfile     = $result_dir . "/" . $result_file;

  if ( defined $config->{$section}{additional_rmd_files} ) {
    my @additional_rmd_files = split( ';', $config->{$section}{additional_rmd_files} );
    for my $rmd_file (@additional_rmd_files){
      my $source_rmd = dirname(__FILE__) . "/" . trim($rmd_file);
      my $target_rmd = $result_dir . "/" . basename(trim($rmd_file));
      copy( $source_rmd, $target_rmd ) or die "Copy failed: $!";
    } 
  }
  
  if ( defined $config->{$section}{function_r_files} ) {
    my @function_r_files = split( ';', $config->{$section}{function_r_files} );

    open( my $rmd, ">$rfile" )     or die "Cannot create $rfile";
    open( my $rin, "<$rtemplate" ) or die "Cannot open $rtemplate";
    my $function_done = 0;
    while (<$rin>) {
      my $rinline = $_;
      $rinline = s/\r|\n//g;
      if ( $rinline =~ /^```\{/ ) {
        if ( !$function_done ) {
          print $rmd '```{r, include=FALSE} \n';
          for my $func_file (@function_r_files) {
            $func_file = dirname(__FILE__) . "/" . trim($func_file);
            open( my $fin, "<$func_file" ) or die "Cannot open $func_file";
            while (<$fin>) {
              my $finline = $_;
              $finline = s/\r|\n//g;
              print $rmd $finline . "\n";
            }
          }
          print $rmd '``` \n';
          $function_done = 1;
        }
      }
    }
  }
  else {
    copy( $rtemplate, $rfile ) or die "Copy failed: $!";
  }

  return($rfile);
}

1;
