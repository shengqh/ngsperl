#!/usr/bin/perl
package CQS::ConfigUtils;

use strict;
use warnings;
use File::Basename;
use CQS::FileUtils;
use CQS::PBS;
use CQS::ClassFactory;
use CQS::StringUtils;
use CQS::CQSDebug;

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(get_parameter get_param_file parse_param_file get_raw_files get_raw_files2 get_run_command get_option_value)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub get_parameter {
  my ( $config, $section ) = @_;

  die "no section $section found!" if !defined $config->{$section};

  my $task_name = $config->{general}{task_name} or die "define general::task_name first";
  my $path_file = get_param_file( $config->{general}{path_file}, "path_file", 0 );
  if ( defined $path_file && -e $path_file ) {
    $path_file = "source $path_file";
  }
  else {
    $path_file = "";
  }

  my $refPbs     = $config->{$section}{pbs}        or die "define ${section}::pbs parameters first";
  my $target_dir = $config->{$section}{target_dir} or die "define ${section}::target_dir parameters first";
  my ( $logDir, $pbsDir, $resultDir ) = init_dir($target_dir);
  my ($pbsDesc) = get_pbs_desc($refPbs);

  die "define ${section}::option first" if ( !defined $config->{$section}{option} );
  my $option    = $config->{$section}{option};
  my $sh_direct = $config->{$section}{sh_direct};
  if ( !defined $sh_direct ) {
    $sh_direct = 0;
  }

  return ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct );
}

#get parameter which indicates a file. If required, not defined or not exists, die. If defined but not exists, die.
#returned file either undefined or exists.
sub get_param_file {
  my ( $file, $name, $required ) = @_;

  my $result = $file;

  if ($required) {
    if ( !defined $file ) {
      die "$name was not defined!";
    }

    if ( !is_debug() && !-e $file ) {
      die "$name $file defined but not exists!";
    }
  }
  else {
    if ( defined($file) ) {
      if ( $file eq "" ) {
        undef($result);
      }
      elsif ( !is_debug() && !-e $file ) {
        die "$name $file defined but not exists!";
      }
    }
  }
  return ($result);
}

sub parse_param_file {
  my ( $config, $section, $key, $required ) = @_;

  die "section $section was not defined!" if !defined $config->{$section};
  die "parameter key must be defined!" if !defined $key;

  if ( defined $config->{$section}{$key} ) {
    return $config->{$section}{$key};
  }

  my $key_ref = $key . "_ref";
  if ( defined $config->{$section}{$key_ref} ) {
    my $refSectionName = $config->{$section}{$key_ref};
    my $pattern;
    if ( ref($refSectionName) eq 'ARRAY' ) {
      my @parts = @{$refSectionName};
      if ( scalar(@parts) == 2 ) {
        $pattern        = $parts[1];
        $refSectionName = $parts[0];
      }
      else {
        $refSectionName = $parts[0];
      }
    }
    die "section $refSectionName was not defined!" if !defined $config->{$refSectionName};
    if ( defined $config->{$refSectionName}{class} ) {
      my $myclass = instantiate( $config->{$refSectionName}{class} );
      my $result = $myclass->result( $config, $refSectionName, $pattern );
      foreach my $k ( sort keys %{$result} ) {
        my @files = @{ $result->{$k} };
        return $files[0];
      }
    }
  }

  if ($required) {
    die "define ${section}::${key} first.";
  }

  return undef;
}

sub do_get_raw_files {
  my ( $config, $section, $returnself, $mapname, $pattern ) = @_;

  die "section $section was not defined!" if !defined $config->{$section};

  if ( !defined $mapname ) {
    $mapname = "source";
  }
  my $mapname_ref = $mapname . "_ref";

  if ( defined $config->{$section}{$mapname} ) {
    return ( $config->{$section}{$mapname}, 1 );
  }

  if ( defined $config->{$section}{$mapname_ref} ) {
    my $refSectionName = $config->{$section}{$mapname_ref};
    if ( ref($refSectionName) eq 'ARRAY' ) {
      my @parts = @{$refSectionName};
      if ( scalar(@parts) == 2 ) {
        $pattern        = $parts[1];
        $refSectionName = $parts[0];
      }
      else {
        $refSectionName = $parts[0];
      }
    }
    die "section $refSectionName was not defined!" if !defined $config->{$refSectionName};
    if ( defined $config->{$refSectionName}{class} ) {
      my $myclass = instantiate( $config->{$refSectionName}{class} );
      return ( $myclass->result( $config, $refSectionName, $pattern ), 0 );
    }
    else {
      my ( $result, $issource ) = do_get_raw_files( $config, $refSectionName, 1 );
      return ( $result, 0 );
    }
  }

  if ($returnself) {
    return ( $config->{$section}, 0 );
  }
  else {
    die "define $mapname or $mapname_ref for $section";
  }
}

sub get_raw_files {
  my ( $config, $section, $mapname, $pattern ) = @_;
  my ( $result, $issource ) = do_get_raw_files( $config, $section, 0, $mapname, $pattern );
  return $result;
}

#return raw files and if the raw files are extracted from source directly
sub get_raw_files2 {
  my ( $config, $section, $mapname, $pattern ) = @_;
  return do_get_raw_files( $config, $section, 0, $mapname, $pattern );
}

sub get_run_command {
  my $sh_direct = shift;
  if ($sh_direct) {
    return ("MYCMD=\"bash\" \n");
  }
  else {
    return ("type -P qsub &>/dev/null && MYCMD=\"qsub\" || MYCMD=\"bash\" \n");
  }
}

sub get_option_value{
  my ($value, $defaultValue) = @_;
  if(!defined $value){
    return ($defaultValue);
  }
  else{
    return ($value);
  }
}

1;

