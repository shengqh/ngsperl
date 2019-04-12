#!/usr/bin/perl
package Imputation::Shapeit;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_si";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my %mapFiles = %{ get_raw_files( $config, $section, "genetic_map_file" ) };

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sample       = $sample_files[0];

    my $map = "";
    if (defined $mapFiles{$sample_name}){
      my @mapFiles = @{ $mapFiles{$sample_name} };
      $map      = $mapFiles[0];
    }else{
      my ( $key ) = $sample_name =~ /_(chr)\S+$/;
      my @mapFiles = @{ $mapFiles{$key} };
      $map      = $mapFiles[0];
    }

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $source_option;

    if ( $sample =~ /bed$/ ) {
      my $fam = change_extension( $sample, ".fam" );
      my $bim = change_extension( $sample, ".bim" );
      $source_option = "-B $sample $bim $fam";
    }
    elsif ( $sample =~ /ped$/ ) {
      my $ped_map = change_extension( $sample, ".map" );
      $source_option = "-P $sample $ped_map";
    }
    else {
      my $gen_sample = change_extension( $sample, ".sample" );
      $source_option = "-G $sample $gen_sample";
    }

    my $phase_file  = "${sample_name}.haps";
    my $sample_file = "${sample_name}.sample";
    my $plink_sample_file = "${sample_name}.plink.sample";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $phase_file );
    print $pbs "shapeit $option $source_option -M $map -O $phase_file $sample_file

if [[ -s $sample_file ]]; then
  cut -d ' ' -f 1-3,6,7  $sample_file > $plink_sample_file
fi
";
    $self->close_pbs( $pbs, $pbs_file );
    
    print $sh "\$MYCMD ./$pbs_name \n";
  }
  print $sh "exit 0\n";
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();

    push( @result_files, "${result_dir}/${sample_name}.haps" );
    push( @result_files, "${result_dir}/${sample_name}.sample" );
    push( @result_files, "${result_dir}/${sample_name}.plink.sample" );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
