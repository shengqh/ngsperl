#!/usr/bin/perl
package Chipseq::Enhancer;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::Task;
use CQS::NGSCommon;
use Data::Dumper;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_ec";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $pipeline_dir = get_directory( $config, $section, "pipeline_dir", 1 );
  my $genome          = get_option( $config, $section, "genome" );
  my $genomeDirectory = get_option( $config, $section, "genome_path" );
  my $gseaPath        = get_option( $config, $section, "gsea_path" );
  my $gmxPath         = get_option( $config, $section, "gmx_path" );
  my $cpgPath         = get_option( $config, $section, "cpg_path" );
  my $activityFile    = get_option( $config, $section, "activity_file", "" );
  if($activityFile ne ""){
    $option = $option . " -a $activityFile";
  }

  my $rawFiles = get_raw_files( $config, $section );
  my %peakFiles;
  if ( defined $config->{$section}{peaks} && -e $config->{$section}{peaks} ) {
    for my $fileName ( keys %$rawFiles ) {
      $peakFiles{$fileName} = $config->{$section}{peaks};
    }
  }
  else {
    %peakFiles = %{ get_raw_files( $config, $section, "peaks" ) };
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $fileName ( keys %$rawFiles ) {
    my @bamFiles = @{ $rawFiles->{$fileName} };
    my $bamFile = join( ' ', @bamFiles );

    my @peaks = @{ $peakFiles{$fileName} };

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $fileName );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $fileName );

    my $log_desc   = $cluster->get_log_description($log);
    my $final_file = "${fileName}/${fileName}_GENE_TABLE.txt";
    my $pbs        = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
    my $peak       = $peaks[0];

    print $pbs "python $pipeline_dir/enhancerPromoter.py $option -b $bamFile -i $peak -g $genome -o . --name $fileName --genomeDirectory $genomeDirectory --gseaPath $gseaPath --gmxPath $gmxPath --cpgPath $cpgPath \n";

    $self->close_pbs( $pbs, $pbs_file );
    print $sh "\$MYCMD ./$pbs_name \n";
  }

  print $sh "exit 0\n";
  close $sh;

  print "!!!shell file $shfile created, you can run this shell file to submit tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $rawFiles = get_raw_files( $config, $section );

  my $result = {};
  for my $fileName ( keys %$rawFiles ) {
    my @result_files = ();
    my $final_file   = "${fileName}/${fileName}_GENE_TABLE.txt";
    push( @result_files, $result_dir . "/$final_file" );
    $result->{$fileName} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
