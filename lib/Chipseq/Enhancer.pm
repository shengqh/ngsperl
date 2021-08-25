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
use CQS::UniqueTask;
use CQS::NGSCommon;
use Data::Dumper;

our @ISA = qw(CQS::UniqueTask);

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $pipeline_dir = get_directory( $config, $section, "pipeline_dir", 1 );
  my $genome          = get_option( $config, $section, "genome" );
  my $genomeDirectory = get_option( $config, $section, "genome_path" );
  my $gseaPath        = get_option( $config, $section, "gsea_path" );
  my $gmxPath         = get_option( $config, $section, "gmx_path" );
  my $cpgPath         = get_option( $config, $section, "cpg_path" );
  my $activityFile    = get_option( $config, $section, "activity_file", "" );
  if ( $activityFile ne "" ) {
    $option = $option . " -a $activityFile";
  }

  my $rawFiles = get_grouped_raw_files( $config, $section, "treatments" );

  my $peakCommand = "";
  my $peakFile;
  if ( defined $config->{$section}{peaks} && -e $config->{$section}{peaks} ) {    #it is a file
    $peakFile = $config->{$section}{peaks};
  }
  else {
    $peakFile = $task_name . "_merged_peaks.bed";
    my $tempFile  = $task_name . ".tmp.bed";
    my $peakfiles = get_raw_files( $config, $section, "peaks" );
    my @bedfiles  = ();
    for my $sample_name ( sort keys %$peakfiles ) {
      my $sample_files = $peakfiles->{$sample_name};
      push @bedfiles, @$sample_files;
    }
    my $fileOption = join( " ", @bedfiles );

    $peakCommand = "if [ ! -s $peakFile ]; then
  cat $fileOption | sort -k1,1 -k2,2n > $tempFile
  bedtools merge -i $tempFile -c 4 -o collapse > $peakFile
  rm $tempFile 
fi

";
  }

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

  print $pbs "$peakCommand \n";

  for my $fileName ( keys %$rawFiles ) {
    my @bamFiles = @{ $rawFiles->{$fileName} };
    my $bamFile = join( ' ', @bamFiles );

    my $final_file = "${fileName}/${fileName}_GENE_TABLE.txt";

    print $pbs "if [ ! -s $final_file ]; then
  python3 $pipeline_dir/enhancerPromoter.py $option -b $bamFile -i ${result_dir}/$peakFile -g $genome -o . --name $fileName --genomeDirectory $genomeDirectory --gseaPath $gseaPath --gmxPath $gmxPath --cpgPath $cpgPath
fi

";
  }
  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $rawFiles = get_grouped_raw_files( $config, $section, "treatments" );

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
