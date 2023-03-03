#!/usr/bin/perl
package QC::ChipseqQC;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::UniqueTask;
use Pipeline::PipelineUtils;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_cqc";
  $self->{_docker_prefix} = "chipqc_";
  bless $self, $class;
  return $self;
}

sub get_qctable {
  my ( $self, $config, $section, $task_name, $treatments ) = @_;
  my $result = {};
  if ( has_raw_files($config, $section, "qctable") ){
    $result = get_raw_files( $config, $section, "qctable" );
  } else {
    my $task = {};
    for my $treatment (keys %$treatments) {
      $task->{$treatment} = {
        Condition => $treatment,
        Replicate => 1
      };
    }
 
    $result->{$task_name} = $task;
  }
  return ($result);
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $bamfiles = get_raw_files( $config, $section );

  my $treatments = $config->{$section}{"groups"};
  my $controls = $config->{$section}{"inputs"};
  if ( !defined $controls ) {
    $controls = $config->{$section}{"controls"};
  }
  my $peaksfiles = get_raw_files( $config, $section, "peaks" );
  my $peakSoftware = get_option( $config, $section, "peak_software" );
  my $genome       = get_option( $config, $section, "genome", "unknown" );
  my $combined     = get_option( $config, $section, "combined" );

  my $chromosomes = get_option( $config, $section, "chromosomes", "" );
  
  my $paired_end = get_option( $config, $section, "is_paired_end");

  my $qctable = $self->get_qctable( $config, $section, $task_name, $treatments);

  my $script = dirname(__FILE__) . "/ChipseqQC.r";
  if ( !-e $script ) {
    die "File not found : " . $script;
  }

  my $target_script_name = $task_name . ".r";
  my $target_script = $result_dir . "/" . $target_script_name;
  `echo "setwd('$result_dir')\n\n" > $target_script`;
  `cat $script >> $target_script`;

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $expectFiles = $self->result( $config, $section );
  my @sortedKeys  = ( sort keys %$expectFiles );
  my $final_file  = $expectFiles->{ $sortedKeys[ scalar(@sortedKeys) - 1 ] }->[0];

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
  
  my $sourceBamFiles = $bamfiles;
  if ($paired_end){
    $sourceBamFiles = {};
    for my $bamName (keys %$bamfiles){
      $sourceBamFiles->{$bamName} = [$result_dir . "/" . $bamName . ".firstread.bam"];
    } 
  }

  my ($mapFiles, $condition_map) = writeDesignTable( $result_dir, $section, $qctable, $sourceBamFiles, $peaksfiles, $peakSoftware, $combined, $task_name, $treatments, $controls );

  $config->{$section}{parameterSampleFile1} = {
    task_name => $task_name,
    genome => $genome,
    chromosomes => $chromosomes,
    consensus => get_option($config, $section, "consensus", 1)
  };
  writeParameterSampleFile( $config, $section, $result_dir, 1, 0 );  

  $config->{$section}{parameterSampleFile2} = $mapFiles;
  writeParameterSampleFile( $config, $section, $result_dir, 2, 0 );  

  if ($combined) {
    my $mapFile=$mapFiles->{$task_name};
    my $mapFileName = basename($mapFile);
    my $rdataFile=$mapFileName . ".rdata";

    if ($paired_end){
      #for paired end data, keep the first read only for ChipQC
      print $pbs "
if [[ ! -s $rdataFile ]]; then
";
      for my $bamName (keys %$bamfiles){
        my $oldFile = $bamfiles->{$bamName}->[0];
        my $newFile = $bamName . ".firstread.bam";
        print $pbs "
  if [[ ! -s $newFile ]]; then
    samtools view -b -f 65 -o $newFile $oldFile
    samtools index $newFile
  fi
";
        $sourceBamFiles->{$bamName} = [$result_dir . "/" . $newFile ];
      } 
      print $pbs "
fi

";
    }
    
    print $pbs "R --vanilla -f $target_script_name\n\n";
    
    if ($paired_end){
      print $pbs "if [[ -s $final_file ]]; then 
";
      for my $bamName (keys %$sourceBamFiles){
        my $newFile = basename($sourceBamFiles->{$bamName}->[0]);
        print $pbs "  rm $newFile ${newFile}.bai
";
      } 
      print $pbs "fi

";
    }
  }
  else {
    for my $qcname ( sort keys %$mapFiles ) {
      my $mapFileName = $mapFiles->{$qcname};
      my $curdir      = $result_dir . "/" . $qcname;
      print $pbs "cd $curdir
R --vanilla -f $target_script_name \n";
    }
  }

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  my $treatments = $config->{$section}{"groups"};
  my $combined = get_option( $config, $section, "combined" );
  my $result = {};

  if ($combined) {
    my @result_files = ();
    push( @result_files, $result_dir . "/${task_name}.config.txt.rdata" );
    push( @result_files, $result_dir . "/${task_name}.config.txt.version" );
    my $targetDir    = $result_dir . "/ChIPQCreport";
    push( @result_files, $targetDir . "/ChIPQC.html" );
    push( @result_files, $targetDir . "/GenomicFeatureEnrichment.png" );
    push( @result_files, $targetDir . "/CCPlot.png" );
    push( @result_files, $targetDir . "/PeakCorHeatmap.png" );
    push( @result_files, $targetDir . "/PeakPCA.png" );
    push( @result_files, $targetDir . "/CoverageHistogramPlot.png" );
    $result->{$task_name} = filter_array( \@result_files, $pattern );
  }
  else {
    my $qctable = $self->get_qctable( $config, $section, $task_name, $treatments);
    for my $qcname ( sort keys %$qctable ) {
      if ( $qcname eq "Tissue" || $qcname eq "Factor" ) {
        next;
      }
      my @result_files = ();
      my $curdir       = $result_dir . "/" . $qcname;
      my $targetDir    = $curdir . "/ChIPQCreport";
      push( @result_files, $targetDir . "/ChIPQC.html" );
      push( @result_files, $targetDir . "/GenomicFeatureEnrichment.png" );
      push( @result_files, $targetDir . "/CCPlot.png" );
      push( @result_files, $targetDir . "/PeakCorHeatmap.png" );
      push( @result_files, $targetDir . "/PeakPCA.png" );
      push( @result_files, $targetDir . "/CoverageHistogramPlot.png" );
      $result->{$qcname} = filter_array( \@result_files, $pattern );
    }
  }
  return $result;
}

1;
