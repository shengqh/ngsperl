#!/usr/bin/perl
package Variants::GlmvcValidate;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::GroupTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_gv";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $glmvcfile = get_param_file( $config->{$section}{execute_file}, "execute_file", 1, not $self->using_docker() );
  my $source_type = $config->{$section}{source_type} or die "source_type is not defined in $section";

  my $raw_files = get_raw_files( $config, $section );

  my $validateFiles = get_raw_files( $config, $section, "validation_files" );

  my %group_sample_map = ();

  my $fafile           = "";
  my $mpileupParameter = "";
  my $isbam            = lc($source_type) eq "bam";
  if ($isbam) {
    $fafile = get_param_file( $config->{$section}{fasta_file}, "fasta_file (for mpileup)", 1 );
    $mpileupParameter = $config->{$section}{mpileup_option};
    if ( defined $mpileupParameter ) {
      if ( $mpileupParameter eq "" ) {
        undef($$mpileupParameter);
      }
    }

    my $groups = get_raw_files( $config, $section, "groups" );
    for my $group_name ( sort keys %{$groups} ) {
      my @samples = @{ $groups->{$group_name} };
      my @gfiles  = ();
      my $index   = 0;
      foreach my $sample_name (@samples) {
        my @samplebam_files = @{ $raw_files->{$sample_name} };
        push( @gfiles, $samplebam_files[0] );
      }
      $group_sample_map{$group_name} = \@gfiles;
    }
  }
  else {
    %group_sample_map = %{$raw_files};
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";
  print $sh "cd $pbs_dir \n";

  for my $group_name ( sort keys %group_sample_map ) {
    my $validateFile;
    if ( is_hash($validateFiles) ) {
      if ( defined $validateFiles->{$group_name} ) {
        $validateFile = $validateFiles->{$group_name}[0];
      }
      else {
        $validateFile = $validateFiles->{$task_name}[0];
      }
    }
    elsif ( is_array($validateFiles) ) {
      $validateFile = $validateFiles->[0];
    }
    else {
      $validateFile = $validateFiles;
    }

    my @sample_files = @{ $group_sample_map{$group_name} };
    my $sampleCount  = scalar(@sample_files);

    my $cur_dir = create_directory_or_die( $result_dir . "/$group_name" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $group_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir );

    my $final = "${group_name}.validation.tsv";
    if ($isbam) {
      if ( $sampleCount != 2 ) {
        die "SampleFile should be normal,tumor paired for " . $group_name . ".";
      }
      my $normal = $sample_files[0];
      my $tumor  = $sample_files[1];

      my $cmd;
      if ( defined $mpileupParameter ) {
        $cmd = "samtools mpileup -f $fafile $mpileupParameter $normal $tumor | mono-sgen $glmvcfile validate -t console $option -o ${cur_dir}/${group_name}.validation -v $validateFile";
      }
      else {
        $cmd = "mono $glmvcfile validate -c $thread -t bam -f $fafile $option --normal $normal --tumor $tumor -o ${cur_dir}/${group_name}.validation -v $validateFile";
      }

      print $pbs "
if [ -s $final ]; then
  echo job has already been done. if you want to do again, delete ${cur_dir}/${final} and submit job again.
  exit 0;
fi      
      
if [ ! -s ${normal}.bai ]; then
  samtools index ${normal}
fi

if [ ! -s ${tumor}.bai ]; then
  samtools index ${tumor}
fi

$cmd

";
    }
    else {
      print $pbs "mono $glmvcfile validate -t mpileup -m $sample_files[0] $option -o ${cur_dir}/${group_name}.validation -v $validateFile \n";
    }
    $self->close_pbs( $pbs, $pbs_file );

  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %group_sample_map = ();
  my $source_type      = $config->{$section}{source_type} or die "source_type is not defined in $section";
  my $isbam            = lc($source_type) eq "bam";
  if ($isbam) {
    %group_sample_map = %{ get_raw_files( $config, $section, "groups" ) };
  }
  else {
    %group_sample_map = %{ get_raw_files( $config, $section ) };
  }

  my $result = {};
  for my $group_name ( keys %group_sample_map ) {
    my @result_files = ();
    my $cur_dir      = $result_dir . "/$group_name";

    my $final = "$cur_dir/${group_name}.validation.tsv";

    push( @result_files, "$final" );
    $result->{$group_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
