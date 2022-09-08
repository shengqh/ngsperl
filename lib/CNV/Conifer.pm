#!/usr/bin/perl
package CNV::Conifer;

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
  $self->{_name}   = "CNV::Conifer";
  $self->{_suffix} = "_conifer";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $conifer = $config->{$section}{conifer} or die "define conifer program location first.\nconifer => \"location\"";
  my $bedfile = $config->{$section}{bedfile} or die "define $section:bedfile first";

  my $isbamsorted = $config->{$section}{isbamsorted};
  if ( !defined($isbamsorted) ) {
    $isbamsorted = 0;
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  create_directory_or_die( $result_dir . "/rpkm" );
  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $bam_file = $sample_files[0];

    if ( !$isbamsorted ) {
      ( $bam_file, my $bamSorted ) = get_sorted_bam($bam_file);

      #print $bam_file . "\n";
    }
    my $final_file = "rpkm/" . $sample_name . ".rpkm";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    print $pbs "   
python3 $conifer rpkm --probes $bedfile --input $bam_file --output $final_file 
";

    $self->close_pbs( $pbs, $pbs_file );

  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  my $pbs_file2 = $self->get_pbs_filename( $pbs_dir, $task_name . "_after_rpkm" );
  my $pbs_name  = basename($pbs_file2);
  my $hdf5File  = "${task_name}.hdf5";
  my $svalsFile = "${task_name}.svals";
  my $plotFile  = "${task_name}.png";
  my $sdFile    = "${task_name}.sd";
  my $callFile  = "${task_name}.call";

  my $log      = "${log_dir}/${task_name}_after_rpkm.log";
  my $log_desc = $cluster->get_log_description($log);

  create_directory_or_die( $result_dir . "/call_images" );

  my $pbs2 = $self->open_pbs( $pbs_file2, $pbs_desc, $log_desc, $path_file, $result_dir, $callFile );

  print $pbs2 "
#2 analysis
echo analyze=`date`
python3 $conifer analyze --probes $bedfile --rpkm_dir rpkm/ --output $hdf5File --svd 6 --write_svals $svalsFile --plot_scree $plotFile --write_sd $sdFile 

#3 call
echo call=`date`
python3 $conifer call --input $hdf5File --output $callFile 

#4 plot
python3 $conifer plotcalls --input $hdf5File --calls $callFile --output call_images 
";

  $self->close_pbs( $pbs2, $pbs_file2 );

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $result = { $task_name => [ $result_dir . "/${task_name}.call" ] };

  return $result;
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %fqFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( sort keys %fqFiles ) {
    $result->{$sample_name} = $self->get_pbs_filename( $pbs_dir, $sample_name );
  }

  $result->{$task_name} = $self->get_pbs_filename( $pbs_dir, $task_name . "_after_rpkm" );

  return $result;
}

1;
