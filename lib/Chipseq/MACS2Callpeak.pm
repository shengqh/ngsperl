#!/usr/bin/perl
package Chipseq::MACS2Callpeak;

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
use Data::Dumper;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_mc";
  bless $self, $class;
  return $self;
}

sub get_raw_files_with_or_without_groups{
  my ( $self, $config, $section ) = @_;
  my $result;
  if ( has_raw_files( $config, $section, "groups" ) ) {
    $result = get_grouped_raw_files( $config, $section, "groups" );
  }
  else {
    $result = get_raw_files( $config, $section );
  }
  
  return $result;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my %raw_files = %{$self->get_raw_files_with_or_without_groups($config, $section)};

  my %control_files;
  my $has_control = 0;
  if ( has_raw_files( $config, $section, "groups" ) ) {
    $has_control = has_raw_files( $config, $section, "controls" );
    if ($has_control) {
      %control_files = %{ get_grouped_raw_files( $config, $section, "controls" ) };
    }
  }

  my $peak_name = ( $option =~ /--broad/ ) ? "broadPeak" : "narrowPeak";

  my $outputBigwig      = get_option( $config, $section, "output_bigwig",   0 );
  my $chr_size_file ="";
  if ($outputBigwig) {
    $chr_size_file = $config->{$section}{chr_size_file} or die "define ${section}::chr_size_file first to be used in output_bigwig";
  }


  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    my @sample_files = @{ $raw_files{$sample_name} };
    my $treatment = "-t " . join( " ", @sample_files );

    my $control = "";
    if ($has_control) {
      my @control_files = @{ $control_files{$sample_name} };
      $control = "-c " . join( " ", @control_files );
    }
    my $final_file = "${sample_name}_peaks.${peak_name}";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file, undef, 1 );

    print $pbs "   
rm -f $sample_name.macs2.failed $sample_name.macs2.succeed 

macs2 callpeak $option $treatment $control -n $sample_name

status=\$?
if [[ \$status -ne 0 ]]; then
  touch $sample_name.macs2.failed
  rm -f ${sample_name}_peaks.${peak_name}
else
  touch $sample_name.macs2.succeed
  cut -f1-6 ${sample_name}_peaks.${peak_name} > ${sample_name}_peaks.${peak_name}.bed
fi

";
      if ($outputBigwig) {
    $chr_size_file = $config->{$section}{chr_size_file} or die "define ${section}::chr_size_file first to be used in output_bigwig";
  }
    if ($outputBigwig) {
    print $pbs "
bedGraphToBigWig ${sample_name}_control_lambda.bdg $chr_size_file ${sample_name}_control_lambda.bigwig
bedGraphToBigWig ${sample_name}_treat_pileup.bdg $chr_size_file ${sample_name}_treat_pileup.bigwig
";
    }

    $self->close_pbs( $pbs, $pbs_file );

    print $sh "\$MYCMD ./$pbs_name \n";
  }

  print $sh "exit 0\n";
  close $sh;
  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my %raw_files = %{$self->get_raw_files_with_or_without_groups($config, $section)};

  my $peak_name = ( $option =~ /--broad/ ) ? "broadPeak" : "narrowPeak";

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my $cur_dir      = $result_dir . "/$sample_name";
    my @result_files = ();
    push( @result_files, $cur_dir . "/${sample_name}_treat_pileup.bdg" );
    push( @result_files, $cur_dir . "/${sample_name}_control_lambda.bdg" );
    push( @result_files, $cur_dir . "/${sample_name}_peaks.${peak_name}" );
    push( @result_files, $cur_dir . "/${sample_name}_peaks.${peak_name}.bed" );
    #push( @result_files, $cur_dir . "/${sample_name}_summits.bed" );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my %raw_files = %{$self->get_raw_files_with_or_without_groups($config, $section)};

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
      $result->{$sample_name} = $self->get_pbs_filename( $pbs_dir, $sample_name );
  }

  #print(Dumper($result));
  return $result;
}

1;
