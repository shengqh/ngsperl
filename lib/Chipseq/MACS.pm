#!/usr/bin/perl
package Chipseq::MACS;

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
  $self->{_suffix} = "_macs";
  $self->{_group_keys} = ["groups", "inputs", "controls"];
  $self->{_use_tmp_folder} = 1;
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my %treatment_files = %{ get_grouped_raw_files( $config, $section, "groups" ) };

  my %control_files;
  if ( has_raw_files( $config, $section, "inputs" ) ) {
    %control_files = %{ get_grouped_raw_files( $config, $section, "inputs" ) };
  }
  elsif ( has_raw_files( $config, $section, "controls" ) ) {
    %control_files = %{ get_grouped_raw_files( $config, $section, "controls" ) };
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %treatment_files ) {
    my $cur_dir    = create_directory_or_die( $result_dir . "/$sample_name" );
    my $final_file = "${sample_name}_peaks.name.bed";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    my @sample_files = @{ $treatment_files{$sample_name} };

    my $localized_files = [];
    @sample_files = @{$self->localize_files_in_tmp_folder($pbs, \@sample_files, $localized_files, [".bai"])};
    my $treatment = "-t " . join( ",", @sample_files );

    my $control = "";
    if (%control_files) {
      if ( defined $control_files{$sample_name} ) {
        my @c_files = @{ $control_files{$sample_name} };

         @c_files = @{$self->localize_files_in_tmp_folder($pbs, \@c_files, $localized_files, [".bai"])};
        $control = "-c " . join( ",", @c_files );
      }
    }

    my $sname = $sample_name;
    $sname =~ s/ /_/g;

    print $pbs "
rm -f $sample_name.macs.failed $sample_name.macs.succeed

macs $option $treatment $control -n $sample_name

status=\$?
if [[ \$status -ne 0 ]]; then
  touch $sample_name.macs.failed
  rm -f ${sample_name}_peaks.bed
else
  touch $sample_name.macs.succeed
  sed 's/\\tMACS_peak_/\\t${sname}_/' ${sample_name}_peaks.bed > ${sample_name}_peaks.name.bed
fi
";

    $self->clean_temp_files($pbs, $localized_files);

    $self->close_pbs( $pbs, $pbs_file );

    print $sh "\$MYCMD ./$pbs_name \n";
  }

  print $sh "exit 0\n";
  close $sh;

  print "!!!shell file $shfile created, you can run this shell file to submit tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %treatment_files = %{ get_grouped_raw_files( $config, $section, "groups" ) };

  my $result = {};
  for my $group_name ( sort keys %treatment_files ) {
    my $cur_dir      = $result_dir . "/$group_name";
    my @result_files = ();
    push( @result_files, $cur_dir . "/${group_name}_peaks.name.bed" );
    push( @result_files, $cur_dir . "/${group_name}_MACS_wiggle/treat" );

    $result->{$group_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
