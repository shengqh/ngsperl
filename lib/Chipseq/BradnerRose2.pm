#!/usr/bin/perl
package Chipseq::BradnerRose2;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::GroupTask;
use CQS::NGSCommon;
use Data::Dumper;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_br";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $genome = get_option( $config, $section, "genome" );
  my $pipeline_dir = get_directory( $config, $section, "pipeline_dir", 1 );
  my $binding_site_file = parse_param_file( $config, $section, "binding_site_file", 1 );

  if ( $option eq "" ) {
    $option = "-s 12500 -t 2500";
  }

  my %raw_files = %{ $self->get_grouped_raw_files( $config, $section, "groups" ) };
  my %control_files = %{ $self->get_grouped_raw_files( $config, $section, "controls" ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $group_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$group_name} };
    my $treatment = "-r " . $sample_files[0];
  
    print Dumper(@sample_files) . "\n";
    
    my @controls = @{ $control_files{$group_name} };

    print Dumper(@controls) . "\n";

    my $control = "-c " . $controls[0];

    my $cur_dir = create_directory_or_die( $result_dir . "/$group_name" );

    my $cpRawFile  = "${group_name}.copynumber";
    my $cpCallFile = "${group_name}.call";
    my $cpSegFile  = "${group_name}.segment";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $group_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc   = $cluster->get_log_description($log);
    my $final_file = "${group_name}_peaks_AllEnhancers.table.txt";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );
    print $pbs "
cd $pipeline_dir
python ROSE2_main.py -g $genome -i $binding_site_file $treatment $control -o $cur_dir $option
";

    $self->close_pbs( $pbs, $pbs_file );

    print $sh "\$MYCMD ./$pbs_name \n";
  }

  print $sh "exit 0\n";
  close $sh;

  print "!!!shell file $shfile created, you can run this shell file to submit tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %raw_files = %{ $self->get_grouped_raw_files( $config, $section, "groups" ) };

  my $result = {};
  for my $group_name ( sort keys %raw_files ) {
    my $cur_dir      = $result_dir . "/$group_name";
    my @result_files = ();
    push( @result_files, $cur_dir . "/${group_name}_peaks_AllEnhancers.table.txt" );
    $result->{$group_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
