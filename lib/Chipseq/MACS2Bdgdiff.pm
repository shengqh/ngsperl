#!/usr/bin/perl
package Chipseq::MACS2Bdgdiff;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::GroupTask;
use CQS::StringUtils;

our @ISA = qw(CQS::GroupTask);

my $directory;

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_mb";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %group_sample_treat   = %{ get_raw_files( $config, $section, "source", "_treat_pileup.bdg" ) };
  my %group_sample_control = %{ get_raw_files( $config, $section, "source", "_control_lambda.bdg" ) };
  my $comparisons = get_raw_files( $config, $section, "groups" );

  #print Dumper(%group_sample_treat);
  #print Dumper(%group_sample_control);

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $comparisonName ( sort keys %{$comparisons} ) {
    my @group_names = @{ $comparisons->{$comparisonName} };
    my $groupCount  = scalar(@group_names);
    if ( $groupCount != 2 ) {
      die "Comparison should be control,treatment paired.";
    }

    my $condition1treat   = $group_sample_treat{ $group_names[0] }->[0];
    my $condition1control = $group_sample_control{ $group_names[0] }->[0];
    my $condition2treat   = $group_sample_treat{ $group_names[1] }->[0];
    my $condition2control = $group_sample_control{ $group_names[1] }->[0];

    my $cur_dir = create_directory_or_die( $result_dir . "/$comparisonName" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $comparisonName );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $comparisonName );

    my $log_desc = $cluster->get_log_description($log);
    my $final_file = "${comparisonName}_c3.0_common.bed";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    print $pbs "
macs2 bdgdiff $option --t1 $condition1treat --t2 $condition2treat --c1 $condition1control --c2 $condition2control --o-prefix $comparisonName  
";

    $self->close_pbs( $pbs, $pbs_file );

    print $sh "\$MYCMD ./$pbs_name \n";
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $comparisons = get_raw_files( $config, $section, "groups" );
  my $result = {};
  for my $comparisonName ( sort keys %{$comparisons} ) {
    my $cur_dir = $result_dir . "/$comparisonName";

    my @result_files = ();
    push( @result_files, $cur_dir . "/${comparisonName}_c3.0_cond1.bed" );
    push( @result_files, $cur_dir . "/${comparisonName}_c3.0_cond2.bed" );
    push( @result_files, $cur_dir . "/${comparisonName}_c3.0_common.bed" );
    $result->{$comparisonName} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
