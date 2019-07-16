#!/usr/bin/perl
package Samtools::DepthStat;

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
  $self->{_suffix} = "_dps";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %group_sample_map = %{ get_group_sample_map( $config, $section ) };

  my $minimum_depth = $config->{$section}{minimum_depth};
  my $cqscommand = "";
  if ( defined $minimum_depth ) {
    $cqscommand = " | cqstools depth_filter -d $minimum_depth";
  }
  else {
    $minimum_depth = 1;
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc = $cluster->get_log_description($log);
  my $final    = "${task_name}.tsv";

  my $template_r_file = dirname(__FILE__) . "/DepthStat.r";
  my $target_r_file   = $result_dir . "/" . $task_name . ".r";
  open( my $targetr,   ">$target_r_file" ) or die "Cannot create file $target_r_file";
  open( my $rtemplate, $template_r_file )  or die "Cannot open file $template_r_file";
  print $targetr "resultDir=\"$result_dir\"
inputFile=\"$final\"
outputFile=\"${task_name}.png\"
minimumDepth=$minimum_depth
";

  while ( my $line = <$rtemplate> ) {
    chomp($line);
    $line =~ s/\r//g;
    if ( $line =~ /^#/ ) {
      next;
    }
    print $targetr $line, "\n";
  }

  close($rtemplate);
  close($targetr);

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final );
  
  print $pbs "echo \"Group\tSample\tCount\" >> $final \n";
  
  my @group_names = sort keys %group_sample_map;
  for my $group_name (@group_names) {
    my @samples = @{ $group_sample_map{$group_name} };
    for my $sample (@samples) {
      my @sample_files = @{$sample};
      my $sample_name  = shift @sample_files;
      my $samples      = join( " ", @sample_files );

      print $pbs "echo processing $sample_name ...\n";
      print $pbs "samtools depth $option $samples $cqscommand | wc | awk '{print \"${group_name}\\t${sample_name}\\t\" \$1;}'>> $final \n";
    }
  }

  if ( $config->{$section}{is_groups_paired} ) {
    my @first_samples = @{ $group_sample_map{ $group_names[0] } };
    my $sample_count  = scalar(@first_samples);

    for ( my $index = 0 ; $index < $sample_count ; $index++ ) {
      my @cursample_names = ();
      my @cursamples      = ();
      for my $group_name (@group_names) {
        my @samples      = @{ $group_sample_map{$group_name} };
        my @sample_files = @{ $samples[$index] };
        my $sample_name  = shift @sample_files;
        my $samples      = join( " ", @sample_files );
        push( @cursample_names, $sample_name );
        push( @cursamples,      $samples );
      }

      my $show_sample_names = join( "-", @cursample_names );
      my $show_samples      = join( " ", @cursamples );
      print $pbs "echo processing $show_sample_names ...\n";
      print $pbs "samtools depth $option $show_samples $cqscommand | wc | awk '{print \"Common\\t${show_sample_names}\\t\" \$1;}'>> $final \n";
    }
  }

  print $pbs "R --vanilla -f $target_r_file \n";
  $self->close_pbs( $pbs, $pbs_file );

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $result       = {};
  my @result_files = ();
  push( @result_files, "$result_dir/${task_name}.tsv" );
  push( @result_files, "$result_dir/${task_name}.png" );
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
