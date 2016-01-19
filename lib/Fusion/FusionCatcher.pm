#!/usr/bin/perl
package Fusion::FusionCatcher;

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
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_fc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $data_dir   = get_directory( $config, $section, "data_dir", 1 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files    = @{ $raw_files{$sample_name} };
    my $files = join(",", @sample_files);

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log     = $self->get_log_filename( $log_dir, $sample_name );

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );
    
    my $final = $cur_dir . "/final-list_candidate-fusion-genes.txt";

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    open( my $out, ">$pbs_file" ) or die $!;
    print $out "$pbs_desc
$log_desc

$path_file

if [ -s $final ]; then
  echo job has already been done. if you want to do again, delete $final and submit job again.
  exit 0
fi

cd $cur_dir

echo FusionCatcher_start=`date` 

fusioncatcher -p $thread -d $data_dir -i $files -o .

echo FusionCatcher_end=`date` 
";

    close $out;

    print "$pbs_file created\n";
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );
    
    my $final =  $result_dir . "/${sample_name}/final-list_candidate-fusion-genes.txt";
    my @result_files = ();
    push( @result_files, $final );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
