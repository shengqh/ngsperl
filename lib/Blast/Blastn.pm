#!/usr/bin/perl
package Blast::Blastn;

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
  $self->{_suffix} = "_bn";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  my %raw_files        = %{ get_raw_files( $config, $section ) };
  my $blastn_interpret = dirname(__FILE__) . "/blastn-interpret.pl";
  
  my $blastdb = get_option($config, $section, "localdb", 0);
  my $blastn  = dirname(__FILE__);
  if($blastdb){
    $option = $option . " --localdb " . $blastdb;
    $blastn = $blastn . "/blastn-short-local.pl";
  }else{
    $blastn = $blastn . "/blastn-short-remote.pl";
  }

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files     = @{ $raw_files{$sample_name} };
    my $sample           = $sample_files[0];
    my $blast_result     = $sample_name . ".blastn.tsv";
    my $interpret_result = $sample_name . ".table.tsv";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $interpret_result );
    print $pbs "
if [ ! -s $blast_result ]; then
  perl $blastn $option -i $sample -o $blast_result -t $thread
fi

perl $blastn_interpret -i $blast_result -o $interpret_result
";
    $self->close_pbs( $pbs, $pbs_file );
    print $sh "\$MYCMD ./$pbs_name \n";
  }

  print $sh "exit 0\n";
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my @result_files = ();

    my $interpret_result = $sample_name . ".table.tsv";
    push( @result_files, "${result_dir}/${sample_name}.blastn.tsv" );
    push( @result_files, "${result_dir}/${sample_name}.table.tsv" );
    push( @result_files, "${result_dir}/${sample_name}.table.sequences.tsv" );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }

  return $result;
}
1;
