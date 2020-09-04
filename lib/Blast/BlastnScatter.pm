#!/usr/bin/perl
package Blast::BlastnScatter;

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
use Bio::SeqIO;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_bs";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  my %raw_files        = %{ get_raw_files( $config, $section ) };
  
  my $blastdb = get_option($config, $section, "localdb", 0);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files     = @{ $raw_files{$sample_name} };
    my $sample           = $sample_files[0];
    my $blast_result     = $sample_name . ".blastn.tsv";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $blast_result );
    print $pbs "
blastn -task blastn-short -db $blastdb/nt -perc_identity 100 -query $sample -num_threads $thread -outfmt '6 qlen nident qacc sallacc salltitles' | awk '\$1 == \$2 {print}' | cut -f3- | sort | uniq > $blast_result
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my $blast_result     = $sample_name . ".blastn.tsv";
    $result->{$sample_name} = filter_array( ["${result_dir}/${sample_name}.blastn.tsv"], $pattern );
  }

  return $result;
}
1;
