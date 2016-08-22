#!/usr/bin/perl
package CNV::cnMops;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::UniqueTask;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_cnmops";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $cqstools = get_cqstools( $config, $section, 0 );

  my $bedfile = $config->{$section}{bedfile};

  my $refSeqNames = $config->{$section}{ref_seq_names};
  if ( !defined $bedfile && !defined $refSeqNames ) {
    die "Either bedfile or ref_seq_names should be defined!";
  }
  
  my $window = get_option($config, $section, "window", 100);

  my $isbamsorted = $config->{$section}{isbamsorted};
  if ( !defined($isbamsorted) ) {
    $isbamsorted = 0;
  }

  my $pairmode = $config->{$section}{pairmode};
  if ( !defined($pairmode) ) {
    $pairmode = "unpaired";
  }

  my $refnames = $config->{$section}{refnames};
  if ( !defined $refnames ) {
    $refnames = [];
  }

  my $rtemplate = dirname(__FILE__) . "/cnMops.r";
  if ( !-e $rtemplate ) {
    die "File not found : " . $rtemplate;
  }

  my $callFile = "${task_name}.call";

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $rfile = $result_dir . "/cnmops_${task_name}.r";
  open( my $r, ">$rfile" ) or die "Cannot create $rfile";
  print $r "setwd(\"$result_dir\")
callfile<-\"$callFile\"
prefix<-\"$task_name\"
pairmode<-\"$pairmode\"
parallel<-$thread
";
  if ( defined $bedfile ) {
    print $r "hasbed<-1
bedfile<-\"$bedfile\"
";
  }
  else {
    print $r "hasbed<-0
window<-$window
";
    
  }

  print $r "sample_names <- c( \n";
  my $isfirst = 1;
  for my $sample_name ( sort keys %raw_files ) {
    if ($isfirst) {
      print $r "\"$sample_name\"\n";
      $isfirst = 0;
    }
    else {
      print $r ",\"$sample_name\"\n";
    }
  }
  print $r ")

bam_files <- c(
";

  $isfirst = 1;
  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $bam_file     = $sample_files[0];

    if ( !$isbamsorted ) {
      ( $bam_file, my $bamSorted ) = get_sorted_bam($bam_file);
    }

    if ($isfirst) {
      print $r "\"$bam_file\"\n";
      $isfirst = 0;
    }
    else {
      print $r ",\"$bam_file\"\n";
    }
  }
  print $r ")
";

  $isfirst = 1;
  print $r "refnames<-c(";
  for my $refName ( @{$refnames} ) {
    if ($isfirst) {
      print $r "\"$refName\"";
      $isfirst = 0;
    }
    else {
      print $r ",\"$refName\"";
    }
  }
  print $r ")
";

  if ( defined $refSeqNames ) {
    print $r "refSeqNames<-c(";
    for my $refSeqName ( @{$refSeqNames} ) {
      if ($isfirst) {
        print $r "\"$refSeqName\"";
        $isfirst = 0;
      }
      else {
        print $r ",\"$refSeqName\"";
      }
    }
    print $r ")
";
  }

  open my $rt, "<$rtemplate" or die $!;
  while (<$rt>) {
    if ( $_ =~ '^#' ) {
      next;
    }
    last;
  }
  while (<$rt>) {
    s/\r|\n//g;
    print $r $_, "\n";
  }
  close($rt);
  close $r;

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc   = $cluster->get_log_description($log);
  my $final_file = "${task_name}.call";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  print $pbs "   
cd $pbs_dir
R --vanilla -f $rfile 
";

  if (defined $cqstools){
  print $pbs "
if [[ -s  $final_file && ! -s ${final_file}.tsv ]]; then 
  mono $cqstools cnmops_merge -i $final_file -o ${final_file}.tsv
fi
";
  }

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $cqstools = get_cqstools( $config, $section, 0 );

  my @result_files = ();
  push( @result_files, $result_dir . "/${task_name}.call" );
  if(defined $cqstools){
    push( @result_files, $result_dir . "/${task_name}.call.tsv" );
  }
  push( @result_files, $result_dir . "/${task_name}.call.bed" );
  push( @result_files, $result_dir . "/${task_name}.cnvr.tsv" );

  my $result = { $task_name => filter_array( \@result_files, $pattern ) };
  return $result;
}

1;
