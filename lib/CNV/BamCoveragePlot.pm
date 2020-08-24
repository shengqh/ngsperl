#!/usr/bin/perl
package CNV::BamCoveragePlot;

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
  $self->{_suffix} = "_bcp";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $genome_length_file = parse_param_file( $config, $section, "genome_length_file", 1 );
  my $window = get_option( $config, $section, "window", 100000 );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $rtemplate = dirname(__FILE__) . "/BamCoveragePlot.r";
  if ( !-e $rtemplate ) {
    die "File not found : " . $rtemplate;
  }

  open( my $length_file, "<$genome_length_file" ) or die;
  my $firstline = <$length_file>;
  my %genome_length = map { chomp; split /\t/ } (<$length_file>);
  close($length_file);

  my $bedfile = $self->get_name( $task_name, ".bed", 0 );
  my $resbedfile = $result_dir . "/" . $bedfile;

  open( my $bed, ">$resbedfile" ) or die "Cannot create $resbedfile";
  for my $chrom ( sort keys %genome_length ) {
    my $chrom_length = $genome_length{$chrom};
    my $start        = 0;
    my $end          = $start + $window;
    while ( $end <= $chrom_length ) {
      print $bed "$chrom\t$start\t$end\n";
      $start = $end;
      $end   = $end + $window;
    }
  }
  close($bed);

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $covfile = $bedfile . ".cov";

  my $rfile = "${result_dir}/${bedfile}.r";
  open( my $r, ">$rfile" ) or die "Cannot create $rfile";
  print $r "setwd(\"$result_dir\")
covfile<-\"$covfile\"
";

  print $r "sample_names <- c( \n";
  my $isfirst      = 1;
  my @sample_names = sort keys %raw_files;
  for my $sample_name (@sample_names) {
    if ($isfirst) {
      print $r "\"$sample_name\"\n";
      $isfirst = 0;
    }
    else {
      print $r ",\"$sample_name\"\n";
    }
  }
  print $r ") \n";

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

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc   = $cluster->get_log_description($log);
  my $final_file = "${covfile}.pdf";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  print $pbs "   
bedtools multicov -bams ";
  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $bam_file     = $sample_files[0];

    print $pbs " \"$bam_file\"";
  }
  print $pbs " -bed $bedfile > $covfile 

R --vanilla -f $rfile 
";

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $chrom = get_option( $config, $section, "chrom" );
  my $bedfile = $self->get_name( $task_name, "." . $chrom . ".bed", 0 );

  my @result_files = ();
  push( @result_files, $result_dir . "/$bedfile" );
  push( @result_files, $result_dir . "/$bedfile.cov" );
  push( @result_files, $result_dir . "/$bedfile.cov.pdf" );
  my $result = { $task_name => filter_array( \@result_files, $pattern ) };
  return $result;
}

1;
