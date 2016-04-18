#!/usr/bin/perl
package CNV::CNVnator;

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
  $self->{_suffix} = "_cnvnator";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $binsize        = $config->{$section}{binsize}        or die "define ${section}::binsize first";
  my $chromosome_dir = $config->{$section}{chromosome_dir} or die "define ${section}::chromosome_dir first";

  my $genome    = $config->{$section}{genome};
  my $genomestr = "";

  if ( defined $genome ) {
    $genomestr = "-genome " . $genome;
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };

    my $bam_file = $sample_files[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $cur_dir  = create_directory_or_die( $result_dir . "/$sample_name" );
    my $rootFile = $sample_name . ".root";
    my $callFile = $sample_name . ".call";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $callFile );

    print $pbs "   
if [ ! -s $rootFile ]; then
  echo \"EXTRACTING READ MAPPING FROM BAM/SAM FILES =\" `date`
  cnvnator $genomestr -unique -root $rootFile -tree $bam_file 
fi

echo \"GENERATING HISTOGRAM =\" `date`
cnvnator $genomestr -root $rootFile -d $chromosome_dir -his $binsize 

echo \"CALCULATING STATISTICS =\" `date`
cnvnator -root $rootFile -stat $binsize 

echo \"RD SIGNAL PARTITIONING =\" `date`
cnvnator -root $rootFile -partition $binsize 

echo \"CNV CALLING =\" `date`
cnvnator -root $rootFile -call $binsize > $callFile
";

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all cnvnator tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my $cur_dir      = create_directory_or_die( $result_dir . "/$sample_name" );
    my $rootFile     = $sample_name . ".root";
    my $callFile     = $sample_name . ".call";
    my @result_files = ();
    push( @result_files, $cur_dir . "/" . $callFile );
    push( @result_files, $cur_dir . "/" . $rootFile );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }

  return $result;
}

1;
