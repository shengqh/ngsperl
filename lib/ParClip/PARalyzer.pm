#!/usr/bin/perl
package ParClip::PARalyzer;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_pz";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $genome2bit = get_param_file( $config->{$section}{genome2bit}, "genome2bit", 1 );
  my $mirna_db   = get_param_file( $config->{$section}{mirna_db},   "mirna_db",   1 );
  my $sorted_by_name = get_option( $config, $section, "sorted_by_name" );
  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    my @sample_files = @{ $raw_files{$sample_name} };
    my $bam_file     = $sample_files[0];

    my $fileType  = "BOWTIE_FILE";
    my $bam2sam   = "";
    my $rmcmd     = "";
    my $inputFile = $bam_file;
    if ( $bam_file =~ /bam$/ ) {
      $fileType  = "SAM_FILE";
      $inputFile = "${sample_name}.sortedName.sam";

      my $sortCmd;
      my $bamInput;
      my $rmList;
      if ($sorted_by_name) {
        $sortCmd  = "";
        $bamInput = $bam_file;
        $rmList = "";
      }
      else {
        $sortCmd  = "ln -s $bam_file ${sample_name}.bam; samtools sort -n -o ${sample_name}.sortedName.bam ${sample_name}.bam; rm ${sample_name}.bam;";
        $bamInput = "${sample_name}.sortedName.bam";
        $rmList = $bamInput; 
      }

      $bam2sam = "if [ ! -s $inputFile ]; then
  $sortCmd
  samtools view -h -F 4 -o $inputFile $bamInput
fi";
      $rmcmd = "if [ -s ${sample_name}.target.csv ]; then
  rm $inputFile $rmList 
fi";
    }

    if ( $bam_file =~ /sam$/ ) {
      $fileType = "SAM_FILE";
    }

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    my $iniFile = "${sample_name}.ini";
    open( INI, ">${cur_dir}/${iniFile}" ) or die "Cannot create ${cur_dir}/${iniFile}";
    print INI "
BANDWIDTH=3
CONVERSION=T>C
MINIMUM_READ_COUNT_PER_GROUP=5
MINIMUM_READ_COUNT_PER_CLUSTER=5
MINIMUM_READ_COUNT_FOR_KDE=5
MINIMUM_CLUSTER_SIZE=10
MINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER=2
MINIMUM_CONVERSION_COUNT_FOR_CLUSTER=1
MINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION=5
MINIMUM_READ_LENGTH=16
MAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES=5

ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL=5

$fileType=$inputFile
GENOME_2BIT_FILE=$genome2bit
FIND_MIRNA_SEEDMATCHES=$mirna_db

OUTPUT_DISTRIBUTIONS_FILE=${sample_name}.distribution.csv
OUTPUT_GROUPS_FILE=${sample_name}.groups.csv
OUTPUT_CLUSTERS_FILE=${sample_name}.cluster.csv
OUTPUT_MIRNA_TARGETS_FILE=${sample_name}.target.csv
    
";
    close(INI);

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $final_file = "${sample_name}.target.csv";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    print $pbs "

$bam2sam

PARalyzer $memory $iniFile

$rmcmd

";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all ParClip::PARalyzer tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    push( @result_files, "${result_dir}/${sample_name}/${sample_name}.cluster.csv" );
    push( @result_files, "${result_dir}/${sample_name}/${sample_name}.distribution.csv" );
    push( @result_files, "${result_dir}/${sample_name}/${sample_name}.groups.csv" );
    push( @result_files, "${result_dir}/${sample_name}/${sample_name}.target.csv" );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
