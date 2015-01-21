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
  $self->{_name} = "ParClip::PARalyzer";
  $self->{_suffix} = "_pz";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $genome2bit = get_param_file( $config->{$section}{genome2bit}, "genome2bit", 1 );
  my $mirna_db = get_param_file( $config->{$section}{mirna_db}, "mirna_db", 1 );
  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my $curDir      = create_directory_or_die( $resultDir . "/$sampleName" );

    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $bamfile     = $sampleFiles[0];

    my $pbsFile = $self->pbsfile($pbsDir, $sampleName);
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    my $iniFile = "${sampleName}.ini";
    open( INI, ">${curDir}/${iniFile}" ) or die "Cannot create ${curDir}/${iniFile}";
    print INI "
BANDWIDTH=3
CONVERSION=T>C
MINIMUM_READ_COUNT_PER_GROUP=5
MINIMUM_READ_COUNT_PER_CLUSTER=5
MINIMUM_READ_COUNT_FOR_KDE=5
MINIMUM_CLUSTER_SIZE=10
MINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER=1
MINIMUM_CONVERSION_COUNT_FOR_CLUSTER=1
MINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION=5
MINIMUM_READ_LENGTH=13
MAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES=0

ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL=5

BOWTIE_FILE=$bamfile
GENOME_2BIT_FILE=$genome2bit
FIND_MIRNA_SEEDMATCHES=$mirna_db

OUTPUT_DISTRIBUTIONS_FILE=${sampleName}.distribution.csv
OUTPUT_GROUPS_FILE=${sampleName}.groups.csv
OUTPUT_CLUSTERS_FILE=${sampleName}.cluster.csv
OUTPUT_MIRNA_TARGETS_FILE=${sampleName}.target.csv
    
";
    close(INI);

    print SH "\$MYCMD ./$pbsName \n";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

cd $curDir

echo PARalyzer_started=`date`

PARalyzer $iniFile

echo PARalyzer_finished=`date`

exit 0 
";

    close OUT;

    print "$pbsFile created \n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all ParClip::PARalyzer tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my @resultFiles = ();
    push( @resultFiles, "${resultDir}/${sampleName}/${sampleName}.cluster.csv" );
    push( @resultFiles, "${resultDir}/${sampleName}/${sampleName}.distribution.csv" );
    push( @resultFiles, "${resultDir}/${sampleName}/${sampleName}.groups.csv" );
    push( @resultFiles, "${resultDir}/${sampleName}/${sampleName}.target.csv" );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
