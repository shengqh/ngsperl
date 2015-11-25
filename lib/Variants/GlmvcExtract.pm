#!/usr/bin/perl
package Variants::GlmvcExtract;

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
use Data::Dumper;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "GlmvcExtract";
  $self->{_suffix} = "_ge";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $glmvcfile = get_param_file( $config->{$section}{execute_file}, "execute_file", 1 );

  my $rawFiles = get_raw_files( $config, $section );
  my $bamFiles = get_raw_files( $config, $section, "bam_files" );

  print Dumper($rawFiles);
  print Dumper($bamFiles);

  my %group_sample_map = ();
  my %group_name_map   = ();

  my $fafile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my $groups = get_raw_files( $config, $section, "groups" );
  for my $groupName ( sort keys %{$rawFiles} ) {
    my @samples = @{ $groups->{$groupName} };
    my @gfiles  = ();
    my @names   = ();
    my $index   = 0;
    foreach my $sampleName (@samples) {
      my @sampleBamFiles = @{ $bamFiles->{$sampleName} };
      push( @gfiles, $sampleBamFiles[0] );
      push( @names,  $sampleName );
    }
    $group_sample_map{$groupName} = \@gfiles;
    $group_name_map{$groupName}   = \@names;
  }

  print Dumper(%group_name_map);
  print Dumper(%group_sample_map);

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $groupName ( sort keys %{$rawFiles} ) {
    my $validateFile = $rawFiles->{$groupName}[0];
    my @sampleFiles  = @{ $group_sample_map{$groupName} };
    my @sampleNames  = @{ $group_name_map{$groupName} };

    my $samples = join( ',', @sampleFiles );
    my $names   = join( ',', @sampleNames );

    my $curDir = create_directory_or_die( $resultDir . "/$groupName" );

    my $pbsFile = $self->pbsfile( $pbsDir, $groupName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $groupName );

    print SH "\$MYCMD ./$pbsName \n";

    my $log_desc = $cluster->get_log_desc($log);
    my $final    = "${groupName}.tsv";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file 

echo Glmvc=`date` 

cd $curDir

if [ -s $final ]; then
  echo job has already been done. if you want to do again, delete ${curDir}/${final} and submit job again.
  exit 0;
fi      
      
mono-sgen $glmvcfile extract $option --bam_files $samples --bam_names $names -f $fafile -o ${curDir}/$final -v $validateFile

echo finished=`date`
";

    close OUT;

    print "$pbsFile created \n";

  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $rawFiles = get_raw_files( $config, $section );

  my $result = {};
  for my $groupName ( keys %{$rawFiles} ) {
    my @resultFiles = ();
    my $curDir      = $resultDir . "/$groupName";

    push( @resultFiles, "$curDir/${groupName}.tsv" );
    $result->{$groupName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
