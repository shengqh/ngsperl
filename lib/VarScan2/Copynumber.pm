#!/usr/bin/perl
package VarScan2::Copynumber;

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
  $self->{_name}   = "VarScan2::Copynumber";
  $self->{_suffix} = "_cp";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $varscan2_jar = get_param_file( $config->{$section}{VarScan2_jar}, "VarScan2_jar", 1 );
  my $faFile       = get_param_file( $config->{$section}{fasta_file},   "fasta_file",   1 );

  my $mpileup_options = get_option( $config, $section, "mpileup_options", "" );
  my $call_options =  get_option( $config, $section, "call_options", "--amp-threshold 0.3 --del-threshold 0.3" );

  my %group_sample_map = %{ $self->get_group_sample_map( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";
  print SH "cd $pbsDir\n";

  my $java_option = get_option( $config, $section, "java_option", "" );

  for my $groupName ( sort keys %group_sample_map ) {
    my @sampleFiles = @{ $group_sample_map{$groupName} };
    my $sampleCount = scalar(@sampleFiles);

    if ( $sampleCount != 2 ) {
      die "SampleFile should be normal,tumor paired.";
    }

    my $curDir = create_directory_or_die( $resultDir . "/$groupName" );

    my $normal = $sampleFiles[0][1];
    my $tumor  = $sampleFiles[1][1];

    my $cpRawFile   = "${groupName}.copynumber";
    my $cpCallFile = "${groupName}.call";

    my $pbsFile = $self->pbsfile( $pbsDir, $groupName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $groupName );

    print SH "\$MYCMD ./$pbsName \n";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file 

echo varscan2_somatic_copynumber=`date` 

cd $curDir

if [ ! -s $cpRawFile ]; then
  if [ ! -s ${normal}.bai ]; then
    samtools index ${normal}
  fi

  if [ ! -s ${tumor}.bai ]; then
    samtools index ${tumor}
  fi

  samtools mpileup -q 1 -f $faFile $normal $tumor | awk 'NF==9 && \$4!=0' | java $java_option -jar $varscan2_jar copynumber $groupName $option 
fi

if [[ -s $cpRawFile && ! -s $cpCallFile ]]; then
  java $java_option -jar $varscan2_jar copyCaller $cpRawFile $call_options --output-file $cpCallFile 
fi

echo finished=`date`
";
    close OUT;

    print "$pbsFile created \n";
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all VarScan2 tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $groupName ( keys %{$groups} ) {
    my @resultFiles = ();
    my $curDir      = $resultDir . "/$groupName";
    push( @resultFiles, "$curDir/${groupName}.call" );
    $result->{$groupName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
