#!/usr/bin/perl
package Variants::GlmvcAnnotation;

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
  $self->{_name}   = "GlmvcAnnotation";
  $self->{_suffix} = "_ga";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $glmvcfile = get_param_file( $config->{$section}{execute_file}, "execute_file", 1 );
  my $source_type = $config->{$section}{source_type} or die "source_type is not defined in $section";

  my $rnaediting_db = get_directory( $config, $section, "rnaediting_db", 0 );
  if ( defined $rnaediting_db ) {
    $option = $option . " --rnaediting_db $rnaediting_db ";
  }

  my $annovar_buildver = $config->{$section}{annovar_buildver};
  if ( defined $annovar_buildver ) {
    $option = $option . " --annovar_buildver $annovar_buildver ";

    my $annovar_protocol = $config->{$section}{annovar_protocol};
    if ( defined $annovar_protocol ) {
      $option = $option . " --annovar_protocol $annovar_protocol ";
    }

    my $annovar_operation = $config->{$section}{annovar_operation};
    if ( defined $annovar_operation ) {
      $option = $option . " --annovar_operation $annovar_operation ";
    }

    my $annovar_db = $config->{$section}{annovar_db};
    if ( defined $annovar_db ) {
      $option = $option . " --annovar_db $annovar_db ";
    }
  }

  my $distance_exon_gtf = get_param_file( $config->{$section}{distance_exon_gtf}, "distance_exon_gtf", 0 );
  if ( defined $distance_exon_gtf ) {
    $option = $option . " --distance_exon_gtf $distance_exon_gtf ";
  }

  my $anno = defined $rnaediting_db || defined $annovar_buildver || defined $annovar_buildver;

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";
  print SH "cd $pbsDir \n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $sampleFile  = $sampleFiles[0];
    my $curDir      = create_directory_or_die( $resultDir . "/$sampleName" );

    my $pbsFile  = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName  = basename($pbsFile);
    my $log      = $self->logfile( $logDir, $sampleName );
    my $log_desc = $cluster->get_log_desc($log);
    my $final    = "${sampleName}.annotation.tsv";

    print SH "\$MYCMD ./$pbsName \n";

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
      
mono $glmvcfile annotation $option -i $sampleFile -o ${final}

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
  for my $sampleName ( keys %{$rawFiles} ) {
    my @resultFiles = ();
    my $curDir      = $resultDir . "/$sampleName";
    push( @resultFiles, "$curDir/${sampleName}.annotation.tsv" );
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}
1;
