#!/usr/bin/perl
package Alignment::STARIndex;

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
  $self->{_name}   = "STARIndex";
  $self->{_suffix} = "_si";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my %sjdbFiles = %{ get_raw_files( $config, $section ) };

  my $transcript_gtf = get_param_file( $config->{$section}{transcript_gtf}, "transcript_gtf", 0 );

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my $pbsFile = $self->pbsfile( $pbsDir, $task_name );
  my $pbsName = basename($pbsFile);
  my $log     = $self->logfile( $logDir, $task_name );

  my $log_desc = $cluster->get_log_desc($log);

  my $final = $resultDir . "/" . $task_name . ".tab";

  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

if [ -s $final ]; then
  echo job has already been done. if you want to do again, delete ${resultDir}/${final} and submit job again.
  exit 0;
fi

echo STARIndex=`date` 
";

  for my $sampleName ( sort keys %sjdbFiles ) {
    my @sjdbs = @{ $sjdbFiles{$sampleName} };
    for my $sjdb (@sjdbs) {
      print OUT "awk 'BEGIN {OFS=\"\t\"; strChar[0]=\".\"; strChar[1]=\"+\"; strChar[2]=\"-\";} {if(\$5>0){print \$1,\$2,\$3,strChar[\$4]}}' $sjdb >> $final";
    }
  }

  print OUT "STAR $option --runThreadN $thread --runMode genomeGenerate --genomeDir . --genomeFastaFiles $faFile --sjdbGTFfile $transcript_gtf --sjdbFileChrStartEnd $final \n";

  print OUT "echo finished=`date`

exit 0
";

  close(OUT);

  print "$pbsFile created\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $result = { $task_name => [ $resultDir . "/merged.gtf" ] };

  return $result;
}

1;
