#!/usr/bin/perl
package CQS::Shrimp2;

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
  $self->{_name} = "Shrimp2";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $shrimp2_index = $config->{$section}{shrimp2_index} or die "define ${section}::shrimp2_index first";
  die "shrimp2_index ${shrimp2_index}.genome not exist" if ( !-e "${shrimp2_index}.genome" );

  my $is_mirna = $config->{$section}{is_mirna} or die "define ${section}::is_mirna first";
  my $mirna = "-M mirna" if $is_mirna or "";

  my $output_bam = $config->{$section}{output_bam} or die "define ${section}::output_bam first";

  my %rawFiles = %{ get_raw_files( $config, $section, "source", ".fastq\$" ) };

  my $shfile = $pbsDir . "/${task_name}_shrimp2.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  if ($sh_direct) {
    print SH "export MYCMD=\"bash\" \n";
  }
  else {
    print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
  }

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{$rawFiles{$sampleName}};
    my $sampleFile =  $sampleFiles[0];

    my $shrimpFile = $sampleName . ".shrimp";
    my $samFile    = $sampleName . ".sam";
    my $bamFile    = $sampleName . ".bam";

    my $pbsName = "${sampleName}_srp2.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log    = "${logDir}/${sampleName}_srp2.log";
    my $curDir = create_directory_or_die( $resultDir . "/$sampleName" );

    open( OUT, ">$pbsFile" ) or die $!;

    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

echo shrimp2=`date`

cd $curDir

";

    if ($output_bam) {
      print OUT "gmapper -L $shrimp2_index $sampleFile $mirna $option --extra-sam-fields > $samFile

if [ -s $samFile ]; then
  samtools view -S -b $samFile | samtools sort - $sampleName
  samtools index $bamFile 
  samtools flagstat $bamFile > ${bamFile}.stat 
fi

echo finished=`date`
";
    }
    else {
      print OUT "gmapper -L $shrimp2_index $sampleFile $mirna $option --pretty >$shrimpFile
      
echo finished=`date` 
";
    }

    close OUT;

    print "$pbsFile created \n";
  }
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $output_bam = $config->{$section}{output_bam} or die "define ${section}::output_bam first";

  my %rawFiles = %{ get_raw_files( $config, $section, "source", ".fastq\$" ) };

  my $result = {};
  for my $sampleName ( sort keys %rawFiles ) {
    my $curDir      = $resultDir . "/$sampleName/";
    my @resultFiles = ();
    if ($output_bam) {
      push( @resultFiles, $curDir . $sampleName . ".bam" );
    }
    else {
      push( @resultFiles, $curDir . $sampleName . ".shrimp" );
    }

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
