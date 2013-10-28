#!/usr/bin/perl
package CQS::Samtools;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::ConfigUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(mpileup)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub mpileup {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

  my $fafile = get_param_file( $config->{$section}{reference_sequence}, "reference_sequence", 0 );
  my $mincount = $config->{$section}{minimum_count};

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbsDir . "/${task_name}.submit";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $pbsName = "${sampleName}_mpileup.pbs";
    my $pbsFile = "${pbsDir}/$pbsName";

    print SH "\$MYCMD ./$pbsName \n";

    my $log = "${logDir}/${sampleName}_mpileup.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file 

cd $resultDir 

echo mpileup=`date` 
";

    my $sampleCount = scalar(@sampleFiles);

    for my $sampleFile (@sampleFiles) {
      my $bamindex = $sampleFile . ".bai";

      print OUT "if [ ! -s $bamindex ]; \n";
      print OUT "then \n";
      print OUT "  samtools index $sampleFile \n";
      print OUT "fi \n\n";
    }
    print OUT "samtools mpileup $option";
    if ( defined $fafile ) {
      print OUT " -f $fafile";
    }
    for my $sampleFile (@sampleFiles) {
      print OUT " $sampleFile";
    }

    if ( defined $mincount ) {
      if ( $mincount > 0 ) {
        print OUT " | awk '(\$4 > $mincount) && (\$7 > $mincount)' "
      }
    }
    print OUT "> ${sampleName}.mpileup \n\n";
    print OUT "echo finished=`date` \n";
    close OUT;

    print "$pbsFile created \n";
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all samtools mpileup tasks.\n";
}

1;
