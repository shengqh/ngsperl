#!/usr/bin/perl
package CQS::TCGA;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(tcga_download tcga_get_coordinate)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub output_header {
  my ( $pbsFile, $pbsDesc, $path_file, $log ) = @_;

  #print "writing file " . $pbsFile . "\n";
  open( OUT, ">$pbsFile" ) or die $!;
  print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file
";
}

sub output_footer() {
  print OUT "echo finished=`date`\n";
  close OUT;

  #print "close file \n";
}

my $bamfilter = sub {
  my $filename = shift;

  return ( $filename =~ /.bam$/ );
};

sub tcga_download {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

  my $idfile = get_param_file( $config->{$section}{idfile}, "analysis id file", 1 );

  my $batchindex = $config->{$section}{batchidindex};
  if ( !defined $batchindex ) {
    die "Define batchidindex at section $section";
  }

  my $batches = $config->{$section}{batches} or die "Define batches at section $section";
  my @batches = @{$batches};

  #print @batches;

  my %batchmap = map { $_ => 1 } @batches;

  my $tcgaidindex = $config->{$section}{tcgaidindex};
  if ( !defined $tcgaidindex ) {
    die "Define tcgaidindex at section $section";
  }

  my $analysisidindex = $config->{$section}{analysisidindex};
  if ( !defined $analysisidindex ) {
    die "Define analysisidindex at section $section";
  }

  open( DAT, $idfile ) || die("Could not open file $idfile!");
  my $line     = <DAT>;
  my @raw_data = <DAT>;
  close(DAT);

  #print @raw_data;

  #print %batchmap;

  my $rawdir = create_directory_or_die( $resultDir . "/raw" );

  my $shfile = $pbsDir . "/${task_name}.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  my $index  = 0;
  my $dindex = 0;
  foreach $line (@raw_data) {
    chomp($line);
    my @parts = split( '\t', $line );

    my $batch = $parts[$batchindex];

    #print $batch . "\n";
    if ( scalar(@batches) > 0 ) {
      if ( !exists( $batchmap{$batch} ) ) {
        next;
      }
    }

    my $partSize = @parts;

    my $tcga       = $parts[$tcgaidindex];
    my $analysisid = $parts[$analysisidindex];
    my $url        = "https://cghub.ucsc.edu/cghub/data/analysis/download/" . $analysisid;

    if ( 0 == $index % 10 ) {
      if ( $dindex != 0 ) {
        output_footer();
      }
      $dindex = $dindex + 1;
      my $pbsName = "${task_name}_${dindex}_download.pbs";
      my $pbsFile = "${pbsDir}/$pbsName";
      my $log     = "${logDir}/${task_name}_${dindex}_download.log";

      print SH "\$MYCMD ./$pbsName \n";

      output_header( $pbsFile, $pbsDesc, $path_file, $log );

      print OUT "echo download=`date` \n";
      print OUT "cd $rawdir \n";

      print "$pbsFile created\n";
    }

    $index = $index + 1;
    print OUT "echo $tcga `date` \n";
    print OUT "GeneTorrent -v -c ~/.ssh/mykey.pem -C ~/pylibs/share/GeneTorrent -d $url \n";
  }
  if ( $dindex != 0 ) {
    output_footer();
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }
}

sub tcga_get_coordinate {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option ) = get_parameter( $config, $section );

  my $idfile = get_param_file( $config->{$section}{idfile}, "analysis id file", 1 );

  my $batchindex = $config->{$section}{batchidindex};
  if ( !defined $batchindex ) {
    die "Define batchidindex at section $section";
  }

  my $batches = $config->{$section}{batches} or die "Define batches at section $section";
  my @batches = @{$batches};

  my %batchmap = map { $_ => 1 } @batches;

  my $tcgaidindex = $config->{$section}{tcgaidindex};
  if ( !defined $tcgaidindex ) {
    die "Define tcgaidindex at section $section";
  }

  my $analysisidindex = $config->{$section}{analysisidindex};
  if ( !defined $analysisidindex ) {
    die "Define analysisidindex at section $section";
  }

  my $coordinateindex = $config->{$section}{coordinateindex};
  if ( !defined $coordinateindex ) {
    die "Define coordinateindex at section $section";
  }

  open( DAT, $idfile ) || die("Could not open file $idfile!");
  my $line     = <DAT>;
  my @raw_data = <DAT>;
  close(DAT);

  my $pbsFile = $pbsDir . "/${task_name}_coordidate.pbs";
  my $log     = $logDir . "/${task_name}_coordidate.log";

  output_header( $pbsFile, $pbsDesc, $path_file, $log );

  my $rawdir        = create_directory_or_die( $resultDir . "/raw" );
  my $coordinatedir = create_directory_or_die( $resultDir . "/coordindates" );

  print OUT "echo download=`date` \n";

  foreach $line (@raw_data) {
    chomp($line);
    my @parts = split( '\t', $line );

    my $batch = $parts[$batchindex];
    if ( scalar(@batches) > 0 ) {
      if ( !exists( $batchmap{$batch} ) ) {
        next;
      }
    }

    my $partSize   = @parts;
    my $tcga       = $parts[$tcgaidindex];
    my $analysisid = $parts[$analysisidindex];
    my $coordinate = $parts[$coordinateindex];

    if ( !defined($analysisid) ) {
      next;
    }

    if ( length($analysisid) == 0 ) {
      next;
    }

    if ( !defined($coordinate) ) {
      next;
    }

    if ( length($coordinate) == 0 ) {
      next;
    }

    #remove the character '"' and ','
    $coordinate =~ s/[", ]//g;
    if ( $coordinate =~ /(.+):(.+)-(.+)$/ ) {
      if ( $2 > $3 ) {
        print "$coordinate => ";
        $coordinate = "$1:$3-$2";
        print "$coordinate \n";
      }
    }

    my $subdir = $rawdir . '/' . $analysisid;

    my $coord = $coordinate;
    $coord =~ s/:/./g;
    my $targetfile = "${coordinatedir}/${tcga}.${coord}.sam";

    my @bamfiles     = list_files( $subdir, $bamfilter );
    my $bamfile      = $subdir . "/" . $bamfiles[0];
    my $bamfileindex = $bamfile . ".bai";

    if ( !-e $bamfileindex ) {
      my $cmd = "samtools index " . $bamfile . " ";
      print OUT $cmd . "\n";
    }

    my $cmd = "samtools view " . $bamfile . " " . $coordinate . " > " . $targetfile . " ";
    print OUT $cmd . "\n";
  }
  output_footer();
  print "$pbsFile created\n";
}

1;
