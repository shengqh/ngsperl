#!/usr/bin/perl
package Proteomics::Engine::Comet;

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
  $self->{_name}   = "Proteomics::Engine::Comet";
  $self->{_suffix} = "_comet";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;
  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $param_file = get_param_file( $config->{$section}{param_file}, "param_file", 1 );
  my $database   = get_param_file( $config->{$section}{database},   "database",   1 );
  my $delete_temp_ms2 = get_option($config, $section, "delete_temp_ms2", 1); 

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir

";
    for my $sampleFile (@sampleFiles) {
      my $sname = basename($sampleFile);
      my $resultFile = change_extension( $sname, ".pep.xml" );

      print OUT "if [ ! -s $resultFile ]; then\n";

      my $ismgf = $sname =~ /\.mgf$/i;
      my $tempFile = $resultDir . "/" . change_extension( $sname, ".ms2" );

      if ($ismgf) {
        my $proteomicstools = get_param_file( $config->{$section}{proteomicstools}, "proteomicstools", 1 );
        my $titleformat = get_option( $config, $section, "titleformat" );
        print OUT "  if [ ! -s $tempFile ]; then
    mono $proteomicstools MGF2MS2 -i $sampleFile -t $titleformat -o $tempFile
  fi
";
        $sampleFile = $tempFile;
      }

      print OUT "  comet -P$param_file -D$database $sampleFile
  if [ -s $resultFile ]; then
    RefreshParser $resultFile $database
  fi
";

      if ($ismgf && $delete_temp_ms2) {
        print OUT "  if [ -s $resultFile ]; then
    rm $tempFile
  fi
";
      }

      print OUT "
fi
";
    }
    print OUT "
echo finished=`date`

exit 0 
";
    close OUT;

    print "$pbsFile created \n";
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all ", $self->{_name}, " tasks.\n";

  #`qsub $pbsFile`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my @resultFiles = ();

    for my $sampleFile (@sampleFiles) {
      my $sname = basename($sampleFile);
      my $resultFile = change_extension( $sname, ".pep.xml" );
      push( @resultFiles, "${resultDir}/${resultFile}" );
    }
    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
