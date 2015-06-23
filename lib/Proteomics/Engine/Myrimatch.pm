#!/usr/bin/perl
package Proteomics::Engine::Myrimatch;

use strict;
use warnings;
use File::Basename;
use String::Util qw(trim);
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
  $self->{_name}   = "Proteomics::Engine::Myrimatch";
  $self->{_suffix} = "_mm";
  bless $self, $class;
  return $self;
}

sub readConfigurationFile {
  my ( $self, $cfgfile ) = @_;

  my $result = {};
  open my $fh, '<', $cfgfile or die "Unable to open file:$!\n";
  while ( my $row = <$fh> ) {
    chomp $row;
    my @parts = split( '=', $row );
    if ( scalar(@parts) > 1 ) {
      my $key   = trim( $parts[0] );
      my $value = trim( $parts[1] );
      $value =~ s/^"(.*)"$/$1/;
      $result->{$key} = $value;
    }
  }

  close $fh;
  return ($result);
}

sub getExtension {
  my ( $self, $cfgfile ) = @_;
  my $mmconfig = $self->readConfigurationFile($cfgfile);
  my $extension = ( ( $mmconfig->{OutputSuffix} ) ? $mmconfig->{OutputSuffix} : "" ) . ( ( $mmconfig->{OutputFormat} eq "pepXML" ) ? ".pepXml" : ".mzid" );
  return $extension;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );
  my $database   = get_param_file( $config->{$section}{database},   "database",   1 );
  my $executable = get_param_file( $config->{$section}{executable}, "executable", 1 );
  my $cfgfile    = get_param_file( $config->{$section}{cfgfile},    "cfgfile",    1 );

  my $extension = $self->getExtension($cfgfile);

  my %mgffiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct);

  for my $sampleName ( sort keys %mgffiles ) {
    my @sampleFiles = @{ $mgffiles{$sampleName} };
    my $pbsFile     = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName     = basename($pbsFile);
    my $log         = $self->logfile( $logDir, $sampleName );
    my $log_desc    = $cluster->get_log_desc($log);
    open( OUT, ">$pbsFile" ) or die $!;

    print OUT "$pbsDesc
$log_desc

$path_file

cd $resultDir 
";

    for my $sampleFile (@sampleFiles) {
      my $sname = basename($sampleFile);
      my $resultFile = change_extension( $sname, $extension );
      print OUT "if [ ! -s $resultFile ]; then
  myrimatch -cfg $cfgfile -workdir $resultDir -cpus $thread -ProteinDatabase $database $sampleFile
fi
";
    }

    print OUT "
echo finished=`date` 

exit 0
";
    close(OUT);

    print SH "\$MYCMD ./$pbsName \n";
    print "$pbsFile created\n";
  }
  print SH "exit 0\n";
  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );
  my $cfgfile = get_param_file( $config->{$section}{cfgfile}, "cfgfile", 1 );
  my $extension = $self->getExtension($cfgfile);

  my %rawFiles = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {
    my @resultFiles = ();
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    for my $sampleFile (@sampleFiles) {
      my $sname = basename($sampleFile);
      my $resultFile = change_extension( $sname, $extension );
      push( @resultFiles, "${resultDir}/${resultFile}" );
    }

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
