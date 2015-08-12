#!/usr/bin/perl
package CQS::Perl;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use File::Spec;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "Perl";
  $self->{_suffix} = "_perl";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $perlFile = get_option_value( $config->{$section}{perlFile}, "perlFile", 1 );
  my $is_absolute = File::Spec->file_name_is_absolute($perlFile);
  if ( !$is_absolute ) {
    $perlFile = dirname(__FILE__) . "/$perlFile";
  }
  if ( !( -e $perlFile ) ) {
    die("perlFile $perlFile defined but not exists!");
  }

  my $output_ext = get_option( $config, $section, "output_ext", 0 );

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  my %parameterFiles2 = ();
  if ( has_raw_files( $config, $section, "source2" ) ) {
    %parameterFiles2 = %{ get_raw_files( $config, $section, "source2" ) };
  }
  my %parameterFiles3 = ();
  if ( has_raw_files( $config, $section, "source3" ) ) {
    %parameterFiles3 = %{ get_raw_files( $config, $section, "source3" ) };
  }
  my %parameterFiles4 = ();
  if ( has_raw_files( $config, $section, "source4" ) ) {
    %parameterFiles4 = %{ get_raw_files( $config, $section, "source4" ) };
  }
  my %parameterFiles5 = ();
  if ( has_raw_files( $config, $section, "source5" ) ) {
    %parameterFiles5 = %{ get_raw_files( $config, $section, "source5" ) };
  }

  my $shfile = $pbsDir . "/${task_name}_perl.sh";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  if ($sh_direct) {
    print SH "export MYCMD=\"bash\" \n";
  }
  else {
    print SH "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
  }

  for my $sampleName ( sort keys %rawFiles ) {
    my @sampleFiles = @{ $rawFiles{$sampleName} };
    my $sampleCount = scalar(@sampleFiles);
    my $samples     = join( ' ', @sampleFiles );

    my $parameterFile2 = "";
    if ( defined $parameterFiles2{$sampleName} ) {
      my @files = @{ $parameterFiles2{$sampleName} };
      my $file  = $files[0];
      $parameterFile2 = "$file";
    }

    my $parameterFile3 = "";
    if ( defined $parameterFiles3{$sampleName} ) {
      my @files = @{ $parameterFiles3{$sampleName} };
      my $file  = $files[0];
      $parameterFile3 = "$file";
    }

    my $parameterFile4 = "";
    if ( defined $parameterFiles4{$sampleName} ) {
      my @files = @{ $parameterFiles4{$sampleName} };
      my $file  = $files[0];
      $parameterFile4 = "$file";
    }

    my $parameterFile5 = "";
    if ( defined $parameterFiles5{$sampleName} ) {
      my @files = @{ $parameterFiles5{$sampleName} };
      my $file  = $files[0];
      $parameterFile5 = "$file";
    }

    #    my $curDir      = create_directory_or_die( $resultDir . "/$sampleName" );

    #    my $pbsName = "${sampleName}_perl.pbs";
    #    my $pbsFile = "${pbsDir}/$pbsName";

    my $pbsFile = $self->pbsfile( $pbsDir, $sampleName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $sampleName );

    print SH "\$MYCMD ./$pbsName \n";

    #    my $log = "${logDir}/${sampleName}_perl.log";

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

echo Perl=`date`
cd $target_dir/result
perl $perlFile $sampleName$output_ext $option $samples $parameterFile2 $parameterFile3 $parameterFile4 $parameterFile5

echo finished=`date`
";
    close OUT;

    print "$pbsFile created \n";
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all Perl tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %rawFiles = %{ get_raw_files( $config, $section ) };
  my $output_ext = get_option( $config, $section, "output_ext", 0 );

  my $result = {};
  for my $sampleName ( keys %rawFiles ) {

    #    my @sampleFiles = @{ $rawFiles{$sampleName} };
    #    my @resultFiles = ();
    #    for my $sampleFile (@sampleFiles) {
    #      my $name = basename($sampleFile);
    #		push( @resultFiles, "${resultDir}/${name}${output_ext}" );
    #    }
    my @resultFiles = ();
    push( @resultFiles, "${resultDir}/${sampleName}${output_ext}" );

    $result->{$sampleName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
