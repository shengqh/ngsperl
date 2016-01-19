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
  $self->{_name}   = __PACKAGE__;
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
  my $extension = ( ( $mmconfig->{OutputSuffix} ) ? $mmconfig->{OutputSuffix} : "" ) . ( ( $mmconfig->{OutputFormat} eq "pepXML" ) ? ".pepXML" : ".mzid" );
  return $extension;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );
  
  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );
  
  my $database   = get_param_file( $config->{$section}{database},   "database",   1 );
  my $cfgfile    = get_param_file( $config->{$section}{cfgfile},    "cfgfile",    1 );

  my $extension = $self->getExtension($cfgfile);

  my %mgffiles = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %mgffiles ) {
    my @sample_files = @{ $mgffiles{$sample_name} };
    my $pbs_file     = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name     = basename($pbs_file);
    my $log         = $self->get_log_filename( $log_dir, $sample_name );
    my $log_desc    = $cluster->get_log_description($log);
    open( my $out, ">$pbs_file" ) or die $!;

    print $out "$pbs_desc
$log_desc

$path_file

cd $result_dir 
";

    for my $sampleFile (@sample_files) {
      my $sname = basename($sampleFile);
      my $result_file = change_extension( $sname, $extension );
      print $out "if [ ! -s $result_file ]; then
  myrimatch -cfg $cfgfile -workdir $result_dir -cpus $thread -ProteinDatabase $database $sampleFile
fi

";
    }

    print $out "
echo finished=`date` 

exit 0
";
    close $out;

    print $sh "\$MYCMD ./$pbs_name \n";
    print "$pbs_file created\n";
  }
  print $sh "exit 0\n";
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );
  my $cfgfile = get_param_file( $config->{$section}{cfgfile}, "cfgfile", 1 );
  my $extension = $self->getExtension($cfgfile);

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    my @sample_files = @{ $raw_files{$sample_name} };
    for my $sampleFile (@sample_files) {
      my $sname = basename($sampleFile);
      my $result_file = change_extension( $sname, $extension );
      push( @result_files, "${result_dir}/${result_file}" );
    }

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
