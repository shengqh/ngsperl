#!/usr/bin/perl
package Proteomics::Engine::Msgf2Pin;

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
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_m2pin";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $decoyPattern = get_option( $config, $section, "decoy_pattern" );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $final_file = $sample_name . ".msgf.pin";
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    if ( scalar(@sample_files) == 1 ) {
      my $sampleFile = $sample_files[0];
      print $pbs "
echo processing $sampleFile
msgf2pin -P $decoyPattern -o $final_file $sampleFile
";
    }
    else {
      my $temp = $final_file . ".tmp";
      print $pbs "
if [ -e $temp ]; then
  rm $temp
fi
";
      my $index = 0;
      for my $index ( 0 .. $#sample_files ) {
        my $curtemp    = $temp . "." . $index;
        my $sampleFile = $sample_files[$index];
        if ( $index == 0 ) {
          print $pbs "
echo processing $sampleFile
msgf2pin -P $decoyPattern -o $temp $sampleFile
";
        }else{
          print $pbs "
echo processing $sampleFile
msgf2pin -P $decoyPattern -o $curtemp $sampleFile
tail -n +3 $curtemp >> $temp
rm $curtemp
";
          
        }
      }
      
      print $pbs "
mv $temp $final_file
"
    }

    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all ", $self->{_name}, " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $final_file = $sample_name . ".msgf.pin";
    my @result_files = ();
    push( @result_files, "${result_dir}/${final_file}" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
