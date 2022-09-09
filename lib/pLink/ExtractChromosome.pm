#!/usr/bin/perl
package pLink::ExtractChromosome;

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
use Pipeline::PipelineUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_ec";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $plinkFiles = get_raw_files( $config, $section );
  my $chromosome = get_option( $config, $section, "chromosome" );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  my $expectFiles = $self->result( $config, $section );
  for my $sampleName ( keys %$plinkFiles ) {
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sampleName );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sampleName );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";

    my $final_file = $sampleName . "." . $chromosome . ".bed";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
    my $bedPrefix = $plinkFiles->{$sampleName};
    if( is_array($bedPrefix) ){
      $bedPrefix = $bedPrefix->[0];
    }

    print $pbs "awk '{ if (\$1 == $chromosome) print \$2 }' ${bedPrefix}.bim > ${sampleName}.${chromosome}.txt 
plink --bfile $bedPrefix --extract ${sampleName}.${chromosome}.txt --make-bed --out ${sampleName}.${chromosome} 
if [ -s $final_file ]; then
  rm  ${sampleName}.${chromosome}.txt ${sampleName}.${chromosome}.nosex ${sampleName}.${chromosome}.log
fi
";
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section, 0 );

  my $plinkFiles = get_raw_files( $config, $section );
  my $chromosome = get_option( $config, $section, "chromosome" );

  my $result = {};
  for my $sampleName ( keys %$plinkFiles ) {
    my @result_files = ();
    push( @result_files, "${result_dir}/${sampleName}.${chromosome}.bed" );
    push( @result_files, "${result_dir}/${sampleName}.${chromosome}.bim" );
    $result->{$sampleName} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
