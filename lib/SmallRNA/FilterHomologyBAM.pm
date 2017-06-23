#!/usr/bin/perl
package SmallRNA::FilterHomologyBAM;

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
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_fhm";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $py_script = dirname(__FILE__) . "/filterHomologyBAM.py";
  if ( !-e $py_script ) {
    die "File not found : " . $py_script;
  }

  my $referencePrefix    = get_option( $config, $section, "reference_prefix" );
  my $homologyPrefix     = get_option( $config, $section, "homology_prefix" );
  my $outputReferneceBAM = get_option( $config, $section, "output_reference_bam", 0 );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sampleFile   = $sample_files[0];
    my $final_file   = $sample_name . "." . $homologyPrefix . ".bam";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
    print $pbs "python $py_script -i $sampleFile -o $final_file $option --referencePrefix $referencePrefix --homologyPrefix $homologyPrefix \n";
    
    if($outputReferneceBAM){
      my $py_script2 = dirname(__FILE__) . "/filterReferenceBAM.py";
      if ( !-e $py_script2 ) {
        die "File not found : " . $py_script2;
      }
      my $refBam = $sample_name . "." . $referencePrefix . ".bam";
      print $pbs "python $py_script2 -i $sampleFile -o $refBam $option --referencePrefix $referencePrefix \n";
    }
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $referencePrefix    = get_option( $config, $section, "reference_prefix" );
  my $homologyPrefix     = get_option( $config, $section, "homology_prefix" );
  my $outputReferneceBAM = get_option( $config, $section, "output_reference_bam", 0 );

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my $final_file = $result_dir . "/" . $sample_name . "." . $homologyPrefix . ".bam";

    my @result_files = ();
    push( @result_files, $final_file );

    if ($outputReferneceBAM) {
      push( @result_files, $result_dir . "/" . $sample_name . "." . $referencePrefix . ".bam" );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
