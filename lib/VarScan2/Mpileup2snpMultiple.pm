#!/usr/bin/perl
package VarScan2::Mpileup2snpMultiple;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::UniqueTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_vs2";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $varscan2_jar = get_param_file( $config->{$section}{VarScan2_jar}, "VarScan2_jar", 1 );
  my $faFile       = get_param_file( $config->{$section}{fasta_file},   "fasta_file",   1 );

  my $mpileup_options = get_option( $config, $section, "mpileup_options", "" );
  my $java_option     = get_option( $config, $section, "java_option",     "" );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );

  my $log_desc   = $cluster->get_log_description($log);
  my $final_file = "${task_name}.snp.vcf";
  my $pbs        = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  my $sampleNameFile = $result_dir . "/sampleNames.txt";
  open( my $snf, ">$sampleNameFile" ) or die "Cannot create $sampleNameFile";
  my $samples = [];
  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    push( @$samples, $sample_files[0] );
    print $snf $sample_name . "\n";
  }
  close($snf);

  my $sample_list = join( ' ', @$samples);

  print $pbs "
samtools mpileup $mpileup_options -f $faFile $sample_list | java $java_option -jar $varscan2_jar mpileup2snp $option --vcf-sample-list sampleNames.txt --output-vcf 1 > $final_file
";
  $self->close_pbs( $pbs, $pbs_file );

  print "!!!pbs file $pbs_file created.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result       = {};
  my @result_files = ();
  my $final_file   = "${task_name}.snp.vcf";
  push( @result_files, "$result_dir/${final_file}" );
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
