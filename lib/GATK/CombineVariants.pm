#!/usr/bin/perl
package GATK::CombineVariants;

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

our @ISA = qw(CQS::UniqueTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_cv";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );
  my $gatk_jar = get_param_file( $config->{$section}{gatk_jar}, "gatk_jar", 1, not $self->using_docker() );

  #my $is_mutect = get_option($config, $section, "is_mutect", 0);

  my $filter_pass = get_option($config, $section, "filter_pass" , 0);

  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option || $java_option eq "" ) {
    $java_option = "-Xmx${memory}";
  }

  my %vcfFiles = %{ get_raw_files( $config, $section ) };

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $final_vcf  = $task_name . ".vcf";
  my $final_file  = $task_name . ".vcf.gz";
  my $merged_file = $task_name . ".tmp.vcf";

  my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  print $pbs "
java $java_option -jar $gatk_jar \\
  -T CombineVariants \\
  -R $faFile \\
";

  for my $sample_name ( sort keys %vcfFiles ) {
    my @sample_files = @{ $vcfFiles{$sample_name} };
    my $vcfFile      = $sample_files[0];
    print $pbs "  --variant $vcfFile \\\n";
  }

  my $additional_filter = $filter_pass?"|awk '{ if (\$7 == \"PASS\") { print } }'":"";

  print $pbs "  -o $merged_file \\
  -genotypeMergeOptions UNIQUIFY
  
status=\$?
if [[ \$status -ne 0 ]]; then
  touch $task_name.combined.failed
  rm $merged_file ${merged_file}.idx
else
  awk '{if(\$1 ~ /^##/){print;}else{exit;}}' $merged_file > $final_vcf
  awk '{if(\$1 ~ /^#/){print;}else{exit;}}' $merged_file | grep -v \"^##\" | sed -e \"s/.variant\\S*//g\" >> $final_vcf
  grep -v \"^#\" $merged_file $additional_filter | sed -e \"s/;set=filterInvariant\\S*//g\"  | sed -e \"s/;set=variant\\S*//g\" >> $final_vcf
  bgzip $final_vcf
  tabix $final_file
  rm -f $merged_file ${merged_file}.idx
  touch $task_name.succeed
fi 
";

  $self->close_pbs( $pbs, $pbs_file );
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my @result_files = ();
  my $final_file  = $task_name . ".vcf.gz";

  push( @result_files, $result_dir . "/" . $final_file );
  my $result = {};
  $result->{$task_name} = filter_array( \@result_files, $pattern );

  return $result;
}

1;
