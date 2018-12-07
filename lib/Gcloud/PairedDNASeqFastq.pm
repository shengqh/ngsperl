#!/usr/bin/perl
package Gcloud::PairedDNASeqFastq;

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
use Alignment::AlignmentUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_gpd";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );

  my $pipeline_file = get_option_file( $config, $section, "pipeline_file" );
  my $wdl_file = get_option_file( $config, $section, "wdl_file");
  my $input_json_file = get_option_file( $config, $section, "input_json_file" );
  my $input_option_file = get_option_file( $config, $section, "input_option_file" );
  my $result_bucket = get_option( $config, $section, "result_bucket");

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $fastq1 = $sample_files[0];
    my $fastq2 = $sample_files[1];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $gslog = $result_bucket . "/logging/";
    my $gsworkspace=$result_bucket . "/workspace/";
    my $gsresult = $result_bucket . "/result/";
    
    my $sample_input_file = $self->get_file( $result_dir, $sample_name, ".inputs.json" );
    
    my $fh;
    my $fo;
    my $row;
    open($fh, '<', $input_json_file);
    open($fo, '>', $sample_input_file);
    while ($row = <$fh>) {
      if ($row =~ /dnaseq.sample_name/){
        print $fo "    \"dnaseq.sample_name\": \"$sample_name\",\n";
      }elsif($row =~ /dnaseq.base_file_name/){
        print $fo "    \"dnaseq.base_file_name\": \"$sample_name\",\n";
      }elsif($row =~ /dnaseq.final_gvcf_base_name/){
        print $fo "    \"dnaseq.final_gvcf_base_name\": \"$sample_name\",\n";
      }elsif($row =~ /dnaseq.raw_fastq1/){
        print $fo "    \"dnaseq.raw_fastq1\": \"$fastq1\",\n";
      }elsif($row =~ /dnaseq.raw_fastq2/){
        print $fo "    \"dnaseq.raw_fastq2\": \"$fastq2\",\n";
      }elsif($row =~ /dnaseq.bwa_commandline/){
        print $fo "    \"dnaseq.bwa_commandline\": \"bwa mem -K 100000000 -p -v 3 -t 16 -R '\@RG\\tID:" . $sample_name ."\\tPU:illumina\\tLB:" . $sample_name . "\\tSM:" . $sample_name ."\\tPL:illumina' -Y \$bash_ref_fasta\",\n";
      }else{
        print $fo $row;
      }
    }
    close($fh);
    close($fo);
    
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

    print $pbs "
gcloud alpha genomics pipelines run --pipeline-file $pipeline_file \\
  --logging $gslog \\
  --inputs-from-file WDL=$wdl_file \\
  --inputs-from-file WORKFLOW_INPUTS=$sample_input_file \\
  --inputs-from-file WORKFLOW_OPTIONS=$input_option_file \\
  --inputs WORKSPACE=$gsworkspace --inputs OUTPUTS=$gsresult $option 2>&1 | tee ${pbs_file}.id
    
";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  
  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;