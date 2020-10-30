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
use Gcloud::GcloudUtils;

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $project = get_option($config, $section, "project");
  my $pipeline_file = get_option_file( $config, $section, "pipeline_file" );
  my $wdl_file = get_option_file( $config, $section, "wdl_file");
  my $input_json_file = get_option_file( $config, $section, "input_json_file" );
  my $input_option_file = get_option_file( $config, $section, "input_option_file" );
  my $result_bucket = get_option( $config, $section, "result_bucket");
  my $calling_interval_file = get_option( $config, $section, "calling_interval_file", "" );

  my $gslog = $result_bucket . "/logging";
  my $gsworkspace=$result_bucket . "/workspace";
  my $gsresult = $result_bucket . "/result";
  my $gslogfailed = $gslog . "/failed/";
  my $gslogsucceed = $gslog . "/succeed/";

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  my $logcheckfile = $pbs_dir . "/" . $task_name . ".checklog.sh";
  open( my $flog, ">$logcheckfile" ) or die "Cannot create $logcheckfile";
  print $flog "gcloud config set project $project 

gslog=$gslog
gsworkspace=$gsworkspace
gsresult=$gsresult
gslogfailed=$gslogfailed
gslogsucceed=$gslogsucceed

arr=\$(gsutil ls \$gslog/*.log)
for logfile in \$arr
do
  sample_name=\$(echo \"\${logfile}\" | sed -e 's/.*\\///g' | sed -e 's/.log//g')
  log_file=\$gslog/\$sample_name.log
  failed_log_file=\$gslogfailed\$sample_name.log
  succeed_log_file=\$gslogsucceed\$sample_name.log
  final_file=\$gsresult/\$sample_name/\$sample_name.g.vcf.gz
  bb=\$(gsutil stat \$log_file)
  echo \"\"

  if [[ \$bb ]]; then
    aa=\$(gsutil stat \$final_file)
    echo \"\"

    if [[ \$aa ]]; then
      as=\$(gsutil stat \$succeed_log_file)
      echo \"\"
      if [[ ! \$as ]]; then
        echo gsutil mv \$log_file \$gslogsucceed
        gsutil mv \$log_file \$gslogsucceed
        continue
      fi
    fi

    for ifail in '' '.2' '.3' '.4' '.5' '.6'
    do
      afail=\$failed_log_file\$ifail
      ag=\$(gsutil stat \$afail)
      echo \"\"
      if [[ ! \$ag ]]; then
        break
      fi
    done

    echo gsutil mv \$log_file \$afail
    gsutil mv \$log_file \$afail
  fi
done

";
  close($flog);

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

    my $final_file = $gsresult . "/" . $sample_name . "/" . $sample_name . ".g.vcf.gz";
    my $sample_input_file = $self->get_file( $result_dir, $sample_name, ".inputs.json" );
    
    my $fh;
    my $fo;
    my $row;
    open($fh, '<', $input_json_file);
    open($fo, '>', $sample_input_file);
    while ($row = <$fh>) {
      if ($row =~ /PreProcessingForVariantDiscovery_GATK4.sample_name/){
        print $fo "    \"PreProcessingForVariantDiscovery_GATK4.sample_name\": \"$sample_name\",\n";
      }elsif($row =~ /PreProcessingForVariantDiscovery_GATK4.base_file_name/){
        print $fo "    \"PreProcessingForVariantDiscovery_GATK4.base_file_name\": \"$sample_name\",\n";
      }elsif($row =~ /PreProcessingForVariantDiscovery_GATK4.final_gvcf_base_name/){
        print $fo "    \"PreProcessingForVariantDiscovery_GATK4.final_gvcf_base_name\": \"$sample_name\",\n";
      }elsif($row =~ /PreProcessingForVariantDiscovery_GATK4.raw_fastq1/){
        print $fo "    \"PreProcessingForVariantDiscovery_GATK4.raw_fastq1\": \"$fastq1\",\n";
      }elsif($row =~ /PreProcessingForVariantDiscovery_GATK4.raw_fastq2/){
        print $fo "    \"PreProcessingForVariantDiscovery_GATK4.raw_fastq2\": \"$fastq2\",\n";
      }elsif($row =~ /PreProcessingForVariantDiscovery_GATK4.bwa_commandline/){
        my $newbwa = replaceSampleNameBWACommand($row, $sample_name);
        print $fo $newbwa;
      }elsif($row =~ /PreProcessingForVariantDiscovery_GATK4.wgs_calling_interval_list/){
        if($calling_interval_file ne ""){
          print $fo "    \"PreProcessingForVariantDiscovery_GATK4.wgs_calling_interval_list\": \"$calling_interval_file\",\n";
        }else{
          print $fo $row;
        }
      }else{
        print $fo $row;
      }
    }
    close($fh);
    close($fo);
    
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

    print $pbs "
gcloud config set project $project

aa=\$(gsutil stat $final_file)
echo \"\"

if [[ \$aa ]]; then
  echo result has been generated $final_file, exit.
  exit 0
fi
  
gcloud \\
  alpha genomics pipelines run \\
  --pipeline-file $pipeline_file \\
  --inputs-from-file WDL=$wdl_file,\\
WORKFLOW_INPUTS=$sample_input_file,\\
WORKFLOW_OPTIONS=$input_option_file \\
  --env-vars WORKSPACE=$gsworkspace,\\
OUTPUTS=${gsresult}/${sample_name} \\
  --logging ${gslog}/${sample_name}.log $option 2>&1 | tee ${pbs_file}.id
    
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  
  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
