#!/usr/bin/perl
package GATK4::GenotypeGVCFsScatter;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use GATK4::GATK4Task;
use Data::Dumper;
use Tie::IxHash;

our @ISA = qw(GATK4::GATK4Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_gg";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $ref_fasta = get_option_file( $config, $section, "fasta_file", 1 );
  my $dbsnp_vcf = get_option_file( $config, $section, "dbsnp_vcf", 1 );

  my $import_folders = get_raw_files( $config, $section );
  my $interval_files = get_interval_file_map($config, $section);
  #print(Dumper($interval_files));

  my $task_interval_files = {};
  tie %$task_interval_files, 'Tie::IxHash';

  for my $scatter_name (sort keys %$interval_files) {
    my $prefix = get_key_name($task_name, $scatter_name);
    $task_interval_files->{$prefix} = $interval_files->{$scatter_name};
    if (not defined $import_folders->{$prefix}) {
      die "$prefix not in import_folders " . Dumper($import_folders);
    }
  }
  #print(Dumper($task_interval_files));

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $scatter_name ( sort keys %$import_folders ) {
    my $WORKSPACE = $import_folders->{$scatter_name}[0];
    my $interval = $task_interval_files->{$scatter_name};

    if (not defined $interval){
      die "$scatter_name not in task_interval_files " . Dumper($task_interval_files);
    }

    my $output_vcf_filename    = $scatter_name . ".g.vcf.gz";
    my $output_vcf_filename_index = $output_vcf_filename . ".tbi";

    my $tmp_file    = $scatter_name . ".tmp.g.vcf.gz";
    my $tmp_index = $tmp_file . ".tbi";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $scatter_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $scatter_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $output_vcf_filename_index );
    print $pbs "
gatk --java-options \"-Xmx8g -Xms8g\" \\
  GenotypeGVCFs \\
  -R ${ref_fasta} \\
  -O ${tmp_file} \\
  -D ${dbsnp_vcf} \\
  -G StandardAnnotation -G AS_StandardAnnotation \\
  --only-output-calls-starting-in-intervals \\
  --use-new-qual-calculator \\
  -V gendb://$WORKSPACE \\
  -L ${interval} \\
  --merge-input-intervals

if [[ -s $tmp_index ]]; then
  mv $tmp_file $output_vcf_filename
  mv $tmp_index $output_vcf_filename_index
fi

";
    
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print " !!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks. \n ";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;
  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );
  my $result = {};

  my $workspaces = get_raw_files( $config, $section );
  for my $scatter_name ( sort keys %$workspaces ) {
    my $snvOut    = $scatter_name . ".g.vcf.gz";
    my @result_files = ();
    push( @result_files, "${result_dir}/${snvOut}" );
    $result->{$scatter_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
