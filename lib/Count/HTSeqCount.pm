#!/usr/bin/perl
package Count::HTSeqCount;

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
  $self->{_suffix} = "_ht";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $gffFile = get_param_file( $config->{$section}{gff_file}, "gff_file", 1 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $ispaired = get_option_value( $config->{$section}{ispairend}, "ispairend", 0 );
  my $ispairoption = $ispaired ? "-f 1" : "";

  my $stranded = get_option_value( $config->{$section}{stranded}, "stranded", "no" );
  my $strandedoption = "-s " . $stranded;

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @bam_files  = @{ $raw_files{$sample_name} };
    my $bam_file   = $bam_files[0];
    my $countFile = "${sample_name}.count";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log     = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    open( my $out, ">$pbs_file" ) or die $!;
    print $out "$pbs_desc
$log_desc

$path_file

cd $result_dir

if [ -s $countFile ]; then
  echo job has already been done. if you want to do again, delete $countFile and submit job again.
  exit 0
fi

echo HTSeqCount=`date`

samtools view $ispairoption $bam_file | htseq-count $option -q -m intersection-nonempty $strandedoption -i gene_id - $gffFile > $countFile

echo finished=`date`

exit 0 
";

    close $out;

    print "$pbs_file created \n";
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all mirna_count tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    my $countFile   = "${result_dir}/${sample_name}.count";
    push( @result_files, $countFile );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
