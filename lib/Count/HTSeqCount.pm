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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $gffFile = parse_param_file( $config, $section, "gff_file", 1 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $ispaired       = get_option_value( $config->{$section}{ispairend},      0 );
  my $sorted_by_name = get_option_value( $config->{$section}{sorted_by_name}, 0 );

  my $stranded = get_option_value( $config->{$section}{stranded}, "no" );
  my $strandedoption = "-s " . $stranded;

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @bam_files  = @{ $raw_files{$sample_name} };
    my $bam_file   = $bam_files[0];
    my $final_file = "${sample_name}.count";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";
    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    my $format = ( $bam_file =~ /.sam$/ ) ? "-f sam" : "-f bam";
    my $count_bam_file = $bam_file;
    if ($ispaired) {
      if ( !$sorted_by_name ) {
        $count_bam_file = "${sample_name}_sortedByName.bam";
        $format         = "-f bam";
        print $pbs "samtools sort -n -@ $thread -o $count_bam_file $bam_file";
      }
      $format = $format . " -r name";
    }

    print $pbs "htseq-count $option $format -q -m intersection-nonempty $strandedoption -i gene_id $count_bam_file $gffFile > $final_file \n";
    
    if($count_bam_file ne $bam_file){
      print $pbs "if [ -s $final_file ]; then
  rm $count_bam_file
fi
";
    }

    $self->close_pbs( $pbs, $pbs_file );

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    my $countFile    = "${result_dir}/${sample_name}.count";
    push( @result_files, $countFile );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
