#!/usr/bin/perl
package Alignment::BWA_WGS;

use strict;
use warnings;
use File::Basename;
use POSIX;
use List::Util qw[min max];
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
  $self->{_suffix} = "_bwa";
  $self->init_docker_prefix(__PACKAGE__);
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my ($sort_memory, $isMB) = getMemoryPerThread($memory, $thread);
  my $max_useful_memory = max(5, ($sort_memory * $thread - 10));
  if ($isMB) {
    $sort_memory = $max_useful_memory . "M";
  }else{
    $sort_memory = $max_useful_memory . "G";
  }

  my $selfname = $self->{_name};

  my $rg_name_regex        = get_option( $config, $section, "rg_name_regex",      "" );
  my $rg_id_regex        = get_option( $config, $section, "rg_id_regex",      "" );

  if ( !( $option =~ /\s-t\s/ ) ) {
    if ( $thread > 1 ) {
      $option = $option . " -t " . $thread;
    }
  }

  my $bwa_index = $config->{$section}{bwa_index};
  if ( !defined $bwa_index ) {
    $bwa_index = $config->{$section}{fasta_file} or die "define ${section}::bwa_index first";
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sample_file_0 = $sample_files[0];
    my $sample_files_str = ( scalar(@sample_files) == 2 ) ? "\"" . $sample_file_0 . "\" \"" . $sample_files[1] . "\"" : "\"" . $sample_file_0 . "\"";

    my $rg_sample_name = $sample_name;
    if ($rg_name_regex ne ""){
      if ($sample_name =~ /$rg_name_regex/ ) {
        $rg_sample_name = $1;
      }
    }
    my $rg_sample_id = "1";
    if ($rg_id_regex ne ""){
      if ($sample_name =~ /$rg_id_regex/ ) {
        $rg_sample_id = $1;
      }
    }

    my $rg = "-R \"\@RG\\tID:${rg_sample_id}\\tSM:${rg_sample_name}\\tPL:ILLUMINA\"";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $final_file =  $sample_name . ".sortedByCoord.bam";
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file, $init_command, 0, $sample_file_0 );

    my $localized_files = [];
    @sample_files = @{$self->localize_files_in_tmp_folder($pbs, \@sample_files, $localized_files)};

    print $pbs "
echo bwa_mem=`date`

rm ${final_file}.failed

bwa mem $option $rg $bwa_index $sample_files_str | samtools view -bhu - | sambamba sort -m $memory -l 0 -u -t $thread --tmpdir tmp -o ${final_file} /dev/stdin

status=\$?
if [[ \$status -eq 0 ]]; then
  touch ${final_file}.succeed
  bwa 2>\&1 | grep Version | cut -d ' ' -f2 | cut -d '-' -f1 | awk '{print \"bwa,v\"\$1}' > ${sample_name}.bwa.version
  sambamba index ${final_file}
else
  rm $final_file
  touch ${final_file}.failed
fi

";

    $self->clean_temp_files($pbs, $localized_files);
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $sortByCoordinate = get_option( $config, $section, "sort_by_coordinate", 1 );
  my $mark_duplicates  = hasMarkDuplicate( $config->{$section} );
  my %raw_files        = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    push( @result_files, "${result_dir}/${sample_name}.sortedByCoord.bam" );
    push( @result_files, "${result_dir}/${sample_name}.sortedByCoord.bam.bai" );
    push( @result_files, "${result_dir}/${sample_name}.bwa.version" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
