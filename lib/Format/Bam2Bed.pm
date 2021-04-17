#!/usr/bin/perl
package Format::Bam2Bed;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_b2b";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );
  
  my $sort_memory = $thread == 1? $memory:"4G";

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $blacklistfile = get_param_file( $config->{$section}{"blacklist_file"}, "blacklist_file", 0 );
  my $isPairedEnd       = get_option( $config, $section, "is_paired_end" );
  my $minFragmentLength = 0;
  my $maxFragmentLength = 0;
  if ($isPairedEnd) {
    $minFragmentLength = get_option( $config, $section, "minimum_fragment_length" );
    $maxFragmentLength = get_option( $config, $section, "maximum_fragment_length" );
    
    $option = $option . " --min-fragment-length $minFragmentLength --max-fragment-length $maxFragmentLength"
  }
  my $shiftForward = get_option( $config, $section, "shift_forward", 0 );
  my $shiftReverse = get_option( $config, $section, "shift_reverse", 0 );
  my $minMAPQ = get_option( $config, $section, "minimum_mapq", 10 );
  $option = $option . " --shift_forward  $shiftForward --shift_reverse $shiftReverse --min-mapq $minMAPQ";
  
  my $python_bam2bed = dirname(__FILE__) . "/bam2bed.py";
  if ( !-e $python_bam2bed ) {
    die "File not found : " . $python_bam2bed;
  }
  
  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $bam_file     = $sample_files[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";

    my $final_file = $sample_name . ".bed";
    
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    my $bed_file = $sample_name . ".converted.bed";
    my $rmlist = "$bed_file";
    print $pbs "
if [ ! -s $bed_file ]; then
  echo bamtobed=`date` 
  python3 $python_bam2bed $option -i $bam_file -o $bed_file 
fi
";
    if ( defined $blacklistfile ) {
      my $confident_file = $sample_name . ".confident.bed";
      print $pbs "
if [[ -s $bed_file && ! -s $confident_file ]]; then
  echo remove_read_in_blacklist=`date` 
  bedtools intersect -v -a $bed_file -b $blacklistfile > $confident_file
fi
";
      $rmlist   = $rmlist . " " . $confident_file;
      $bed_file = $confident_file;
    }

    print $pbs "
if [ -s $bed_file ]; then
  echo sort_bed=`date` 
  sort -k1,1V -k2,2n -k3,3n $bed_file > $final_file
  rm $rmlist;
fi
";
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();

    my $final_file = $sample_name . ".bed";
    push( @result_files, "${result_dir}/${final_file}" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
