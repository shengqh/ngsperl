#!/usr/bin/perl
package QC::FastQC;

use strict;
use warnings;
use File::Basename;
use List::Util qw(min);
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
  $self->{_suffix} = "_fq";
  $self->{_use_tmp_folder} = 1;
  bless $self, $class;
  $self->init_docker_prefix(__PACKAGE__);
  return $self;
}

sub parsePairedSamples {
  my ($samples)    = @_;
  my @sample_files = @{$samples};
  my @result       = ();
  for my $sample (@sample_files) {
    if ( $sample =~ /,/ ) {
      my @files = split( ',', $sample );
      for my $file (@files) {
        push( @result, $file );
      }
    }
    else {
      push( @result, $sample );
    }
  }

  return @result;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory, $init_command ) = $self->init_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $fastqc=get_option($config, $section, "fastqc", "fastqc");

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  my $result = $self->result( $config, $section );

  for my $sample_name ( sort keys %raw_files ) {
    my @originalFiles = @{ $raw_files{$sample_name} };
    my @sample_files  = parsePairedSamples( \@originalFiles );
    my $first_sample_file = $sample_files[0];

    my $sampleCount = scalar(@sample_files);
    my $cur_dir     = create_directory_or_die( $result_dir . "/$sample_name" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $log_desc  = $cluster->get_log_description($log);

    my @expectresult = @{ $result->{$sample_name} };
    my $expectname   = $expectresult[0];

    my $curThreadCount = min($thread, $sampleCount);

    if (scalar(@expectresult) == 1){
      print $sh "
if [[ ! -s $expectname ]]; then
  \$MYCMD ./$pbs_name 
fi

";
    }else{
      my $expectname2   = $expectresult[1];
      print $sh "
if [[ ! -s $expectname || ! -s $expectname2 ]]; then
  \$MYCMD ./$pbs_name 
fi

";
    }

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, \@expectresult, "", 0, $first_sample_file );

    my $localized_files = [];
    @sample_files = @{$self->localize_files_in_tmp_folder($pbs, \@sample_files, $localized_files)};

    my $samples     = '"' . join( '" "', @sample_files ) . '"';
    print $pbs "
rm -f $sample_name.fastqc.failed

$fastqc $option --extract -t $curThreadCount -o `pwd` $samples 2> >(tee ${sample_name}.fastqc.stderr.log >\&2)

status=\$?
if [[ \$status -ne 0 ]]; then
  touch $sample_name.fastqc.failed
else
  touch $sample_name.fastqc.succeed
fi

$fastqc --version | cut -d ' ' -f2 | awk '{print \"FastQC,\"\$1}' > `pwd`/fastqc.version

rm -rf .cache .java

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

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @originalFiles = @{ $raw_files{$sample_name} };
    my @sample_files  = parsePairedSamples( \@originalFiles );
    my @result_files  = ();
    for my $sampleFile (@sample_files) {
      my $name = basename($sampleFile);
      if($name =~ /\.gz$/ ) {
        $name = change_extension( $name, "" );
      }
      if($name =~ /\.fq$/){
        $name = change_extension( $name, "" );
      }
      if($name =~ /\.fastq$/){
        $name = change_extension( $name, "" );
      }
      $name = $name . "_fastqc";
      push( @result_files, "${result_dir}/${sample_name}/${name}/fastqc_data.txt" );
    }
    push( @result_files, "${result_dir}/${sample_name}/fastqc.version" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
