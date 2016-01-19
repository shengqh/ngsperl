#!/usr/bin/perl
package QC::FastQC;

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
  $self->{_suffix} = "_fq";
  bless $self, $class;
  return $self;
}

sub parsePairedSamples {
  my ($samples)   = @_;
  my @sample_files = @{$samples};
  my @result      = ();
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $summaryfile = $self->get_task_filename( $pbs_dir, $task_name . "_summary" );
  open( my $sh1, ">$summaryfile" ) or die "Cannot create $summaryfile";
  print $sh1 "cd $result_dir
qcimg2pdf.sh -o $task_name
";
  close $sh1;

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  my $result = $self->result( $config, $section );

  for my $sample_name ( sort keys %raw_files ) {
    my @originalFiles = @{ $raw_files{$sample_name} };
    my @sample_files   = parsePairedSamples( \@originalFiles );

    my $sampleCount = scalar(@sample_files);
    my $samples     = join( ' ', @sample_files );
    my $cur_dir      = create_directory_or_die( $result_dir . "/$sample_name" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log     = $self->get_log_filename( $log_dir, $sample_name );
    my $logdesc = $cluster->get_log_description($log);

    my @expectresult = @{ $result->{$sample_name} };
    my $expectname   = $expectresult[0];

    print $sh "\$MYCMD ./$pbs_name \n";

    open( my $out, ">$pbs_file" ) or die $!;
    print $out "$pbs_desc
$logdesc

$path_file

echo fastqc=`date`

if [ -e $expectname ]; then
  echo job has already been done. if you want to do again, delete ${expectname} and submit job again.
  exit 0;
fi

fastqc $option --extract -t $sampleCount -o $cur_dir $samples

echo finished=`date`
";
    close $out;

    print "$pbs_file created \n";
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @originalFiles = @{ $raw_files{$sample_name} };
    my @sample_files   = parsePairedSamples( \@originalFiles );
    my @result_files   = ();
    for my $sampleFile (@sample_files) {
      my $name = basename($sampleFile);
      if ( $name =~ /gz$/ ) {
        $name = change_extension( $name, "" );
      }
      $name = change_extension( $name, "_fastqc" );
      push( @result_files, "${result_dir}/${sample_name}/${name}" );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
