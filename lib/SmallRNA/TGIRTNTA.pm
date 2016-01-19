#!/usr/bin/perl
package SmallRNA::TGIRTNTA;

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
  $self->{_suffix} = "_tt";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $cqstools = get_cqstools( $config, $section, 1 );
  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $extension = get_option( $config, $section, "extension" );

  my %ccaFiles = %{ get_raw_files( $config, $section, "ccaFile" ) };

  my %seqcount_files = ();
  if ( has_raw_files( $config, $section, "seqcount" ) ) {
    %seqcount_files = %{ get_raw_files( $config, $section, "seqcount" ) };
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my $sampleFile = $raw_files{$sample_name}->[0];
    my $ccaFile    = $ccaFiles{$sample_name}->[0];

    my $final_file   = $sample_name . $extension;
    my $summaryFile = $sample_name . $extension . ".summary";

    my $seqcountFile = "";
    if ( defined $seqcount_files{$sample_name} ) {
      my $seqcount = $seqcount_files{$sample_name}->[0];
      $seqcountFile = " -c $seqcount";
    }

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

if [ -s $final_file ]; then
  echo job has already been done. if you want to do again, delete $final_file and submit job again.
  exit 0
fi

echo FastqTrna=`date` 

mono $cqstools tgirt_nta $option -i $sampleFile --ccaFile $ccaFile -o $final_file -s $summaryFile $seqcountFile

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

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $extension = get_option( $config, $section, "extension" );

  my %seqcount_files = ();
  if ( defined $config->{$section}{"seqcount"} || defined $config->{$section}{"seqcount_ref"} ) {
    %seqcount_files = %{ get_raw_files( $config, $section, "seqcount" ) };
  }

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my $final_file   = $result_dir . "/" . $sample_name . $extension;
    my $summaryFile = $result_dir . "/" . $sample_name . $extension . ".summary";

    my @result_files = ();
    push( @result_files, $final_file );
    if ( defined $seqcount_files{$sample_name} ) {
      push( @result_files, $final_file . ".dupcount" );
    }
    push( @result_files, $summaryFile );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
