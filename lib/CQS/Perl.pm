#!/usr/bin/perl
package CQS::Perl;

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
use File::Spec;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_perl";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $perlFile = get_option_value( $config->{$section}{perlFile}, "perlFile", 1 );
  my $is_absolute = File::Spec->file_name_is_absolute($perlFile);
  if ( !$is_absolute ) {
    $perlFile = dirname(__FILE__) . "/$perlFile";
  }
  if ( !( -e $perlFile ) ) {
    die("perlFile $perlFile defined but not exists!");
  }

  my $output_ext = get_option( $config, $section, "output_ext", 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my %parameterFiles2 = ();
  if ( has_raw_files( $config, $section, "source2" ) ) {
    %parameterFiles2 = %{ get_raw_files( $config, $section, "source2" ) };
  }
  my %parameterFiles3 = ();
  if ( has_raw_files( $config, $section, "source3" ) ) {
    %parameterFiles3 = %{ get_raw_files( $config, $section, "source3" ) };
  }
  my %parameterFiles4 = ();
  if ( has_raw_files( $config, $section, "source4" ) ) {
    %parameterFiles4 = %{ get_raw_files( $config, $section, "source4" ) };
  }
  my %parameterFiles5 = ();
  if ( has_raw_files( $config, $section, "source5" ) ) {
    %parameterFiles5 = %{ get_raw_files( $config, $section, "source5" ) };
  }

  my $shfile = $pbs_dir . "/${task_name}_perl.sh";
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  if ($sh_direct) {
    print $sh "export MYCMD=\"bash\" \n";
  }
  else {
    print $sh "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
  }

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sampleCount = scalar(@sample_files);
    my $samples     = join( ' ', @sample_files );

    my $parameterFile2 = "";
    if ( defined $parameterFiles2{$sample_name} ) {
      my @files = @{ $parameterFiles2{$sample_name} };
      my $file  = $files[0];
      $parameterFile2 = "$file";
    }

    my $parameterFile3 = "";
    if ( defined $parameterFiles3{$sample_name} ) {
      my @files = @{ $parameterFiles3{$sample_name} };
      my $file  = $files[0];
      $parameterFile3 = "$file";
    }

    my $parameterFile4 = "";
    if ( defined $parameterFiles4{$sample_name} ) {
      my @files = @{ $parameterFiles4{$sample_name} };
      my $file  = $files[0];
      $parameterFile4 = "$file";
    }

    my $parameterFile5 = "";
    if ( defined $parameterFiles5{$sample_name} ) {
      my @files = @{ $parameterFiles5{$sample_name} };
      my $file  = $files[0];
      $parameterFile5 = "$file";
    }

    #    my $cur_dir      = create_directory_or_die( $result_dir . "/$sample_name" );

    #    my $pbs_name = "${sample_name}_perl.pbs";
    #    my $pbs_file = "${pbs_dir}/$pbs_name";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log     = $self->get_log_filename( $log_dir, $sample_name );
    my $final_file=$sample_name.$output_ext;

    print $sh "\$MYCMD ./$pbs_name \n";

    #    my $log = "${log_dir}/${sample_name}_perl.log";

    open( my $out, ">$pbs_file" ) or die $!;
    print $out "$pbs_desc
#PBS -o $log
#PBS -j oe

$path_file

echo Perl=`date`
cd $target_dir/result

if [ -s $final_file ]; then
  echo job has already been done. if you want to do again, delete ${result_dir}/${final_file} and submit job again.
  exit 0;
fi

perl $perlFile $sample_name$output_ext $option $samples $parameterFile2 $parameterFile3 $parameterFile4 $parameterFile5

echo finished=`date`
";
    close $out;

    print "$pbs_file created \n";
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all Perl tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $output_ext = get_option( $config, $section, "output_ext", 0 );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {

    #    my @sample_files = @{ $raw_files{$sample_name} };
    #    my @result_files = ();
    #    for my $sampleFile (@sample_files) {
    #      my $name = basename($sampleFile);
    #		push( @result_files, "${result_dir}/${name}${output_ext}" );
    #    }
    my @result_files = ();
    push( @result_files, "${result_dir}/${sample_name}${output_ext}" );

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
