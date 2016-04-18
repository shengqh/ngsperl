#!/usr/bin/perl
package Annotation::Annovar;

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
  $self->{_suffix} = "_ann";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $buildver = $config->{$section}{buildver} or die "buildver is not defined in $section";
  $option = "-buildver $buildver $option";

  my $annovarDB = $config->{$section}{annovar_db} or die "annovar_db is not defined in $section";
  my $isvcf = $config->{$section}{isvcf};
  if ( !defined $isvcf ) {
    $isvcf = 0;
  }

  my $cqstools = get_cqstools( $config, $section, 0 );
  my $affyFile = get_param_file( $config->{$section}{affy_file}, "affy_file", 0 );

  my $raw_files = get_raw_files( $config, $section );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  my $listfile = $self->get_file( $result_dir, $task_name, ".list", 0 );
  open( my $lt, ">$listfile" ) or die "Cannot create $listfile";

  for my $sample_name ( sort keys %{$raw_files} ) {
    my @sample_files = @{ $raw_files->{$sample_name} };

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log_file = $self->get_log_filename( $log_dir, $sample_name );

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    my $log_desc = $cluster->get_log_description($log_file);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir );

    for my $sampleFile (@sample_files) {
      my ( $filename, $dir ) = fileparse($sampleFile);

      if ( $dir eq $cur_dir ) {
        $sampleFile = $filename;
      }

      my $annovar = change_extension( $filename, ".annovar" );
      my $result  = "${annovar}.${buildver}_multianno.txt";
      my $final   = $annovar . ".final.tsv";
      my $excel   = $final . ".xls";

      my $vcf;
      my $passinput;
      if ($isvcf) {
        $passinput = change_extension( $filename, ".avinput" );
        $vcf = "convert2annovar.pl -format vcf4old ${sampleFile} | cut -f1-7 > $passinput ";
      }
      else {
        $passinput = $sampleFile;
        $vcf       = "";
      }

      print $pbs "
if [[ ! -s $result && ! -s $final ]]; then 
  $vcf
  table_annovar.pl $passinput $annovarDB $option --outfile $annovar --remove
fi

if [[ -s $result && ! -s $final ]]; then
  grep \"^##\" ${sampleFile} > ${final}.header
  grep -v \"^##\" ${sampleFile} | cut -f7- > ${sampleFile}.clean
  grep -v \"^##\" ${result} > ${result}.clean
  paste ${result}.clean ${sampleFile}.clean > ${final}.data
  cat ${final}.header ${final}.data > $final
  rm ${sampleFile}.clean ${result}.clean ${final}.header ${final}.data
fi
";

      if ( defined $cqstools ) {
        my $affyoption = defined($affyFile) ? "-a $affyFile" : "";
        print $pbs "
if [ -s $final ]; then
  rm $passinput $result
fi

if [[ -s $final && ! -s $excel ]]; then
  mono-sgen $cqstools annovar_refine -i $final $affyoption -o $excel
fi
";
      }

      print $lt "${cur_dir}/${result}\n";
    }
    $self->close_pbs( $pbs, $pbs_file );

    print $sh "\$MYCMD ./$pbs_name \n";

  }
  close $lt;

  print $sh "exit 0\n";
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }
  print "!!!shell file $shfile created, you can run this shell file to submit Annovar tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $buildver = $config->{$section}{buildver} or die "buildver is not defined in $section";
  my $raw_files = get_raw_files( $config, $section );
  my $cqstools = get_cqstools( $config, $section, 0 );

  my $result = {};
  for my $sample_name ( sort keys %{$raw_files} ) {
    my @sample_files = @{ $raw_files->{$sample_name} };
    my $cur_dir      = $result_dir . "/$sample_name";
    my @result_files = ();
    for my $sampleFile (@sample_files) {
      my $annovar = change_extension( $sampleFile, ".annovar" );
      my $final   = $annovar . ".final.tsv";
      if ( defined $cqstools ) {
        my $excel = $final . ".xls";
        push( @result_files, $cur_dir . "/$excel" );
      }
      push( @result_files, $cur_dir . "/$final" );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
