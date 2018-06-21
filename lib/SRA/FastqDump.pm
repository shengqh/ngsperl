#!/usr/bin/perl
package SRA::FastqDump;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;
use LWP::Simple;
use LWP::UserAgent;
use URI::Escape;
use List::Util qw(first);
use String::Util qw(trim);
use Utils::CollectionUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_fd";
  bless $self, $class;
  return $self;
}

sub getGsmSrrMap {
  my $listFile = shift;
  open my $list_handle, "<$listFile";
  my $first_line = <$list_handle>;
  close $list_handle;
  my @columns = split( '\t', $first_line );
  my $srrIndex = first { $columns[$_] eq 'Run_s' } 0 .. $#columns;
  my $gsmIndex = first { $columns[$_] eq 'Sample_Name_s' } 0 .. $#columns;
  my $taMap = readDictionaryByIndex( $listFile, $gsmIndex, $srrIndex, 1 );
  return $taMap;
}

sub getSraFiles {
  my ( $config, $section ) = @_;
  my $result;
  if ( defined $config->{$section}{"list_file"} ) {
    my $taMap = getGsmSrrMap( $config->{$section}{"list_file"} );
    for my $gsm ( keys %$taMap ) {
      $result->{$gsm} = [ $taMap->{$gsm} ];
    }
  }
  else {
    if ( defined $config->{$section}{"source"} ) {
      my $res = $config->{$section}{"source"};
      if ( ref($res) eq 'ARRAY' ) {
        for my $gsm (@$res) {
          $result->{$gsm} = [$gsm];
        }
      }
      else {
        $result = $res;
      }
    }
    else {
      my $fileSection = $config->{$section}{"source_ref"}[0];
      my $files       = $config->{$fileSection};
      if ( ref($files) eq 'ARRAY' ) {
        for my $gsm (@$files) {
          $result->{$gsm} = [$gsm];
        }
      }
      else {
        $result = get_raw_files( $config, $section );
      }
    }
  }
  return $result;
}

sub GsmToSrr {
  my ($gsm, $sraTable) = @_;

  my $cmd = "grep $gsm $sraTable | grep -e \"^SRR\" | cut -f1";
  my $res = ` $cmd `;
  $res = trim($res);
  return $res;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $ispaired = $config->{$section}{ispaired};
  if ( !defined $ispaired ) {
    $ispaired = 0;
  }

  my $raw_files = getSraFiles( $config, $section );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %$raw_files ) {
    my $sample_file = $raw_files->{$sample_name}->[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";

    my $final_file = $ispaired ? $sample_name . "_1.fastq.gz" : $sample_name . ".fastq.gz";
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    if ( $sample_file =~ /GSM/ ) {
      my $sraTable  = get_option_file($config, $section, "sra_table");
      $sample_file = GsmToSrr($sample_file, $sraTable);
    }
    if ( $sample_file =~ /\.sra/ ) {
      print $pbs "ln -s $sample_file ${sample_name}.sra \n";
    }
    elsif ( $sample_file =~ /SRR/ ) {
      my $six = substr( $sample_file, 0, 6 );
      print $pbs "wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${six}/${sample_file}/${sample_file}.sra -O ${sample_name}.sra\n";
    }
    else {
      die "I don't know what it is " . $sample_file;
    }
    print $pbs "fastq-dump --split-3 --gzip --origfmt --helicos ${sample_name}.sra
rm ${sample_name}.sra
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $ispaired = $config->{$section}{ispaired};
  if ( !defined $ispaired ) {
    $ispaired = 0;
  }

  my $raw_files = getSraFiles( $config, $section );

  my $result = {};
  for my $sample_name ( keys %$raw_files ) {

    my $final_file = $ispaired ? $sample_name . "_1.fastq.gz" : $sample_name . ".fastq.gz";

    my @result_files = ();
    if ($ispaired) {
      push( @result_files, $result_dir . "/" . $sample_name . "_1.fastq.gz" );
      push( @result_files, $result_dir . "/" . $sample_name . "_2.fastq.gz" );
    }
    else {
      push( @result_files, $result_dir . "/" . $sample_name . ".fastq.gz" );
    }

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $raw_files = getSraFiles( $config, $section );

  my $result = {};

  for my $sample_name ( sort keys %$raw_files ) {
    $result->{$sample_name} = $self->get_pbs_filename( $pbs_dir, $sample_name );
  }

  return $result;
}

1;
