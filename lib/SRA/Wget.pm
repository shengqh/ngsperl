#!/usr/bin/perl
package SRA::Wget;

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
use SRA::SRAUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_fd";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $ispaired = get_is_paired_end_option($config, $section, 0);
  
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

    my $current_dir = $result_dir;

    my $final_file = $ispaired ? $sample_name . "_1.fastq.gz" : $sample_name . ".fastq.gz";
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $current_dir, $final_file );

    if ( $sample_file =~ /^GSM/ ) {
      $sample_file = GsmToSrr( $sample_file );
    }

    $sample_file =~ /^SRR/ or die "$sample_file is not SRR id";

    my $urls = SrrToUrl( $sample_file );
    
    if ( scalar(@$urls) == 1 ) {
      my $url = $urls->[0];
      print $pbs "
rm -f $sample_name.failed $sample_name.succeed

wget -c ftp://$url -O $final_file
        
status=\$?
if [[ \$status -ne 0 ]]; then
  touch $sample_name.failed
  rm -f $final_file
else
  touch $sample_name.succeed
fi
";
    }else{
      my $f1 = $urls->[0];
      my $f2 = $urls->[1];
      print $pbs "
rm -f $sample_name.failed $sample_name.succeed

wget -c ftp://$f1 -O ${sample_name}_1.fastq.gz
       
status=\$?
if [[ \$status -ne 0 ]]; then
  touch $sample_name.failed
  rm -f ${sample_name}_1.fastq.gz
else
  wget -c ftp://$f2 -O ${sample_name}_2.fastq.gz

  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sample_name.failed
    rm -f ${sample_name}_1.fastq.gz ${sample_name}_2.fastq.gz
  else
    touch $sample_name.succeed
  fi
fi

";
    }

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

  my $ispaired = get_is_paired_end_option($config, $section, 0);

  my $raw_files = getSraFiles( $config, $section );

  my $result = {};
  for my $sample_name ( keys %$raw_files ) {

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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $raw_files = getSraFiles( $config, $section );

  my $result = {};

  for my $sample_name ( sort keys %$raw_files ) {
    $result->{$sample_name} = $self->get_pbs_filename( $pbs_dir, $sample_name );
  }

  return $result;
}

1;
