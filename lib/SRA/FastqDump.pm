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
  my $is_restricted_data = get_option($config, $section, "is_restricted_data" , 0);

  if($ispaired){
    $option = $option . " --split-3 ";
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

    my $current_dir = create_directory_or_die($result_dir . "/" . $sample_name);

    my $final_file = $ispaired ? $sample_name . "_1.fastq.gz" : $sample_name . ".fastq.gz";
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $current_dir, $final_file );

    if ( $sample_file =~ /GSM/ ) {
      $sample_file = GsmToSrr( $sample_file );
    }
    if ( $sample_file =~ /\.sra/ ) {
      if ($is_restricted_data){
        print $pbs "
ln -s $sample_file ${sample_name}.sra 
rm -f $sample_name.failed

fastq-dump $option --gzip --origfmt --helicos ${sample_name}.sra

status=\$?
if [[ \$status -ne 0 ]]; then
  touch $sample_name.failed
  rm -f $final_file
else
  touch $sample_name.succeed
fi
rm -f ${sample_name}.sra
";
      } else {
        print $pbs "
if [ -z \${SLURM_JOBID+x} ]; then 
  echo \"in bash mode\"; 
  ln -s $sample_file ${sample_name}.sra 
  fastq-dump --split-e --gzip --origfmt --helicos ${sample_name}.sra
  rm ${sample_name}.sra
else 
  echo \"in cluster mode\"; 
  localdir=/tmp/myjob_\${SLURM_JOBID}
  tmp_cleaner()
  {
  rm -rf \${localdir}
  exit -1
  }
  trap 'tmp_cleaner' TERM

  echo creating local directory \${localdir}
  mkdir \${localdir} # create unique directory on compute node
  cd \${localdir}

  echo copying $sample_file to \${localdir}
  cp $sample_file ${sample_name}.sra

  echo performing fastq-dump on ${sample_name}.sra
  fastq-dump --split-3 --gzip --origfmt --helicos ${sample_name}.sra

  rm ${sample_name}.sra

  if [[ -s $final_file ]]; then
    echo copying result back to $result_dir
    mv -f ${sample_name}* $result_dir
  fi
fi

";
      }
    }
    elsif ( $sample_file =~ /[SD]RR/ ) {
      my $sratoolkit_setting_file = get_option_file($config, $section, "sratoolkit_setting_file");

      print $pbs "
if [[ ! -s \${HOME}/.ncbi ]]; then
  echo mkdir \${HOME}/.ncbi
  mkdir \${HOME}/.ncbi
fi

if [[ ! -s \${HOME}/.ncbi/user-settings.mkfg ]]; then
  echo cp user-settings.mkfg
  cp $sratoolkit_setting_file \${HOME}/.ncbi
fi

echo dump $sample_name

";

      my @sample_files = split( '\s+', $sample_file );
      if ( scalar(@sample_files) == 1 ) {
        print $pbs "
rm -f $sample_name.failed

fastq-dump --split-e --gzip --origfmt --helicos $sample_file 
        
status=\$?
if [[ \$status -ne 0 ]]; then
  touch $sample_name.failed
  rm -f ${sample_file}_1.fastq.gz ${sample_file}_2.fastq.gz ${sample_file}.fastq.gz
else
  touch $sample_name.succeed
";
        if($ispaired){
          print $pbs "  mv ${sample_file}_1.fastq.gz ${sample_name}_1.fastq.gz \n";
          print $pbs "  mv ${sample_file}_2.fastq.gz ${sample_name}_2.fastq.gz \n";
        }else{
          print $pbs "  mv ${sample_file}.fastq.gz ${sample_name}.fastq.gz \n";
        }
        print $pbs "fi
";
      }
      else {
        print $pbs "if [[ -s ${sample_name}.fastq ]]; then\n  rm ${sample_name}.fastq \nfi \n";
        for my $sf (@sample_files) {
          print $pbs "fastq-dump --split-e --origfmt --helicos $sf -Z >> ${sample_name}.fastq \n";
        }
        print $pbs "gzip ${sample_name}.fastq \n";
      }
    }
    else {
      die "I don't know what it is " . $sample_file;
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
      push( @result_files, $result_dir . "/$sample_name/" . $sample_name . "_1.fastq.gz" );
      push( @result_files, $result_dir . "/$sample_name/" . $sample_name . "_2.fastq.gz" );
    }
    else {
      push( @result_files, $result_dir . "/$sample_name/" . $sample_name . ".fastq.gz" );
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
