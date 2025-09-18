#!/usr/bin/perl
package SRA::FastqDumpPaired;

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
  $self->{_docker_prefix} = "sratools_";
  $self->{_docker_shell} = "sh";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $ispaired = get_is_paired_end_option($config, $section, 0);
  if(!$ispaired){
    die "SRA::FastqDumpPaired should be used for paired end only.";
  }
  my $is_restricted_data = get_option($config, $section, "is_restricted_data" , 0);
  my $prefetch_option = get_option($config, $section, "prefetch_option", "");

  my $single_cell_data_type = get_option($config, $section, "single_cell_data_type", 0);

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

    my $current_dir = create_directory_or_die($result_dir . "/" . $sample_name);

    my $final_file_1;
    my $final_file_2;
    if($single_cell_data_type == 2){ #R1 and R2
      $final_file_1 = $sample_name . "_S1_L001_R1_001.fastq.gz";
      $final_file_2 = $sample_name . "_S1_L001_R2_001.fastq.gz";
    }else{
      $final_file_1 = $sample_name . "_1.fastq.gz";
      $final_file_2 = $sample_name . "_2.fastq.gz";
    }

    print $sh "if [[ ! -s $current_dir/$final_file_1 ]]; then
  \$MYCMD ./$pbs_name 
fi
";

    # if( -s $pbs_file ){
    #   next;
    # }

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $current_dir, $final_file_1, "", 0, undef, 'sh' );

    if ( $sample_file =~ /GSM/ ) {
      my $srr_file = GsmToSrr( $sample_file );
      if ($srr_file eq "") {
        die "Cannot get SRR for $sample_file";
      }
      $sample_file = $srr_file;
    }
    if ( $sample_file =~ /\.sra/ ) {
      my $dump_file_1 = "${sample_name}_1.fastq.gz";
      my $dump_file_2 = "${sample_name}_2.fastq.gz";
      if ($is_restricted_data){
        print $pbs "
ln -s $sample_file ${sample_name}.sra 
rm -f $sample_name.failed

fastq-dump $option ${sample_name}.sra

status=\$?
if [[ \$status -ne 0 ]]; then
  touch $sample_name.failed
  rm -f $dump_file_1 $dump_file_2
else
  touch $sample_name.succeed
  if [[ $dump_file_1 != $final_file_1 ]]; then #for single cell data
    mv $dump_file_1 $final_file_1
    mv $dump_file_2 $final_file_2
  fi
fi
rm -f ${sample_name}.sra
";
      } else {
        print $pbs "
if [ -z \${SLURM_JOBID+x} ]; then 
  echo \"in bash mode\"; 
  ln -s $sample_file ${sample_name}.sra 
  fastq-dump $option ${sample_name}.sra

  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sample_name.failed
    rm -f $dump_file_1 $dump_file_2
  else
    touch $sample_name.succeed
    if [[ $dump_file_1 != $final_file_1 ]]; then #for single cell data
      mv $dump_file_1 $final_file_1
      mv $dump_file_2 $final_file_2
    fi
  fi
  rm -f ${sample_name}.sra
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
  fastq-dump $option ${sample_name}.sra

  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sample_name.failed
    rm -f $dump_file_1 $dump_file_2
  else
    touch $sample_name.succeed
    if [[ '$dump_file_1' != '$final_file_1' ]]; then #for single cell data
      mv $dump_file_1 $final_file_1
      mv $dump_file_2 $final_file_2
    fi
  fi
  rm -f ${sample_name}.sra

  if [[ -s $final_file_1 ]]; then
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
        my $dump_file_1 = "${sample_file}_1.fastq.gz";
        my $dump_file_2 = "${sample_file}_2.fastq.gz";
        print $pbs "
rm -f $sample_name.failed

if [[ ! -s ${sample_file}.sra ]]; then
  rm -f $sample_name.prefetch.failed $sample_name.prefetch.succeed

  echo prefetch $sample_file $prefetch_option --max-size u --check-rs no --progress-o ${sample_file}.tmp.sra
  prefetch $sample_file $prefetch_option --max-size u --check-rs no --progress -o ${sample_file}.tmp.sra
  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sample_name.prefetch.failed
    rm -f ${sample_file}.tmp.sra
  else
    touch $sample_name.prefetch.succeed
    mv ${sample_file}.tmp.sra ${sample_file}.sra
  fi
fi

if [[ -s ${sample_file}.sra ]]; then
  echo fastq-dump $option ${sample_file}.sra 
  fastq-dump $option ${sample_file}.sra 
  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sample_name.fastq-dump.failed
    rm -f $dump_file_1 $dump_file_2 ${sample_file}.fastq.gz
  else
    touch $sample_name.fastq-dump.succeed
    if [[ '$dump_file_1' != '$final_file_1' ]]; then
      mv $dump_file_1 $final_file_1
      mv $dump_file_2 $final_file_2
    fi
    rm -rf ${sample_file}.sra
  fi
fi
";
      }
      else {
        print $pbs "rm -rf ${sample_name}_1.fastq.gz ${sample_name}_2.fastq.gz \n";
        for my $sf (@sample_files) {
          print $pbs "

if [[ ! -s ${sf}.sra ]]; then
  echo prefetch $sf $prefetch_option -o ${sf}.tmp.sra
  prefetch $sf $prefetch_option --check-rs no -o ${sf}.tmp.sra
  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sf.prefetch.failed
    rm -f ${sf}.tmp.sra
    exit 1
  else
    touch $sf.prefetch.succeed
    rm -f $sf.prefetch.failed
    mv ${sf}.tmp.sra ${sf}.sra
  fi
fi

if [[ -s ${sf}.sra ]]; then
  echo fastq-dump $option ${sf}.sra 
  fastq-dump $option ${sf}.sra 
  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sample_name.fasterq.failed
    rm -f ${sf}_1.fastq.gz ${sf}_2.fastq.gz ${sf}.fastq.gz
    exit 1
  fi
fi

cat ${sf}_1.fastq.gz >> ${sample_name}_1.fastq.gz
cat ${sf}_2.fastq.gz >> ${sample_name}_2.fastq.gz
rm -rf ${sf}_1.fastq.gz ${sf}_2.fastq.gz ${sf}.sra

";
        }
        print $pbs "
if [[ '${sample_name}_1.fastq.gz' != '$final_file_1' ]]; then #for single cell data
  mv ${sample_name}_1.fastq.gz $final_file_1
  mv ${sample_name}_2.fastq.gz $final_file_2
fi
";
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

  my $raw_files = getSraFiles( $config, $section );
  my $is_single_cell_data = get_option($config, $section, "is_single_cell_data", 0);

  my $result = {};
  for my $sample_name ( keys %$raw_files ) {

    my @result_files = ();
    push( @result_files, $result_dir . "/$sample_name/" . $sample_name . "_1.fastq.gz" );
    push( @result_files, $result_dir . "/$sample_name/" . $sample_name . "_2.fastq.gz" );

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
