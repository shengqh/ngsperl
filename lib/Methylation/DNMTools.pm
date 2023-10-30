#!/usr/bin/perl
package Methylation::DNMTools;

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
  $self->{_suffix} = "_dnmtools";
  $self->{_docker_prefix} = "dnmtools_";
  $self->{_docker_shell} = "sh";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $selfname = $self->{_name};

  my $dnmtools_command = get_option($config, $section, "dnmtools_command", "dnmtools");

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $chr_fasta = get_option($config, $section, "chr_fasta");
  # my $chrDir=$config->{$section}{chr_dir};
  # if ( !defined $chrDir ) {
  #   die "define ${section}::chr_dir first";
  # }
  my $chrSizeFile=$config->{$section}{chr_size_file};
  if ( !defined $chrSizeFile ) {
    die "define ${section}::chr_size_file first";
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $sampleFile=$sample_files[0];
    my $sampleFileBase=basename($sampleFile);
 
    my $result_file          = "${sample_name}.cpg.meth";
    my $tag               = get_bam_tag($sample_name);

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $result_file );

    print $pbs "
echo DNMTools=`date`

if [ ! -s ${sample_name}.uniq.bam ]; then
  if [ ! -s ${sample_name}.formatted.bam ]; then
    echo dnmtools format=`date`
    $dnmtools_command format -t $thread -B -f abismal -stdout $sampleFile | samtools sort -@ thread -o ${sample_name}.formatted.sorted.bam
  fi

  echo dnmtools uniq=`date`
  $dnmtools_command uniq -t $thread -B -S ${sample_name}.uniq.bam.dupstats ${sample_name}.formatted.sorted.bam ${sample_name}.uniq.bam

  if [[ -s ${sample_name}.uniq.bam ]]; then
    rm ${sample_name}.formatted.sorted.bam
  fi
fi

if [ ! -s ${sample_name}.bsrate ]; then
  echo dnmtools bsrate=`date`
  $dnmtools_command bsrate -t $thread -c $chr_fasta -o ${sample_name}.bsrate ${sample_name}.uniq.bam 
fi

if [ ! -s ${sample_name}.cpg.meth ]; then
  if [ ! -s ${sample_name}.all.meth ]; then 
    #output cpg only
    echo dnmtools counts=`date`
    $dnmtools_command counts -t $thread -c $chr_fasta -o ${sample_name}.all.meth ${sample_name}.uniq.bam
  fi

  if [ ! -s ${sample_name}.levels ]; then
    echo dnmtools levels=`date`
    $dnmtools_command levels -o ${sample_name}.levels ${sample_name}.all.meth
  fi

  echo dnmtools sym=`date`
  $dnmtools_command sym -o ${sample_name}.cpg.meth ${sample_name}.all.meth

  if [[ -s ${sample_name}.cpg.meth ]]; then
    rm ${sample_name}.all.meth
  fi
fi

if [ ! -s ${sample_name}.hmr ]; then
  echo dnmtools hmr=`date`
  $dnmtools_command hmr -o ${sample_name}.hmr -p ${sample_name}.hmrparams ${sample_name}.cpg.meth
fi

# No pmr in dnmtools 1.4.1
# if [ ! -s ${sample_name}.pmr ]; then
#   echo dnmtools pmr=`date`
#   $dnmtools_command hmr -partial -o ${sample_name}.pmr -p ${sample_name}.pmrparams ${sample_name}.cpg.meth
# fi

if [ ! -s ${sample_name}.pmd ]; then
  echo dnmtools pmd=`date`
  $dnmtools_command pmd -o ${sample_name}.pmd -p ${sample_name}.pmdparams ${sample_name}.cpg.meth
fi

if [ ! -s ${sample_name}.epiread ]; then
  echo dnmtools states=`date`
  $dnmtools_command states -t $thread -c $chr_fasta -o ${sample_name}.epiread ${sample_name}.uniq.bam
fi

if [ ! -s ${sample_name}.allelic ]; then
  echo dnmtools allelic=`date`
  $dnmtools_command allelic -c $chr_fasta -o ${sample_name}.allelic ${sample_name}.epiread
fi

if [ ! -s ${sample_name}.amr ]; then
  echo dnmtools amrfinder=`date`
  $dnmtools_command amrfinder -c $chr_fasta -o ${sample_name}.amr ${sample_name}.epiread
fi

if [ ! -s ${sample_name}.read.bw ]; then
  echo DNMTools To Tracks=`date`
  awk '{OFS=\"\\t\"; print \$1,\$2,\$2+1,\$6}' < ${sample_name}.cpg.meth > ${sample_name}.read.bw.tmp 
  wigToBigWig ${sample_name}.read.bw.tmp $chrSizeFile ${sample_name}.read.bw
  rm ${sample_name}.read.bw.tmp
fi

if [ ! -s ${sample_name}.cpg.meth.bw ]; then
  awk '{OFS=\"\\t\"; print \$1,\$2,\$2+1,\$5}' < ${sample_name}.cpg.meth > ${sample_name}.cpg.meth.bw.tmp 
  wigToBigWig ${sample_name}.cpg.meth.bw.tmp $chrSizeFile ${sample_name}.cpg.meth.bw
  rm ${sample_name}.cpg.meth.bw.tmp 
fi

if [ ! -s ${sample_name}.allelic.bw ]; then
  awk '{OFS=\"\\t\"; print \$1,\$2,\$2+1,\$5}' < ${sample_name}.allelic > ${sample_name}.alle.bw.tmp 
  wigToBigWig ${sample_name}.alle.bw.tmp $chrSizeFile ${sample_name}.allelic.bw
  rm ${sample_name}.alle.bw.tmp 
fi

if [ ! -s ${sample_name}.hmr.bb ]; then
  cut -f 1-3 ${sample_name}.hmr > ${sample_name}.hmr.tmp
  bedToBigBed ${sample_name}.hmr.tmp $chrSizeFile ${sample_name}.hmr.bb
  rm  ${sample_name}.hmr.tmp
fi

# if [ ! -s ${sample_name}.pmr.bb ]; then
#   cut -f 1-3 ${sample_name}.pmr > ${sample_name}.pmr.tmp
#   bedToBigBed ${sample_name}.pmr.tmp $chrSizeFile ${sample_name}.pmr.bb
#   rm  ${sample_name}.pmr.tmp
# fi

if [ ! -s ${sample_name}.pmd.bb ]; then
  cut -f 1-3 ${sample_name}.pmd > ${sample_name}.pmd.tmp
  bedToBigBed ${sample_name}.pmd.tmp $chrSizeFile ${sample_name}.pmd.bb
  rm  ${sample_name}.pmd.tmp
fi

$dnmtools_command | grep Version | cut -d ' ' -f 2 | awk '{print \"dnmtools,v\"\$1}' > ${sample_name}.dnmtools.version
";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all DNMTools tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    push( @result_files, "${result_dir}/${sample_name}.levels" );
    #push( @result_files, "${result_dir}/${sample_name}.all.meth" );
    push( @result_files, "${result_dir}/${sample_name}.cpg.meth" );
    push( @result_files, "${result_dir}/${sample_name}.hmr" );
    #no pmr in dnmtools 1.4.1+
    #push( @result_files, "${result_dir}/${sample_name}.pmr" );
    push( @result_files, "${result_dir}/${sample_name}.pmd" );
    push( @result_files, "${result_dir}/${sample_name}.amr" );
    
    
    push( @result_files, "${result_dir}/${sample_name}.bsrate" );
    push( @result_files, "${result_dir}/${sample_name}.uniq.bam.dupstats" );
    
    push( @result_files, "${result_dir}/${sample_name}.dnmtools.version" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
