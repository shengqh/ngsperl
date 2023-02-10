#!/usr/bin/perl
package Methylation::Methpipe5;

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
  $self->{_suffix} = "_methpipe5";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $selfname = $self->{_name};

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $chrDir=$config->{$section}{chr_dir};
  if ( !defined $chrDir ) {
    die "define ${section}::chr_dir first";
  }
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
 
    my $result_file          = $sample_name . ".meth";
    my $tag               = get_bam_tag($sample_name);

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $rmlist = "";
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $result_file );


#methpip part
    print $pbs "
echo Methpipe5=`date`
if [ ! -s ${sampleFileBase}_f.sam ]; then
   echo format_reads=`date`
   format_reads -o ${sampleFileBase}_f.sam -f abismal $sampleFile
fi
if [ ! -s ${sampleFileBase}_fs.sam ]; then
   echo format_reads=`date`
   samtools sort -o ${sampleFileBase}_fs.sam -O SAM ${sampleFileBase}_f.sam
fi
if [ ! -s ${sampleFileBase}_fsd.sam ]; then
   echo duplicate-remover=`date`
   duplicate-remover -D -s -A -S ${sampleFileBase}.dupstats  -o ${sampleFileBase}_fsd.sam ${sampleFileBase}_fs.sam
fi
if [ ! -s ${sampleFileBase}.bsrate ]; then
bsrate -c $chrDir ${sampleFileBase}_fsd.sam -o ${sampleFileBase}.bsrate
fi
if [ ! -s ${sampleFileBase}.epiread ]; then
methstates -c $chrDir  ${sampleFileBase}_fsd.sam -o ${sampleFileBase}.epiread
fi
if [ ! -s ${sampleFileBase}.all.meth ]; then
methcounts -c $chrDir -o ${sampleFileBase}.all.meth ${sampleFileBase}_fsd.sam
fi
if [ ! -s ${sampleFileBase}.levels ]; then
levels -o ${sampleFileBase}.levels  ${sampleFileBase}.all.meth
fi
if [ ! -s ${sampleFileBase}.meth ]; then
symmetric-cpgs -m -o ${sampleFileBase}.meth ${sampleFileBase}.all.meth
fi
if [ ! -s ${sampleFileBase}.hmr ]; then
hmr -o ${sampleFileBase}.hmr  -p ${sampleFileBase}.hmrparams  ${sampleFileBase}.meth
fi
if [ ! -s ${sampleFileBase}.pmr ]; then
hmr -o ${sampleFileBase}.pmr  -p ${sampleFileBase}.pmrparams -partial ${sampleFileBase}.meth
fi
if [ ! -s ${sampleFileBase}.pmd ]; then
pmd -o ${sampleFileBase}.pmd  -p ${sampleFileBase}.pmdparams ${sampleFileBase}.meth
fi
if [ ! -s ${sampleFileBase}.allelic ]; then
  echo allelicmeth=`date`
  allelicmeth  -o ${sampleFileBase}.allelic -c $chrDir ${sampleFileBase}.epiread
fi
if [ ! -s ${sampleFileBase}.amr ]; then
  echo amrfinder=`date`
  amrfinder -o ${sampleFileBase}.amr  ${sampleFileBase}.epiread  -c $chrDir
fi
";

#make tracks
    print $pbs "
echo Methpipe5 To Tracks=`date`
if [ ! -s ${sampleFileBase}.read.bw ]; then
awk '{OFS=\"\\t\"; print \$1,\$2,\$2+1,\$6}' < ${sampleFileBase}.meth | wigToBigWig /dev/stdin $chrSizeFile ${sampleFileBase}.read.bw
fi
if [ ! -s ${sampleFileBase}.meth.bw ]; then
awk '{OFS=\"\\t\"; print \$1,\$2,\$2+1,\$5}' < ${sampleFileBase}.meth | wigToBigWig /dev/stdin $chrSizeFile ${sampleFileBase}.meth.bw
fi
if [ ! -s ${sampleFileBase}.allelic.bw ]; then
awk '{OFS=\"\\t\"; print \$1,\$2,\$2+1,\$5}' < ${sampleFileBase}.allelic | wigToBigWig /dev/stdin $chrSizeFile ${sampleFileBase}.allelic.bw
fi
if [ ! -s ${sampleFileBase}.hmr.bb ]; then
cut -f 1-3 ${sampleFileBase}.hmr > ${sampleFileBase}.hmr.tmp
bedToBigBed  ${sampleFileBase}.hmr.tmp $chrSizeFile ${sampleFileBase}.hmr.bb && rm  ${sampleFileBase}.hmr.tmp
fi
if [ ! -s ${sampleFileBase}.pmr.bb ]; then
cut -f 1-3 ${sampleFileBase}.pmr > ${sampleFileBase}.pmr.tmp
bedToBigBed  ${sampleFileBase}.pmr.tmp $chrSizeFile ${sampleFileBase}.pmr.bb && rm  ${sampleFileBase}.pmr.tmp
fi
if [ ! -s ${sampleFileBase}.pmd.bb ]; then
cut -f 1-3 ${sampleFileBase}.pmd > ${sampleFileBase}.pmd.tmp
bedToBigBed  ${sampleFileBase}.pmd.tmp $chrSizeFile ${sampleFileBase}.pmd.bb && rm  ${sampleFileBase}.pmd.tmp
fi
";

    if ($rmlist ne "") {
    	    print $pbs "
if [ -s $result_file ]; then
  rm $rmlist 
fi
";
    }
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all Walt tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my $meth_all_file     = "${result_dir}/${sample_name}.all.meth";
    my $meth_file     = "${result_dir}/${sample_name}.meth";
    my $meth_hmr_file     = "${result_dir}/${sample_name}.hmr";
    my $meth_pmr_file     = "${result_dir}/${sample_name}.pmr";
    my $meth_pmd_file     = "${result_dir}/${sample_name}.pmd";
    my $meth_amr_file     = "${result_dir}/${sample_name}.amr";
    my @result_files = ();
    push( @result_files, $meth_all_file );
    push( @result_files, $meth_file );
    push( @result_files, $meth_hmr_file );
    push( @result_files, $meth_pmr_file );
    push( @result_files, $meth_pmd_file );
    push( @result_files, $meth_amr_file );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;