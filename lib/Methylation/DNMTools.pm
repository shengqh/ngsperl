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
  my $chrDir=$config->{$section}{chr_fasta};
  if ( !defined $chrDir ) {
    die "define ${section}::chr_fasta first";
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
 
    my $result_file          = $sample_name . ".sam";
    my $tag               = get_bam_tag($sample_name);

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $rmlist = "";
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $result_file );

#dnmtools part
    print $pbs "
echo DNMTools=`date`

if [ ! -s ${sample_name}.formatted.sam ]; then
  echo dnmtools format=`date`
  $dnmtools_command format -f abismal -o ${sample_name}.formatted.tmp.bam $sampleFile 
  samtools sort -@ $thread -o ${sample_name}.formatted.sam ${sample_name}.formatted.tmp.bam
  rm -rf ${sample_name}.formatted.tmp.bam
fi

if [ ! -s ${sample_name}.sam ]; then
  echo dnmtools uniq=`date`
  $dnmtools_command uniq -D -S ${sample_name}.sam.dupstats ${sample_name}.formatted.sam ${sample_name}.sam

  if [[ -s ${sample_name}.sam ]]; then
    rm -rf ${sample_name}.formatted.sam
  fi
fi

if [ ! -s ${sample_name}.bsrate ]; then
  echo dnmtools bsrate=`date`
  $dnmtools_command bsrate -c $chrDir ${sample_name}.sam -o ${sample_name}.bsrate
fi


if [ ! -s ${sample_name}.all.meth ]; then
  echo dnmtools counts=`date`
  $dnmtools_command counts -c $chrDir -o ${sample_name}.all.meth ${sample_name}.sam
fi

if [ ! -s ${sample_name}.levels ]; then
  echo dnmtools levels=`date`
  $dnmtools_command levels -o ${sample_name}.levels ${sample_name}.all.meth
fi

if [ ! -s ${sample_name}.meth ]; then
  echo dnmtools sym=`date`
  $dnmtools_command sym -m -o ${sample_name}.meth ${sample_name}.all.meth
fi

if [ ! -s ${sample_name}.hmr ]; then
  echo dnmtools hmr=`date`
  $dnmtools_command hmr -o ${sample_name}.hmr -p ${sample_name}.hmrparams ${sample_name}.meth
fi

if [ ! -s ${sample_name}.pmr ]; then
  echo dnmtools pmr=`date`
  $dnmtools_command hmr -partial -o ${sample_name}.pmr -p ${sample_name}.pmrparams ${sample_name}.meth
fi

if [ ! -s ${sample_name}.pmd ]; then
  echo dnmtools pmd=`date`
  $dnmtools_command pmd -o ${sample_name}.pmd -p ${sample_name}.pmdparams ${sample_name}.meth
fi


CHRM=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX' 'chrY' 'chrMT')
for chr in \"\${CHRM[\@]}\"
do
  if [[ ! -s ${sample_name}.\${chr}.epiread && ! -s ${sample_name}.epiread ]]; then
    echo states_\${chr}=`date`
    p1='/\\b';p2='\\b/p'
    pat=\"\$p1\$chr\$p2\"
    sed -n \$pat ${sample_name}.sam > ${sample_name}.\${chr}.sam
    $dnmtools_command states -c $chrDir -o ${sample_name}.\${chr}.epiread ${sample_name}.\${chr}.sam
    rm ${sample_name}.\${chr}.sam
  fi
done

if [ ! -s ${sample_name}.epiread ]; then
  echo dnmtools states=`date`
  cat ${sample_name}.chr*.epiread > ${sample_name}.epiread
  rm ${sample_name}.chr*.epiread
fi

if [ ! -s ${sample_name}.allelic ]; then
  echo dnmtools allelic=`date`
  $dnmtools_command allelic -c $chrDir -o ${sample_name}.allelic ${sample_name}.epiread
fi
if [ ! -s ${sample_name}.amr ]; then
  echo dnmtools amrfinder=`date`
  $dnmtools_command amrfinder -c $chrDir -o ${sample_name}.amr ${sample_name}.epiread
fi

if [[ -s ${sample_name}.meth ]]; then
  samtools view -h -b -o ${sample_name}.bam ${sample_name}.sam
  samtools index ${sample_name}.bam ${sample_name}.bam.bai
  rm ${sample_name}.sam
fi

";

#make tracks
    print $pbs "
echo DNMTools To Tracks=`date`
if [ ! -s ${sample_name}.read.bw ]; then
  awk '{OFS=\"\\t\"; print \$1,\$2,\$2+1,\$6}' < ${sample_name}.meth > ${sample_name}.read.bw.tmp 
  wigToBigWig ${sample_name}.read.bw.tmp $chrSizeFile ${sample_name}.read.bw
  rm ${sample_name}.read.bw.tmp
fi

if [ ! -s ${sample_name}.meth.bw ]; then
  awk '{OFS=\"\\t\"; print \$1,\$2,\$2+1,\$5}' < ${sample_name}.meth > ${sample_name}.meth.bw.tmp 
  wigToBigWig ${sample_name}.meth.bw.tmp $chrSizeFile ${sample_name}.meth.bw
  rm ${sample_name}.meth.bw.tmp 
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

if [ ! -s ${sample_name}.pmr.bb ]; then
  cut -f 1-3 ${sample_name}.pmr > ${sample_name}.pmr.tmp
  bedToBigBed ${sample_name}.pmr.tmp $chrSizeFile ${sample_name}.pmr.bb
  rm  ${sample_name}.pmr.tmp
fi

if [ ! -s ${sample_name}.pmd.bb ]; then
  cut -f 1-3 ${sample_name}.pmd > ${sample_name}.pmd.tmp
  bedToBigBed ${sample_name}.pmd.tmp $chrSizeFile ${sample_name}.pmd.bb
  rm  ${sample_name}.pmd.tmp
fi

$dnmtools_command | grep Version | cut -d ' ' -f 2 | awk '{print \"dnmtools,v\"\$1}' > ${sample_name}.dnmtools.version
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

  print "!!!shell file $shfile created, you can run this shell file to submit all DNMTools tasks.\n";
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
    push( @result_files, "${result_dir}/${sample_name}.dnmtools.version" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
