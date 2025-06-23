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
 
    my $result_file          = "${sample_name}.cpg.meth.gz";
    my $tag               = get_bam_tag($sample_name);

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $result_file );

    print $pbs "
echo DNMTools=`date`

echo dnmtools bsrate=`date`
rm -f ${sample_name}.bsrate.failed ${sample_name}.bsrate

$dnmtools_command bsrate -t $thread -c $chr_fasta -o ${sample_name}.bsrate $sampleFile
status=\$?
if [[ \$status -ne 0 ]]; then
  touch $sample_name.bsrate.failed
  rm ${sample_name}.bsrate
fi

#output cpg only
echo dnmtools counts=`date`
rm -f ${sample_name}.all.meth.failed ${sample_name}.all.meth

$dnmtools_command counts -t $thread -c $chr_fasta -o ${sample_name}.all.meth $sampleFile
status=\$?
if [[ \$status -ne 0 ]]; then
  touch $sample_name.all.meth.failed
  rm ${sample_name}.all.meth
fi

if [[ -s ${sample_name}.all.meth ]]; then
  echo dnmtools levels=`date`
  rm -f ${sample_name}.levels.failed ${sample_name}.levels

  $dnmtools_command levels -o ${sample_name}.levels ${sample_name}.all.meth
  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sample_name.levels.failed
    rm ${sample_name}.levels
  fi

  echo dnmtools sym=`date`
  rm -f ${sample_name}.cpg.meth.failed ${sample_name}.cpg.meth

  $dnmtools_command sym -o ${sample_name}.cpg.meth ${sample_name}.all.meth
  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sample_name.cpg.meth.failed
    rm ${sample_name}.cpg.meth
  else
    rm ${sample_name}.all.meth
  fi
fi

if [[ -s ${sample_name}.cpg.meth ]]; then
  echo cpg.read.bw=`date`
  rm -f ${sample_name}.cpg.read.bw.tmp.failed ${sample_name}.cpg.read.bw.failed

  awk '{OFS=\"\\t\"; print \$1,\$2,\$2+1,\$6}' < ${sample_name}.cpg.meth > ${sample_name}.cpg.read.bw.tmp 
  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch ${sample_name}.cpg.read.bw.tmp.failed
  else
    wigToBigWig ${sample_name}.cpg.read.bw.tmp $chrSizeFile ${sample_name}.cpg.read.bw
    status=\$?
    if [[ \$status -ne 0 ]]; then
      touch ${sample_name}.cpg.read.bw.failed
      rm ${sample_name}.cpg.read.bw
    fi
  fi
  rm -f ${sample_name}.cpg.read.bw.tmp
fi

if [[ -s ${sample_name}.cpg.meth ]]; then
  echo cpg.meth.bw=`date`
  rm -f ${sample_name}.cpg.meth.bw.tmp.failed ${sample_name}.cpg.meth.bw.failed
  
  awk '{OFS=\"\\t\"; print \$1,\$2,\$2+1,\$5}' < ${sample_name}.cpg.meth > ${sample_name}.cpg.meth.bw.tmp 
  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch ${sample_name}.cpg.meth.bw.tmp.failed
  else
    wigToBigWig ${sample_name}.cpg.meth.bw.tmp $chrSizeFile ${sample_name}.cpg.meth.bw
    status=\$?
    if [[ \$status -ne 0 ]]; then
      touch ${sample_name}.cpg.meth.bw.failed
      rm ${sample_name}.cpg.meth.bw
    fi
  fi

  rm -f ${sample_name}.cpg.meth.bw.tmp 
fi

# No pmr in dnmtools 1.4.1
# if [ ! -s ${sample_name}.cpg.pmr ]; then
#   echo dnmtools pmr=`date`
#   $dnmtools_command pmr -partial -o ${sample_name}.cpg.pmr -p ${sample_name}.cpg.pmrparams ${sample_name}.cpg.meth
# fi
#
# if [ ! -s ${sample_name}.cpg.pmr.bb ]; then
#   cut -f 1-3 ${sample_name}.cpg.pmr > ${sample_name}.cpg.pmr.tmp
#   bedToBigBed ${sample_name}.cpg.pmr.tmp $chrSizeFile ${sample_name}.cpg.pmr.bb
#   rm  ${sample_name}.cpg.pmr.tmp
# fi

if [[ -s ${sample_name}.cpg.meth ]]; then
  echo dnmtools hmr=`date`
  rm -f ${sample_name}.cpg.hmr.failed ${sample_name}.cpg.hmr
  
  $dnmtools_command hmr -o ${sample_name}.cpg.hmr -p ${sample_name}.cpg.hmrparams ${sample_name}.cpg.meth
  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sample_name.cpg.hmr.failed
    rm ${sample_name}.cpg.hmr
  else
    echo pg.hmr.bb=`date`
    rm -f ${sample_name}.cpg.hmr.tmp.failed ${sample_name}.cpg.hmr.bb.failed
    
    cut -f 1-3 ${sample_name}.cpg.hmr > ${sample_name}.cpg.hmr.tmp
    status=\$?
    if [[ \$status -ne 0 ]]; then
      touch ${sample_name}.cpg.hmr.tmp.failed
    else
      bedToBigBed ${sample_name}.cpg.hmr.tmp $chrSizeFile ${sample_name}.cpg.hmr.bb
      status=\$?
      if [[ \$status -ne 0 ]]; then
        touch ${sample_name}.cpg.hmr.bb.failed
        rm ${sample_name}.cpg.hmr.bb
      fi
    fi

    rm -f ${sample_name}.cpg.hmr.tmp
  fi
fi


if [[ -s ${sample_name}.cpg.meth ]]; then
  echo dnmtools pmd=`date`
  rm -f ${sample_name}.cpg.pmd.failed

  $dnmtools_command pmd -o ${sample_name}.cpg.pmd -p ${sample_name}.cpg.pmdparams ${sample_name}.cpg.meth
  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sample_name.cpg.pmd.failed
    rm ${sample_name}.cpg.pmd
  else
    echo cpg.pmd.bb=`date`
    rm -f ${sample_name}.cpg.pmd.tmp.failed ${sample_name}.cpg.pmd.bb.failed

    cut -f 1-3 ${sample_name}.cpg.pmd > ${sample_name}.cpg.pmd.tmp
    status=\$?
    if [[ \$status -ne 0 ]]; then
      touch ${sample_name}.cpg.pmd.tmp.failed
    else
      bedToBigBed ${sample_name}.cpg.pmd.tmp $chrSizeFile ${sample_name}.cpg.pmd.bb
      status=\$?
      if [[ \$status -ne 0 ]]; then
        touch ${sample_name}.cpg.pmd.bb.failed
        rm ${sample_name}.cpg.pmd.bb
      fi
    fi

    rm -f ${sample_name}.cpg.pmd.tmp
  fi
fi

echo dnmtools states=`date`
rm -f ${sample_name}.epiread.failed

$dnmtools_command states -t $thread -c $chr_fasta -o ${sample_name}.epiread $sampleFile
status=\$?
if [[ \$status -ne 0 ]]; then
  touch $sample_name.epiread.failed
  rm ${sample_name}.epiread
fi

if [[ -s ${sample_name}.epiread ]]; then
  echo dnmtools allelic=`date`
  rm -f ${sample_name}.allelic.bw.failed ${sample_name}.allelic.bw.tmp.failed

  $dnmtools_command allelic -c $chr_fasta ${sample_name}.epiread | awk '{OFS=\"\\t\"; print \$1,\$2,\$2+1,\$5}' > ${sample_name}.alle.bw.tmp 
  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch ${sample_name}.allelic.bw.tmp.failed
    rm ${sample_name}.alle.bw.tmp
  else
    echo wigToBigWig alle.bw=`date`
    wigToBigWig ${sample_name}.alle.bw.tmp $chrSizeFile ${sample_name}.allelic.bw
    status=\$?
    if [[ \$status -ne 0 ]]; then
      touch ${sample_name}.allelic.bw.failed
      rm ${sample_name}.allelic.bw
    fi
  fi
  rm -f ${sample_name}.alle.bw.tmp
fi

if [[ -s ${sample_name}.epiread ]]; then
  echo dnmtools amrfinder=`date`
  rm -f ${sample_name}.amr.failed

  $dnmtools_command amrfinder -c $chr_fasta -o ${sample_name}.amr ${sample_name}.epiread
  status=\$?
  if [[ \$status -ne 0 ]]; then
    touch $sample_name.amr.failed
    rm ${sample_name}.amr
  fi
fi

rm -f ${sample_name}.epiread

if [[ -s ${sample_name}.cpg.meth ]]; then
  rm -f ${sample_name}.cpg.meth.gz
  gzip ${sample_name}.cpg.meth
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
    push( @result_files, "${result_dir}/${sample_name}.cpg.meth.gz" );
    push( @result_files, "${result_dir}/${sample_name}.cpg.hmr" );
    #no pmr in dnmtools 1.4.1+
    #push( @result_files, "${result_dir}/${sample_name}.cpg.pmr" );
    push( @result_files, "${result_dir}/${sample_name}.cpg.pmd" );
    push( @result_files, "${result_dir}/${sample_name}.amr" );
    
    
    push( @result_files, "${result_dir}/${sample_name}.bsrate" );
    
    push( @result_files, "${result_dir}/${sample_name}.dnmtools.version" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
