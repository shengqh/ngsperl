#!/usr/bin/perl
package Methylation::DNMToolsDiff;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::GroupTask;
use CQS::NGSCommon;
use CQS::StringUtils;
use Data::Dumper;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_dnmtoolsdiff";
  $self->{_group_keys} = ["source"];
  bless $self, $class;
  return $self;
}

sub get_pbs_key {
  my ($self, $config, $section) = @_;
  return "source";
}

sub get_pbs_source {
  my ( $self, $config, $section ) = @_;

  my $method = lc( get_option( $config, $section, "method", "diff" ) );
  my $comparisons = get_raw_files( $config, $section );
  my $groups = $method eq "radmeth" ? get_raw_files( $config, $section, "groups" ) : undef;
  my $pbsFiles = $self->get_pbs_files( $config, $section );

  my $result = {};
  for my $comparison_name ( keys %$pbsFiles ) {
    my ( $ispaired, $group_names ) = get_pair_groups( $comparisons, $comparison_name );
    if ( $method eq "radmeth" ) {
      my @samples = ();
      for my $group_name (@$group_names) {
        die "Cannot find group $group_name for radmeth comparison $comparison_name." if !defined $groups->{$group_name};
        push( @samples, @{ $groups->{$group_name} } );
      }
      $result->{ $pbsFiles->{$comparison_name} } = \@samples;
    }
    else {
      $result->{ $pbsFiles->{$comparison_name} } = $group_names;
    }
  }

  return $result;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $comparisons = get_raw_files( $config, $section);
  my @comparison_names = keys %{$comparisons};

  my $methfiles = get_raw_files( $config, $section, "methfile" );
  my $dnmtools_command = get_option( $config, $section, "dnmtools_command", "dnmtools" );
  my $method = lc( get_option( $config, $section, "method", "diff" ) );
  die "method should be diff or radmeth." if ( $method ne "diff" && $method ne "radmeth" );
  my $hmrfiles = $method eq "diff" ? get_raw_files( $config, $section, "hmrfile" ) : undef;
  my $groups = $method eq "radmeth" ? get_raw_files( $config, $section, "groups" ) : undef;
  my $minCpgInDMR   = get_option( $config, $section, "minCpgInDMR", 10 );
  my $minSigCpgInDMR   = get_option( $config, $section, "minSigCpgInDMR", 5 );
  my $perc_cut   = get_option( $config, $section, "perc_cut", 0.25 );
  my $fdr   = get_option( $config, $section, "FDR", 0.05 );
  my $mincov = get_option( $config, $section, "mincov", 4);
  my $radmeth_factor = get_option( $config, $section, "radmeth_factor", "case" );
  my $radadjust_bins = get_option( $config, $section, "radadjust_bins", "1:200:1" );
  my $radmerge_p = get_option( $config, $section, "radmerge_p", $fdr );
  my $radmeth_option = get_option( $config, $section, "radmeth_option", "" );
  my $radadjust_option = get_option( $config, $section, "radadjust_option", "" );
  my $radmerge_option = get_option( $config, $section, "radmerge_option", "" );
  $thread = get_option( $config, $section, "thread", 1 ) if !defined $thread;

  my $chrSizeFile=$config->{$section}{chr_size_file}; #to make tracks
  if ( !defined $chrSizeFile ) {
    die "define ${section}::chr_size_file first";
  }
  
  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";
  print $sh "cd $pbs_dir\n";

  for my $group_name (@comparison_names) {
    my ( $ispaired, $group_names ) = get_pair_groups( $comparisons, $group_name );
    my @sampleNames = @{ $group_names };
    my $sampleCount = scalar(@sampleNames);

    my $cur_dir = create_directory_or_die( $result_dir . "/$group_name" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
    my $pbs_name = basename($pbs_file);
    my $log = $self->get_log_filename( $log_dir, $group_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir );

    if ( $method eq "radmeth" ) {
      if ( $sampleCount != 2 ) {
        die "Comparison of $group_name should contain two group names for radmeth.";
      }
      if ( !defined $groups || !defined $groups->{ $sampleNames[0] } || !defined $groups->{ $sampleNames[1] } ) {
        die "Define ${section}::groups with $sampleNames[0] and $sampleNames[1] for radmeth.";
      }

      my @controlSamples = @{ $groups->{ $sampleNames[0] } };
      my @treatmentSamples = @{ $groups->{ $sampleNames[1] } };
      if ( scalar(@controlSamples) == 0 || scalar(@treatmentSamples) == 0 ) {
        die "Both groups in comparison $group_name should contain samples for radmeth.";
      }

      my @allSamples = ( @controlSamples, @treatmentSamples );
      my @linkCommands = ();
      my @localMethFiles = ();
      for my $sampleName (@allSamples) {
        if ( !defined $methfiles->{$sampleName} ) {
          die "Cannot find methfile for sample $sampleName in comparison $group_name.";
        }
        my $sourceMethFile = ${ $methfiles->{$sampleName} }[0];
        my $localMethFile = "${sampleName}.cpg.meth.gz";
        push( @localMethFiles, $localMethFile );
        push( @linkCommands, "ln -sf $sourceMethFile $localMethFile" );
      }

      my $dataTable = "${group_name}.radmeth.table.txt";
      my $designFile = "${group_name}.radmeth.design.txt";
      my $radmethFile = "${group_name}.radmeth";
      my $radadjustFile = "${group_name}.radmeth.adjusted";
      my $significantFile = "${group_name}.radmeth.significant";
      my $radmergeFile = "${group_name}.radmeth.dmr";
      my $localMethFilesStr = join( " ", @localMethFiles );
      my $linkCommandsStr = join( "\n", @linkCommands );

      print $pbs "
$linkCommandsStr

if [ ! -s $dataTable ]; then
  echo dnmtools merge radmeth=`date`
  $dnmtools_command merge -t -radmeth -remove .cpg.meth.gz -o $dataTable $localMethFilesStr
fi

if [ ! -s $designFile ]; then
  echo design matrix=`date`
  echo -e \"base\\t$radmeth_factor\" > $designFile
";
      for my $sampleName (@controlSamples) {
        print $pbs "  echo -e \"$sampleName\\t1\\t0\" >> $designFile\n";
      }
      for my $sampleName (@treatmentSamples) {
        print $pbs "  echo -e \"$sampleName\\t1\\t1\" >> $designFile\n";
      }
      print $pbs "fi\n";

      print $pbs "
if [ ! -s $radmethFile ]; then
  echo dnmtools radmeth=`date`
  $dnmtools_command radmeth -t $thread -f $radmeth_factor -o $radmethFile $radmeth_option $designFile $dataTable
fi

if [ ! -s $radadjustFile ]; then
  echo dnmtools radadjust=`date`
  $dnmtools_command radadjust -bins $radadjust_bins -o $radadjustFile $radadjust_option $radmethFile
fi

if [ ! -s $significantFile ]; then
  echo dnmtools radadjust significant CpGs=`date`
  awk '\$7 <= $fdr' $radadjustFile > $significantFile
fi

if [ ! -s $radmergeFile ]; then
  echo dnmtools radmerge=`date`
  $dnmtools_command radmerge -p $radmerge_p -o $radmergeFile $radmerge_option $radadjustFile
fi

echo dnmtools radmerge To Tracks=`date`
if [ ! -s ${radmergeFile}.bb ]; then
  cut -f 1-3 ${radmergeFile} > ${radmergeFile}.tmp
  bedToBigBed ${radmergeFile}.tmp $chrSizeFile ${radmergeFile}.bb
  rm ${radmergeFile}.tmp
fi
";
    }
    else {
      if ( $sampleCount != 2 ) {
        die "SampleFile should be 2 paired samples.";
      }

      my $controlMethFile   = ${$methfiles->{ $sampleNames[0] }}[0];
      my $treatmentMethFile = ${$methfiles->{ $sampleNames[1] }}[0];
      my $methdiffFile      = "${group_name}.methdiff";
      my $controlHmrFile   = ${$hmrfiles->{ $sampleNames[0] }}[0];
      my $treatmentHmrFile = ${$hmrfiles->{ $sampleNames[1] }}[0];
      my $dmrFile1      = ${group_name}."_".basename($controlHmrFile).".DMR";
      my $dmrFile2      = ${group_name}."_".basename($treatmentHmrFile).".DMR";
      my $dmcpgsFile1   = ${group_name}."_".basename($controlHmrFile).".dmcpgs";
      my $dmcpgsFile2   = ${group_name}."_".basename($treatmentHmrFile).".dmcpgs";

      print $pbs "
if [ ! -s $methdiffFile ]; then
  echo dnmtools diff=`date`
  $dnmtools_command diff -o $methdiffFile $controlMethFile $treatmentMethFile
fi

if [[ ! -s $dmrFile1 || ! -s $dmrFile2 ]]; then
  echo dnmtools diff dmr=`date`
  $dnmtools_command dmr $methdiffFile $controlHmrFile $treatmentHmrFile $dmrFile1 $dmrFile2
fi
";
	if ($minCpgInDMR>0 and $minSigCpgInDMR>0) {
		my $dmrFile1Filtered=$dmrFile1.".filtered";
		my $dmrFile2Filtered=$dmrFile2.".filtered";
		
		    print $pbs "
if [ ! -s $dmrFile1Filtered ]; then
  echo methdiff dmr filter=`date`
   awk -F \"[:\\t]\" '\$5 >= $minCpgInDMR && \$6 >= $minSigCpgInDMR {print \$0}' $dmrFile1 > $dmrFile1Filtered
fi
if [ ! -s $dmrFile2Filtered ]; then
  echo methdiff dmr filter=`date`
   awk -F \"[:\\t]\" '\$5 >= $minCpgInDMR && \$6 >= $minSigCpgInDMR {print \$0}' $dmrFile2 > $dmrFile2Filtered
fi
";

#make tracks
    print $pbs "
echo methdiff To Tracks=`date`
if [ ! -s ${dmrFile1Filtered}.bb ]; then
  cut -f 1-3 ${dmrFile1Filtered} > ${dmrFile1Filtered}.tmp
  bedToBigBed  ${dmrFile1Filtered}.tmp $chrSizeFile ${dmrFile1Filtered}.bb
  rm  ${dmrFile1Filtered}.tmp
fi
if [ ! -s ${dmrFile2Filtered}.bb ]; then
  cut -f 1-3 ${dmrFile2Filtered} > ${dmrFile2Filtered}.tmp
  bedToBigBed  ${dmrFile2Filtered}.tmp $chrSizeFile ${dmrFile2Filtered}.bb
  rm  ${dmrFile2Filtered}.tmp
fi

";

# report different methylated CpGs for both directions
    print $pbs "
echo methdiff CpGs=`date`
if [[ ! -s $dmcpgsFile1 &&  ! -s $dmcpgsFile2 ]]; then
  R --vanilla -f /data/cqs/softwares/ngsperl/lib/Methylation/dmcpgs.r --args ${result_dir}/${group_name}/${methdiffFile} $sampleNames[0],$sampleNames[1] $perc_cut $fdr $mincov
fi

";
	}
    }
    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all DNMToolsDiff tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern, $removeEmpty ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $comparisons = get_raw_files( $config, $section );
  my $method = lc( get_option( $config, $section, "method", "diff" ) );
  my $hmrfiles = $method eq "diff" ? get_raw_files( $config, $section, "hmrfile" ) : undef;

  my $result = {};
  for my $group_name ( keys %{$comparisons} ) {
    my @result_files = ();
    my $cur_dir      = $result_dir . "/$group_name";

    if ( $method eq "radmeth" ) {
      push( @result_files, "$cur_dir/${group_name}.radmeth.table.txt" );
      push( @result_files, "$cur_dir/${group_name}.radmeth.design.txt" );
      push( @result_files, "$cur_dir/${group_name}.radmeth" );
      push( @result_files, "$cur_dir/${group_name}.radmeth.adjusted" );
      push( @result_files, "$cur_dir/${group_name}.radmeth.significant" );
      push( @result_files, "$cur_dir/${group_name}.radmeth.dmr" );
      push( @result_files, "$cur_dir/${group_name}.radmeth.dmr.bb" );
      my $filtered = filter_array( \@result_files, $pattern, $removeEmpty );
      if ( scalar(@$filtered) > 0 || !$removeEmpty ) {
        $result->{$group_name} = $filtered;
      }
      next;
    }

    push( @result_files, "$cur_dir/${group_name}.methdiff" );
    my $filtered = filter_array( \@result_files, $pattern, $removeEmpty );
    if ( scalar(@$filtered) > 0 || !$removeEmpty ) {
      $result->{$group_name} = $filtered;
    }
    
    my ( $ispaired, $group_names ) = get_pair_groups( $comparisons, $group_name );
    my @sampleNames = @{ $group_names };

    my $controlHmrFile = defined $hmrfiles->{$sampleNames[0]} ? ${group_name}."_".basename(${$hmrfiles->{$sampleNames[0]}}[0]).".DMR" : ${group_name}."_".$sampleNames[0].".hmr.DMR";
    my $treatmentHmrFile = defined $hmrfiles->{$sampleNames[1]} ? ${group_name}."_".basename(${$hmrfiles->{$sampleNames[1]}}[0]).".DMR" : ${group_name}."_".$sampleNames[1].".hmr.DMR";
    my $controlHmrFileFiltered=$controlHmrFile.".filtered";
    my $treatmentHmrFileFiltered=$treatmentHmrFile.".filtered";
    my $dmcpgsFile1=$controlHmrFile.".dmcpgs";
    my $dmcpgsFile2=$treatmentHmrFile.".dmcpgs";

    @result_files = ();
    push( @result_files, "$cur_dir/${controlHmrFile}" );
    push( @result_files, "$cur_dir/${controlHmrFileFiltered}" );
    push( @result_files, "$cur_dir/${dmcpgsFile1}" );
    $filtered = filter_array( \@result_files, $pattern, $removeEmpty );
    if ( scalar(@$filtered) > 0 || !$removeEmpty ) {
      $result->{$group_name . "_" . $sampleNames[0]} = $filtered;
    }

    @result_files = ();
    push( @result_files, "$cur_dir/${treatmentHmrFile}" );
    push( @result_files, "$cur_dir/${treatmentHmrFileFiltered}" );
    push( @result_files, "$cur_dir/${dmcpgsFile2}" );
    $filtered = filter_array( \@result_files, $pattern, $removeEmpty );
    if ( scalar(@$filtered) > 0 || !$removeEmpty ) {
      $result->{$group_name . "_" . $sampleNames[1]} = $filtered;
    }
  }

  return $result;
}

1;
