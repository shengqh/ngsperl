#!/usr/bin/perl
package Pipeline::WGBS;

use strict;
use warnings;
use List::Util qw(first);
use File::Basename;
use Storable qw(dclone);
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performWGBS)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  fix_task_name($def);

  initDefaultValue( $def, "emailType", "FAIL" );
  initDefaultValue( $def, "cluster",   "slurm" );
  initDefaultValue( $def, "perform_preprocessing", 1 );
  initDefaultValue( $def, "perform_age_estimation", 0 );

  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  $def = initializeDefaultOptions($def);

  my $task_name = $def->{task_name};

  my $email = $def->{email};

  #$def->{perform_cutadapt} = 0;

  my ( $config, $individual, $summary, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);
  my $tasks = [@$individual, @$summary];

  my $target_dir      = $def->{target_dir};

  if ($def->{perform_preprocessing_only}) {

    $config->{sequencetask} = {
      class      => getSequenceTaskClassname($cluster),
      perform    => 1,
      target_dir => "${target_dir}/sequencetask",
      option     => "",
      source     => {
        tasks => $tasks,
      },
      sh_direct => 0,
      cluster   => $cluster,
      pbs       => {
        "nodes"     => "1:ppn=" . $def->{max_thread},
        "walltime"  => $def->{sequencetask_run_time},
        "mem"       => "40gb"
      },
    };
    return $config;
  }

  my $thread = getValue($def, "thread", 4);
  my $abismal_index = getValue($def, "abismal_index");
  my $chr_fasta = getValue($def, "chr_fasta");
  my $chr_size_file = getValue($def, "chr_size_file");
  my $annovar_buildver = getValue($def, "annovar_buildver");
  my $annovar_db = getValue($def, "annovar_db");
  my $annovar_param = getValue($def, "annovar_param");
  my $HOMER_perlFile = getValue($def, "HOMER_perlFile");
  my $addqual_pythonFile = dirname(__FILE__) . "/../Methylation/add_qual.py";
  my $picard = getValue($def, "picard");
  my $interval_list = getValue($def, "interval_list");

  # my $trimgalore_task = "trimgalore";
  # $config->{$trimgalore_task} = {
  #   class => "Trimmer::TrimGalore",
  #   perform => 1,
  #   target_dir => "${target_dir}/" . getNextFolderIndex($def) . "trimgalore",
  #   option => getValue($def, "trimgalore_option", ""),
  #   source_ref => $untrimed_ref,
  #   extension => "_val.fq.gz",
  #   pairend => is_paired_end($def),
  #   do_fastqc => getValue($def, "trimgalore_do_fastqc", 1),
  #   init_command => "",
  #   pbs        => {
  #     "nodes"     => "1:ppn=1",
  #     "walltime"  => "48",
  #     "mem"       => "10gb"
  #   },
  # };
  # push(@$tasks, $trimgalore_task);
  my $methylation_key = "methylation_key";

  my $abismal_task = "abismal";
  $config->{$abismal_task} = {
    class      => "Alignment::Abismal",
    perform    => 1,
    option     => "",
    thread     => $thread,
    target_dir => "${target_dir}/" . getNextFolderIndex($def) . "abismal",
    chr_fasta    => $chr_fasta,
    picard     => $picard,
    interval_list => $interval_list,
    addqual_pythonFile => $addqual_pythonFile,
    abismal_index => $abismal_index,
    source_ref => $source_ref,
    dnmtools_command => getValue($def, "dnmtools_command", "dnmtools"),
    preseq_command => getValue($def, "preseq_command", "preseq"),
    docker_prefix => "dnmtools_",
    pbs        => {
      "nodes"     => "1:ppn=8",
      "walltime"  => getValue($def, "abismal_walltime", "128"),
      "mem"       => "80gb"
    },
  };
  push(@$tasks, $abismal_task);

  my $dnmtools_task = "DNMTools";
  $config->{$dnmtools_task} = {
    class => "Methylation::DNMTools",
    perform       => 1,
    target_dir    => "${target_dir}/" . getNextFolderIndex($def) . "dnmtools",
    option        => "",
    chr_fasta       => $chr_fasta,
    chr_size_file => $chr_size_file,
    source_ref    => ["abismal", ".uniq.bam"],
    dnmtools_command => getValue($def, "dnmtools_command", "dnmtools"),
    docker_prefix => "dnmtools_",
    pbs           => {
      "nodes"     => "1:ppn=8",
      "walltime"  => "72",
      "mem"       => "80gb"
    },
  };
  push(@$tasks, $dnmtools_task);

  my $dnmtoolsannovar_task = "dnmtoolsAnnovar";
  $config->{$dnmtoolsannovar_task} = {
    class      => "Annotation::Annovar",
    perform    => 1,
    target_dir => "${target_dir}/" . getNextFolderIndex($def) . "dnmtoolsAnnovar",
    option     => $annovar_param,
    annovar_db => $annovar_db,
    buildver   => $annovar_buildver,
    perform_splicing => 0,
    docker_prefix => "annovar_",
    isBed      => 1,
    source_ref => [ "DNMTools", ".hmr\$|.amr\$|.pmd\$" ],
    pbs        => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "40gb"
    },
  };
  push(@$tasks, $dnmtoolsannovar_task);

  my $methylkitprep_task = "MethylKitPreparation";
  $config->{$methylkitprep_task} = {
    class                    => "CQS::IndividualR",
    perform                  => 1,
    target_dir               => "${target_dir}/" . getNextFolderIndex($def) . "MethylKitPreparation",
    option                   => "",
    docker_prefix           => "wgbs_r_",
    rtemplate                => "../Methylation/prepare_CpG_input.R",
    parameterSampleFile1_ref   => [ "DNMTools", ".cpg.meth.gz\$" ],
    output_file_ext          => ".CpG.txt.gz",
    sh_direct                => 1,
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "24",
      "mem"       => "100gb"
    },
  };
  push(@$tasks, $methylkitprep_task);

  my $methylkitcorr_task = "MethylKitCorr";
  $config->{$methylkitcorr_task} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => "${target_dir}/" . getNextFolderIndex($def) . "MethylKitCorr",
    #option                   => " --args ${task_name} hg19 group 4 ",
    docker_prefix           => "wgbs_r_",
    rtemplate                => "../Methylation/methylkit_corr.R",
    rReportTemplate          => "../Methylation/methylkit_corr.Rmd;../CQS/reportFunctions.R;../CQS/countTableVisFunctions.R",
    run_rmd_independent => 1,
    rmd_ext => ".methylation.corr.html",
    output_file_ext          => ".methylation.corr.html;.filtered.cpg.meth.rds",
    parameterSampleFile1_ref => $methylkitprep_task,
    parameterSampleFile2 => {
      task_name => $task_name,
      email => getValue($def, "email"),
      affiliation => getValue($def, "affiliation", "CQS/Biostatistics, VUMC"),
      org       => getValue($def, "genome"),
      var       => "group",
      control_group => $def->{control_group},
      mincov     => 4,
      corr_dim1_cutoff => $def->{corr_dim1_cutoff},
      mds_legendPos => getValue($def, "mds_legendPos", "bottomleft"),
    },
    parameterSampleFile3_ref => "groups",
    #parameterFile1 => getValue($def, "meta_file"),
    sh_direct                => 1,
    pbs                      => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "2",
      "mem"       => "80gb"
    },
  };
  push(@$tasks, $methylkitcorr_task);

  if($def->{perform_age_estimation}){
    my $methy_age_task = "MethylKit_age";
    $config->{$methy_age_task} = {
      class              => "CQS::UniqueRmd",
      perform            => 1,
      target_dir         => $target_dir . "/" . getNextFolderIndex($def) . "$methy_age_task",
      report_rmd_file => "../Methylation/methylation_age.rmd",
      additional_rmd_files => "../CQS/reportFunctions.R",
      option => "",
      parameterSampleFile1_ref => [ $methylkitcorr_task, 'meth.rds$' ],
      parameterSampleFile2 => {
        task_name => getValue($def, "task_name"),
        email => getValue($def, "email"),
        affiliation => getValue($def, "affiliation", "CQS/Biostatistics, VUMC"),
        probe_locus_file => getValue($def, "probe_locus_file"),
      },
      parameterSampleFile3 => $def->{"age_dict"},
      suffix => ".age",
      output_file_ext => ".age.html",
      can_result_be_empty_file => 0,
      sh_direct   => 1,
      no_docker => 1,
      pbs => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "10",
        "mem"       => "40gb"
      },
    };
    push( @$tasks, "bamplot_html" );
  }

  my $methylkitdiff_task=undef;
  my $methylkitdiffannovar_task=undef;
  my $MethylKitDiffAnnovarGenes_task=undef;
  my $webgestalt_task=undef;

  if(defined $def->{pairs}){
    my $ncore = getValue($def, "MethylKitDiff_ncore", "8");
    $methylkitdiff_task = "MethylKitDiff";
    $config->{$methylkitdiff_task} = {
      class                    => "Methylation::MethylKitDiff",
      target_dir               => "${target_dir}/" . getNextFolderIndex($def) . "MethylKitDiff",
      docker_prefix           => "wgbs_r_",
      rtemplate                => "../Methylation/methylkit_diff.R",
      source_ref => "pairs",
      parameterSampleFile2 => {
        task_name => $task_name,
        difference => getValue($def, "methylDiff_difference", 25),
        qvalue => getValue($def, "methylDiff_qvalue", 0.01),
        ncore => $ncore,
        overdispersion => getValue($def, "methylDiff_overdispersion", "MN"), #MN for overdispersion, Chisq-test for no overdispersion
        test_method => getValue($def, "methylDiff_test_method", "dss"), #F, Chisq, fast.fisher and dss
        adjust => getValue($def, "methylDiff_adjust", "BH"), #SLIM, holm, hochberg, hommel, bonferroni, BH, BY, fdr, none, qvalue
      },
      parameterSampleFile3_ref => "pairs",
      parameterSampleFile4_ref => "groups",
      parameterFile1_ref => [ "MethylKitCorr", ".filtered.cpg.meth.rds\$" ],
      pbs                      => {
        "nodes"     => "1:ppn=" . $ncore,
        "walltime"  => getValue($def, "MethylKitDiff_walltime", "24"),
        "mem"       => getValue($def, "MethylKitDiff_mem", "80gb")
      },
    };
    push(@$tasks, $methylkitdiff_task);


  #  my $dnmtoolsdiff_task = "DNMToolsDiff";
  #  $config->{$dnmtoolsdiff_task} = {
  #    class      => "Methylation::DNMToolsDiff",
  #    perform    => 1,
  #    target_dir => "${target_dir}/DNMToolsDiff",
  #    option     => "",
  #    source_ref    => "pairs",
  #    methfile_ref  => [ "DNMTools", "^(?!.*?all).*\.meth\$" ],
  #    hmrfile_ref   => [ "DNMTools", ".hmr\$" ],
  #    parameterSampleFile1_ref => "pairs",
  #    parameterSampleFile2_ref => "groups",
  #    minCpG        => 10,
  #    minSigCpG     => 5,
  #    perc_cut      => 0.25,
  #    FDR           => 0.05,
  #    mincov        => 4,
  #    chr_size_file => $chr_size_file,
  #    pbs => {
  #      "nodes"     => "1:ppn=1",
  #      "walltime"  => "12",
  #      "mem"       => "60gb"
  #    },
  #  };
  #  push(@$tasks, "DNMToolsDiff");

    $methylkitdiffannovar_task = "MethylKitDiffAnnovar";
    $config->{$methylkitdiffannovar_task} = {
      class      => "CQS::ProgramWrapperOneToOne",
      perform    => 1,
      program => "",
      check_program => 0,
      target_dir => "${target_dir}/" . getNextFolderIndex($def) . "MethylKitDiffAnnovar",
      option     => "
perl -lane 'my \$fileColNum=scalar(\@F);my \$fileColPart=join(\"  \",\@F[3..(\$fileColNum-1)]);print \"\$F[0]	\$F[1]	\$F[2]	0	-	\$fileColPart\"' __FILE__ | tail -n +2 > __NAME__.avinput

table_annovar.pl __NAME__.avinput $annovar_db -buildver $annovar_buildver --otherinfo -protocol refGene, -operation g --remove --thread 1 --outfile __NAME__.annovar --remove

echo -e \"Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tstrand\tpvalue\tqvalue\tmeth.diff\tdirection\" > __NAME__.dmcpgs.annovar.final.tsv
tail -n +2 __NAME__.annovar.${annovar_buildver}_multianno.txt | perl -pe 's/[ ]+/\\t/g' >> __NAME__.dmcpgs.annovar.final.tsv

rm -rf __NAME__.avinput __NAME__.annovar.${annovar_buildver}_multianno.txt 

",
      docker_prefix => "annovar_",
      output_ext => ".dmcpgs.annovar.final.tsv",
      source_ref => [ "MethylKitDiff", ".dmcpgs\$" ],
      output_to_same_folder => 1,
      sh_direct          => 1,
      pbs        => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "2",
        "mem"       => "10gb"
      },
    };
    push(@$tasks, $methylkitdiffannovar_task);

    $MethylKitDiffAnnovarGenes_task = "MethylKitDiffAnnovarGenes";
    $config->{$MethylKitDiffAnnovarGenes_task} = {
      class              => "CQS::ProgramWrapperOneToOne",
      perform            => 1,
      target_dir         => "$target_dir/" . getNextFolderIndex($def) . "MethylKitDiffAnnovarGenes",
      interpretor => "perl",
      program => "../Methylation/get_gene_names.pl",
      source_ref => ["MethylKitDiffAnnovar", ".dmcpgs.annovar.final.tsv\$" ],
      output_file_prefix => ".dmcpgs.annovar.final.tsv.genename.txt",
      output_ext => ".dmcpgs.annovar.final.tsv.genename.txt",
      output_by_file => 0,
      sh_direct          => 1,
      pbs                => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    };
    push(@$tasks, $MethylKitDiffAnnovarGenes_task);

    $webgestalt_task = addWebgestalt($config, $def, $tasks, $target_dir, $MethylKitDiffAnnovarGenes_task,  [ $MethylKitDiffAnnovarGenes_task, ".genename.txt\$" ]);

  #  my $homer_task = "HOMER_DMR";
  #  $config->{$homer_task} = {
  #    class        => "Homer::FindMotifs",
  #    perform      => 1,
  #    target_dir   => "${target_dir}/HOMER_DMR",
  #    option       => "-nomotif",
  #    homer_genome => "hg19",
  #    source_ref   => [ "DNMToolsDiff", ".DMR.filtered\$" ],
  #    remove_empty_source => 1,
  #    sh_direct    => 0,
  #    pbs          => {
  #      "nodes"     => "1:ppn=1",
  #      "walltime"  => "2",
  #      "mem"       => "40gb"
  #    },
  #  };
  #  push(@$tasks, "HOMER_DMR");
  }

  my @report_files = ();
  my @report_names = ();
  my @copy_files   = ();
  if ( defined $config->{fastqc_raw_summary} ) {
    push( @report_files, "fastqc_raw_summary",                   ".FastQC.baseQuality.tsv.png" );
    push( @report_files, "fastqc_raw_summary",                   ".FastQC.sequenceGC.tsv.png" );
    push( @report_files, "fastqc_raw_summary",                   ".FastQC.adapter.tsv.png" );
    push( @report_names, "fastqc_raw_per_base_sequence_quality", "fastqc_raw_per_sequence_gc_content", "fastqc_raw_adapter_content" );
  }

  if ( defined $config->{fastqc_post_trim_summary} ) {
    push( @report_files, "fastqc_post_trim_summary",                   ".FastQC.baseQuality.tsv.png" );
    push( @report_files, "fastqc_post_trim_summary",                   ".FastQC.sequenceGC.tsv.png" );
    push( @report_files, "fastqc_post_trim_summary",                   ".FastQC.adapter.tsv.png" );
    push( @report_names, "fastqc_post_trim_per_base_sequence_quality", "fastqc_post_trim_per_sequence_gc_content", "fastqc_post_trim_adapter_content" );
  }

  if(( defined $dnmtools_task ) && ( defined $config->{$dnmtools_task} )) {
    push( @copy_files, $dnmtools_task, ".amr\$" );
    push( @copy_files, $dnmtools_task, ".hmr\$" );
    push( @copy_files, $dnmtools_task, ".pmd\$" );
  }
  if(( defined $dnmtoolsannovar_task ) && ( defined $config->{$dnmtoolsannovar_task} )) {
    push( @copy_files, $dnmtoolsannovar_task, ".annovar.final.tsv\$" );
  }
  if (( defined $methylkitcorr_task ) && ( defined $config->{$methylkitcorr_task} )) {
    push( @copy_files, $methylkitcorr_task, ".pdf\$|.png\$|.rds\$" );
  }
  if (( defined $methylkitdiff_task ) && ( defined $config->{$methylkitdiff_task} )) {
    push( @copy_files, $methylkitdiff_task, ".dmcpgs\$" );
  }
  if (( defined $methylkitdiffannovar_task ) && ( defined $config->{$methylkitdiffannovar_task} )) {
    push( @copy_files, $methylkitdiffannovar_task, ".annovar.final.tsv\$" );
  }
  if ( defined $webgestalt_task ) {
    push( @copy_files, $webgestalt_task, "_geneontology_Biological_Process\$" );
    push( @copy_files, $webgestalt_task, "_geneontology_Cellular_Component\$" );
    push( @copy_files, $webgestalt_task, "_geneontology_Molecular_Function\$" );
    push( @copy_files, $webgestalt_task, "_pathway_KEGG\$" );
  }

  $config->{report} = {
    class                      => "CQS::BuildReport",
    perform                    => 1,
    target_dir                 => "$target_dir/report",
    #docker_command             => "singularity exec -c -e -B /home,/gpfs51,/gpfs52,/panfs,/data,/dors,/nobackup,/tmp -H `pwd` /data/cqs/softwares/singularity/wgbs_r.1.1.sif ",
    docker_prefix           => "wgbs_r_",
    report_rmd_file            => "../Pipeline/WGBS.Rmd",
    additional_rmd_files       => "../Pipeline/Pipeline.R;reportFunctions.R",
    parameterSampleFile1_ref   => \@report_files,
    parameterSampleFile1_names => \@report_names,
    parameterSampleFile2 => {
      task_name => $task_name,
      meta_data => "../../" . $task_name . "_meta.tsv",
      abismal_path  => $config->{abismal}{target_dir} . "/result/",
      dnmtools_path => $config->{DNMTools}{target_dir} . "/result/",
      MethylKitCorr_path => $config->{MethylKitCorr}{target_dir} . "/result/",
      MethylKitDiff_path => defined $config->{MethylKitDiff} ? $config->{MethylKitDiff}{target_dir} . "/result/" : undef,
    },
    parameterSampleFile3_ref   => \@copy_files,
    parameterSampleFile4_ref   => $webgestalt_task,
    parameterSampleFile5_ref   => $abismal_task,
    parameterSampleFile6_ref   => $dnmtools_task,
    parameterSampleFile7_ref   => [ $methylkitcorr_task, ".pdf\$|.png\$" ],
    parameterSampleFile8_ref   => $methylkitdiff_task,
    parameterSampleFile9_ref   => $methylkitdiffannovar_task,
    sh_direct                  => 1,
    pbs                        => {
      "nodes"     => "1:ppn=1",
      "walltime"  => "1",
      "mem"       => "40gb"
    },
  };
  push( @$tasks, "report" );

  $config->{"sequencetask"} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      tasks => $tasks,
    },
    sh_direct => 0,
    pbs       => {
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  return($config);
};


sub performWGBS {
  my ( $def, $perform ) = @_;
  if ( !defined $perform ) {
    $perform = 1;
  }

  my $config = getConfig($def);

  if ($perform) {
    saveConfig( $def, $config );

    performConfig($config);
  }

  return $config;
}

sub performWGBSTask {
  my ( $def, $task ) = @_;

  my $config = performWGBS($def, 0);

  performTask( $config, $task );

  return $config;
}


1;

