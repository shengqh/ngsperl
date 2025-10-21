#!/usr/bin/perl
package Pipeline::MethylationUtils;

use strict;
use warnings;
use CQS::StringUtils;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Data::Dumper;
use List::MoreUtils qw(first_index);
use Hash::Merge     qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(
        add_MethylKitCorr
        add_MethylAgeEstimation
        add_MethylDiffAnalysis
    )
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';


sub add_MethylKitCorr {
  my ( $config, $def, $tasks, $target_dir, $methylkitcorr_task, $methylkitprep_ref, $pipeline ) = @_;
  $config->{$methylkitcorr_task} = {
    class                    => "CQS::UniqueR",
    perform                  => 1,
    target_dir               => "${target_dir}/" . getNextFolderIndex($def) . "MethylKitCorr",
    docker_prefix            => "wgbs_r_",
    rtemplate                => "../Methylation/methylkit_corr.R",
    rReportTemplate          => "../Methylation/methylkit_corr.Rmd;../CQS/reportFunctions.R;../CQS/countTableVisFunctions.R",
    run_rmd_independent      => 1,
    rmd_ext                  => ".methylation.corr.html",
    output_file_ext          => ".methylation.corr.html;.methylKit.version;.filtered.cpg.meth.rds;.CpG.bvalue_hist.png;.CpG.all.PCA.png;.CpG.top10000.PCA.png;.CpG.euclidean_distance.all.MDS.png;.CpG.euclidean_distance.top10000.MDS.png;.CpG.pearson_corr.all.MDS.png;.CpG.pearson_corr.top10000.MDS.png",
    parameterSampleFile1_ref => $methylkitprep_ref,
    parameterSampleFile2     => {
      task_name        => getValue( $def, "task_name" ),
      email            => getValue( $def, "email" ),
      affiliation      => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
      org              => getValue( $def, "genome" ),
      var              => "group",
      control_group    => $def->{control_group},
      mincov           => getValue( $def, "methylation_mincov" ),
      corr_dim1_cutoff => $def->{corr_dim1_cutoff},
      mds_legendPos    => getValue( $def, "mds_legendPos", "bottomleft" ),
      pipeline         => $pipeline,
    },
    parameterSampleFile3_ref => "groups",
    sh_direct                => 1,
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "4",
      "mem"      => "80gb"
    },
  };
  push( @$tasks, $methylkitcorr_task );

  return ($methylkitcorr_task);
} ## end sub add_MethylKitCorr


sub add_MethylAgeEstimation {
  my ( $config, $def, $tasks, $target_dir, $methy_age_task, $methylkitcorr_task ) = @_;
  $config->{$methy_age_task} = {
    class                    => "CQS::UniqueRmd",
    perform                  => 1,
    target_dir               => $target_dir . "/" . getNextFolderIndex($def) . "$methy_age_task",
    report_rmd_file          => "../Methylation/methylation_age.rmd",
    additional_rmd_files     => "../CQS/reportFunctions.R",
    option                   => "",
    parameterSampleFile1_ref => [ $methylkitcorr_task, 'meth.rds$' ],
    parameterSampleFile2     => {
      task_name        => getValue( $def, "task_name" ),
      email            => getValue( $def, "email" ),
      affiliation      => getValue( $def, "affiliation", "CQS/Biostatistics, VUMC" ),
      probe_locus_file => getValue( $def, "probe_locus_file" ),
    },
    parameterSampleFile3     => $def->{"age_dict"},
    suffix                   => ".age",
    output_file_ext          => ".age.html",
    can_result_be_empty_file => 0,
    sh_direct                => 1,
    no_docker                => 1,
    pbs                      => {
      "nodes"    => "1:ppn=1",
      "walltime" => "10",
      "mem"      => "40gb"
    },
  };
  push( @$tasks, $methy_age_task );

  return ($methy_age_task);
} ## end sub add_MethylAgeEstimation


sub add_MethylDiffAnalysis {
  my ( $config, $def, $tasks, $target_dir, $methylkitcorr_task ) = @_;

  my $ncore            = getValue( $def, "MethylKitDiff_ncore", "8" );
  my $annovar_buildver = getValue( $def, "annovar_buildver" );
  my $annovar_db       = getValue( $def, "annovar_db" );
  my $annovar_param    = getValue( $def, "annovar_param" );

  my $test_method = getValue( $def, "methylDiff_test_method", "dss" );

  my $methylkitdiff_task = "MethylKitDiff";
  $config->{$methylkitdiff_task} = {
    class                => "CQS::IndividualR",
    target_dir           => "${target_dir}/" . getNextFolderIndex($def) . "$methylkitdiff_task",
    docker_prefix        => "wgbs_r_",
    rtemplate            => "../Methylation/methylkit_diff.R",
    source_ref           => "pairs",
    parameterSampleFile2 => {
      task_name      => getValue( $def, "task_name" ),
      difference     => getValue( $def, "methylDiff_difference", 25 ),
      qvalue         => getValue( $def, "methylDiff_qvalue",     0.01 ),
      ncore          => $ncore,
      overdispersion => getValue( $def, "methylDiff_overdispersion", "MN" ),     #MN for overdispersion, Chisq-test for no overdispersion
      test_method    => getValue( $def, "methylDiff_test_method",    "dss" ),    #F, Chisq, fast.fisher and dss
      adjust         => getValue( $def, "methylDiff_adjust",         "BH" ),     #SLIM, holm, hochberg, hommel, bonferroni, BH, BY, fdr, none, qvalue
    },
    parameterSampleFile3_ref => "pairs",
    parameterSampleFile4_ref => "groups",
    parameterFile1_ref       => [ $methylkitcorr_task, ".filtered.cpg.meth.rds\$" ],
    output_file_ext          => ".${test_method}.dmcpgs.tsv",
    pbs                      => {
      "nodes"    => "1:ppn=" . $ncore,
      "walltime" => getValue( $def, "MethylKitDiff_walltime", "24" ),
      "mem"      => getValue( $def, "MethylKitDiff_mem",      "80gb" )
    },
  };
  push( @$tasks, $methylkitdiff_task );

  my $methylkitdiffannovar_task = "MethylKitDiffAnnovar";
  $config->{$methylkitdiffannovar_task} = {
    class         => "CQS::ProgramWrapperOneToOne",
    perform       => 1,
    program       => "",
    check_program => 0,
    target_dir    => "${target_dir}/" . getNextFolderIndex($def) . "$methylkitdiffannovar_task",
    option        => "
if [[ ! -s __FILE__ ]]; then
  echo \"Input file __FILE__ not found or empty.\"
  exit 1
fi

perl -lane 'my \$fileColNum=scalar(\@F);my \$fileColPart=join(\"___\",\@F[3..(\$fileColNum-1)]);print \"\$F[0]	\$F[1]	\$F[2]	0	-	\$fileColPart\"' __FILE__ | tail -n +2 > __NAME__.avinput

table_annovar.pl __NAME__.avinput $annovar_db -buildver $annovar_buildver --otherinfo -protocol refGene -operation g --remove --thread 1 --outfile __NAME__.annovar --remove

status=\$?
if [ \$status -ne 0 ]; then
  echo \"Annovar annotation failed with exit code \$status.\"
  exit \$status
fi

echo -e \"chr\tstart\tend\tFunc.refGene\tGene.refGene\tstrand\tpvalue\tqvalue\tmeth.diff\tdirection\" > __NAME__.dmcpgs.annovar.final.tsv
tail -n +2 __NAME__.annovar.${annovar_buildver}_multianno.txt | cut -f1,2,3,6,7,11 | perl -pe 's/___/\\t/g' >> __NAME__.dmcpgs.annovar.final.tsv

rm -rf __NAME__.avinput __NAME__.annovar.${annovar_buildver}_multianno.txt 

version=\$(table_annovar.pl 2>&1 | grep Version | grep -oP '\\d{4}-\\d{2}-\\d{2}')
echo \"Annovar,v\$version\" > __NAME__.annovar.version

",
    docker_prefix         => "annovar_",
    output_file_ext       => ".dmcpgs.annovar.final.tsv;.annovar.version",
    source_ref            => [ $methylkitdiff_task, ".dmcpgs.tsv\$" ],
    output_to_same_folder => 1,
    no_output             => 1,
    sh_direct             => 1,
    pbs                   => {
      "nodes"    => "1:ppn=1",
      "walltime" => "2",
      "mem"      => "10gb"
    },
  };
  push( @$tasks, $methylkitdiffannovar_task );

  my $MethylKitDiffAnnovarGenes_task = "MethylKitDiffAnnovarGenes";
  $config->{$MethylKitDiffAnnovarGenes_task} = {
    class              => "CQS::ProgramWrapperOneToOne",
    perform            => 1,
    target_dir         => "$target_dir/" . getNextFolderIndex($def) . "$MethylKitDiffAnnovarGenes_task",
    interpretor        => "perl",
    program            => "../Methylation/get_gene_names.pl",
    source_ref         => [ $methylkitdiffannovar_task, ".dmcpgs.annovar.final.tsv\$" ],
    output_file_prefix => ".dmcpgs.annovar.final.tsv.genename.txt",
    output_ext         => ".dmcpgs.annovar.final.tsv.genename.txt",
    output_by_file     => 0,
    sh_direct          => 1,
    pbs                => {
      "nodes"    => "1:ppn=1",
      "walltime" => "1",
      "mem"      => "10gb"
    },
  };
  push( @$tasks, $MethylKitDiffAnnovarGenes_task );

  my $webgestalt_task = addWebgestalt( $config, $def, $tasks, $target_dir, $MethylKitDiffAnnovarGenes_task, [ $MethylKitDiffAnnovarGenes_task, ".genename.txt\$" ] );

  my $task_map = {
    methylkitdiff_task             => $methylkitdiff_task,
    methylkitdiffannovar_task      => $methylkitdiffannovar_task,
    MethylKitDiffAnnovarGenes_task => $MethylKitDiffAnnovarGenes_task,
    webgestalt_task                => $webgestalt_task,
  };
  return ($task_map);
} ## end sub add_MethylDiffAnalysis

1;
