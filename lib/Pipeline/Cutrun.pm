#!/usr/bin/perl
package Pipeline::Cutrun;

#this is a wrapper for cutruntools2, except for peak calling, which supports using igg or input as control.
#https://github.com/fl-yu/CUT-RUNTools-2.0/blob/cf72ca5d0801ccab46ca7cfdb6810628797f4c9b/src/bulk/bulk-pipeline.sh

use strict;
use warnings;
use File::Basename;
use CQS::StringUtils;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::ClassFactory;
use Pipeline::PipelineUtils;
use Pipeline::Preprocession;
use Pipeline::PeakPipelineUtils;
use Data::Dumper;
use Hash::Merge qw( merge );

require Exporter;
our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(performCutrun)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub initializeDefaultOptions {
  my $def = shift;

  initDefaultValue( $def, "sra_to_fastq", 0 );
  initDefaultValue( $def, "remove_duplicates", 1 );

  initDefaultValue( $def, "min_read_length", 25);

  initDefaultValue( $def, "bowtie2_option", "--dovetail --phred33" );

  my $cutrun_type = $def->{cutrun_type};
  if(!defined $cutrun_type){
    die "cutrun_type is not defined. It should be either TF or histone.";
  }else{
    if($cutrun_type ne "TF" && $cutrun_type ne "histone"){
      die "cutrun_type should be either TF or histone.";
    }
  }

  initDefaultValue( $def, "frag_120bp", $cutrun_type eq "TF" ? 1 : 0 );

  initDefaultValue( $def, "perform_macs2_broad", $cutrun_type eq "TF" ? 0 : 1 );
  initDefaultValue( $def, "perform_macs2_narrow", $cutrun_type eq "TF" ? 1 : 0 );
  initDefaultValue( $def, "perform_seacr", 0 );

  initDefaultValue( $def, "macs2_broad_option", "-f BAMPE --broad --broad-cutoff 0.1 -B --SPMR --keep-dup all");
  initDefaultValue( $def, "macs2_narrow_option", "-f BAMPE -q 0.01 -B --SPMR --keep-dup all");

  initDefaultValue( $def, "annotate_nearest_gene", 1 );
  initDefaultValue( $def, "perform_homer", 1 );

  initDefaultValue( $def, "perform_trimmomatic", 1);
  initDefaultValue( $def, "trimmomatic_option", ":2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25");

  if(defined $def->{treatments}){
    $def->{treatments_auto} = 0;
  }else{
    initDefaultValue( $def, "treatments_auto",     1 );
  }

  initDefaultValue( $def, "perform_report", 1 );
  
  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my $cutrun_type = $def->{cutrun_type};
  my $frag_120bp = $def->{frag_120bp};

  my ( $config, $individual_ref, $summary_ref, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);

  init_treatments_design_table($def);

  #merge summary and individual 
  push @$individual_ref, @$summary_ref;
  $summary_ref = $individual_ref;

  #https://github.com/fl-yu/CUT-RUNTools-2.0/blob/master/src/bulk/bulk-pipeline.sh

  my $genome_sequence = getValue($def, "fasta_file");

  my $bamTLEN_script = dirname(__FILE__) . "/../Alignment/bamTemplateLength.py";

  my $bowtie2_option = getValue( $def, "bowtie2_option" );
  my $bowtie2_index = getValue( $def, "bowtie2_index" );

  my $cutruntools2_path = getValue($def, "cutruntools2_path");
  my $extratoolsbin = $cutruntools2_path . "/install";
  my $awk1 = "$extratoolsbin/filter_below.awk";
  my $picard_jar = "$extratoolsbin/picard-2.8.0.jar";

  my $remove_duplicates = getValue( $def, "remove_duplicates", 1 );
  my $sorted_bam = $remove_duplicates ? "__NAME__.rmdup.bam" : "__NAME__.sorted.bam";

  my $rmdup_pbs;
  if($remove_duplicates){
    $rmdup_pbs = "
echo MarkDuplicates=`date`
java -jar $picard_jar MarkDuplicates \\
  INPUT=__NAME__.sorted.bam \\
  OUTPUT=__NAME__.rmdup.bam \\
  VALIDATION_STRINGENCY=SILENT \\
  TMP_DIR=. \\
  REMOVE_DUPLICATES=true \\
  METRICS_FILE=__NAME__.rmdup.metrics.txt

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.rmdup.succeed
  rm __NAME__.sorted.bam __NAME__.sorted.bam.bai
else
  echo \$status > __NAME__.rmdup.failed
  rm -f __NAME__.rmdup.bam
  exit \$status
fi

echo index=`date`
samtools index __NAME__.rmdup.bam
";
  }else{
    $rmdup_pbs = "";
  }

  my $filtered_suffix = $remove_duplicates ? ".rmdup" : "";
  my $filtered_name = "__NAME__" . $filtered_suffix;
  my $frag_120bp_pbs = "";

  if($frag_120bp){
    $filtered_suffix = $remove_duplicates ? ".rmdup.120bp" : ".120bp";
    $filtered_name = "__NAME__" . $filtered_suffix;
    $frag_120bp_pbs = "
echo filter_120bp=`date`
samtools view -h $sorted_bam | LC_ALL=C awk -f $awk1 | samtools view -Sb - > $filtered_name.bam

status=\$?
if [[ \$status -eq 0 ]]; then
  touch $filtered_name.succeed
  rm $sorted_bam
else
  echo \$status > $filtered_name.failed
  rm -f $filtered_name.bam
  exit \$status
fi

echo index=`date`
samtools index $filtered_name.bam
";  
  }

  my $remove_chromosome = getValue( $def, "remove_chromosome", "chrM" );
  if ( $remove_chromosome !~ /^\s*$/ ) {
    $remove_chromosome = "| grep -v '^$remove_chromosome'";
  }

  my $keep_chromosome = getValue( $def, "keep_chromosome", "chr" );
  if ( $keep_chromosome !~ /^\s*$/ ) {
    $keep_chromosome = "| grep '^$keep_chromosome'";
  }

  my $bowtie2_task = "bowtie2";
  $config->{ $bowtie2_task } = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "${target_dir}/" . getNextFolderIndex($def) . "$bowtie2_task",
    interpretor => "",
    program => "",
    check_program => 0,
    option => "
if [[ ! -s $awk1 ]]; then
  echo \"awk script not found : $awk1\"
  exit 1
fi

rm -f *.failed *.succeed

echo bowtie2=`date` 
bowtie2 -p 8 $bowtie2_option \\
  -x $bowtie2_index \\
  -1 __FILE__ \\
  --sam-RG ID:__NAME__ \\
  --sam-RG LB:__NAME__ \\
  --sam-RG SM:__NAME__ \\
  --sam-RG PL:ILLUMINA \\
  --no-unal \\
  2> __NAME__.log | samtools view -f 3 -F 4 -F 8 -bS - > __NAME__.unsorted.bam

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.bowtie2.succeed
else
  echo \$status > __NAME__.bowtie2.failed
  rm -f __NAME__.unsorted.bam
  exit \$status
fi

echo sort=`date` 
samtools sort -@ 8 -o __NAME__.sorted.bam -T __NAME__ __NAME__.unsorted.bam

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.sort.succeed
  rm __NAME__.unsorted.bam
else
  echo \$status > __NAME__.sort.failed
  rm -f __NAME__.sorted.bam
  exit \$status
fi

echo index=`date`
samtools index __NAME__.sorted.bam

$rmdup_pbs

$frag_120bp_pbs

echo clean_chromosome=`date`
samtools idxstats $filtered_name.bam | cut -f 1 $remove_chromosome $keep_chromosome | xargs samtools view -b -o $filtered_name.clean.bam $filtered_name.bam

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.clean.succeed
  rm $filtered_name.bam $filtered_name.bam.bai
else
  echo \$status > __NAME__.clean.failed
  rm -f $filtered_name.clean.bam
  exit \$status
fi

echo index=`date`
samtools index $filtered_name.clean.bam

echo idxstats=`date`
samtools idxstats $filtered_name.clean.bam $remove_chromosome $keep_chromosome > $filtered_name.clean.bam.chromosome.count

echo flagstat=`date`
samtools flagstat $filtered_name.clean.bam > $filtered_name.clean.bam.stat 

echo tlen=`date`
python $bamTLEN_script -i $filtered_name.clean.bam -o $filtered_name.clean.tlen.txt

bowtie2 --version | grep -a bowtie2 | grep -a version | cut -d ' ' -f3 | awk '{print \"bowtie2,v\"\$1}' > __NAME__.bowtie2.version

",
    source_ref            => $source_ref,
    source_join_delimiter => " -2 ",
    output_to_same_folder => 0,
    no_output => 1,
    output_file_ext => "$filtered_suffix.clean.bam,$filtered_suffix.clean.bam.chromosome.count,$filtered_suffix.clean.bam.stat,.log,.bowtie2.version,$filtered_suffix.clean.tlen.txt,$filtered_suffix.clean.tlen.txt.png",
    sh_direct             => 0,
    docker_prefix => "cutruntools2_",
    pbs                   => {
      "nodes"    => "1:ppn=8",
      "walltime" => getValue($def, "bowtie2_walltime", "24"),
      "mem"      => getValue($def, "bowtie2_mem", "40gb")
    },
  };
  my $bam_ref = [$bowtie2_task, ".bam\$" ];
  push @$summary_ref, ( $bowtie2_task );

  my $summary_task = "${bowtie2_task}_summary";
  add_alignment_summary($config, $def, $summary_ref, $target_dir, $summary_task, "countTableVisFunctions.R;../Alignment/AlignmentUtils.r;../Alignment/Bowtie2Summary.r", ".reads.csv;.reads.png;.chromosome.csv;.chromosome.png", [ "bowtie2", ".log" ], ["bowtie2", ".chromosome.count"] );

  my $bam_files = get_result_file($config, $bowtie2_task, ".bam\$");
  my $treatment_samples = do_get_group_samplefile_map(getValue($def, "treatments"), $bam_files);
  
  my $control_samples = undef;
  if(defined $def->{controls}){
    $control_samples = do_get_group_samplefile_map(getValue($def, "controls"), $bam_files);
  }

  my $peak_calling_task = undef;

  if($def->{perform_macs2_broad}){
    my $macs2_genome = getValue($def, "macs2_genome");
    my $macs2_broad_task = "macs2_broad";
    $config->{ $macs2_broad_task } = {
      class                 => "Chipseq::MACS2Callpeak",
      perform               => 1,
      target_dir            => "${target_dir}/" . getNextFolderIndex($def) . "$macs2_broad_task",
      option     => "-g $macs2_genome --outdir . ". getValue($def, "macs2_broad_option"),
      source_ref => $bam_ref,
      groups     => $def->{"treatments"},
      controls   => $def->{"controls"},
      sh_direct  => 0,
      docker_prefix => "cutruntools2_",
      pbs        => {
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    push @$summary_ref, ( $macs2_broad_task );
    $peak_calling_task = $macs2_broad_task;

    add_peak_count($config, $def, $summary_ref, $target_dir, $macs2_broad_task . "_count", $macs2_broad_task);
  }
   
  if($def->{perform_macs2_narrow}){
    my $macs2_genome = getValue($def, "macs2_genome");
    my $macs2_narrow_task = "macs2_narrow";
    $config->{ $macs2_narrow_task } = {
      class                 => "Chipseq::MACS2Callpeak",
      perform               => 1,
      target_dir            => "${target_dir}/" . getNextFolderIndex($def) . "$macs2_narrow_task",
      option     => "-g $macs2_genome --outdir . ". getValue($def, "macs2_narrow_option"),
      source_ref => $bam_ref,
      groups     => $def->{"treatments"},
      controls   => $def->{"controls"},
      sh_direct  => 0,
      docker_prefix => "cutruntools2_",
      pbs        => {
        "nodes"    => "1:ppn=1",
        "walltime" => "72",
        "mem"      => "40gb"
      },
    };
    push @$summary_ref, ( $macs2_narrow_task );
    $peak_calling_task = $macs2_narrow_task;

    add_peak_count($config, $def, $summary_ref, $target_dir, $macs2_narrow_task . "_count", $macs2_narrow_task);
  }
 
  if($def->{perform_seacr}){
    my $has_control = defined $def->{controls};
    my $seacr_task;

    if($has_control){
      my $genomecov_task = "genomecov";
      my $genome_sizes = getValue($def, "genome_sizes");
      $config->{ $genomecov_task } = {
        class                 => "CQS::ProgramWrapperOneToOne",
        perform               => 1,
        target_dir            => "${target_dir}/" . getNextFolderIndex($def) . "$genomecov_task",
        interpretor => "",
        program => "",
        check_program => 0,
        option     => "
  bedtools bamtobed -bedpe -i __FILE__ > __NAME__.bed

  awk '\$1==\$4 && \$6-\$2 < 1000 {print \$0}' __NAME__.bed > __NAME__.clean.bed
  rm __NAME__.bed

  cut -f 1,2,6 __NAME__.clean.bed | sort -k1,1 -k2,2n -k3,3n > __NAME__.fragments.bed
  rm __NAME__.clean.bed

  bedtools genomecov -bg -i __NAME__.fragments.bed -g $genome_sizes > __NAME__.fragments.bedgraph
  rm __NAME__.fragments.bed

  ",
        source_ref => $bam_ref,
        sh_direct  => 0,
        no_output => 1,
        output_file_ext => ".fragments.bedgraph",
        docker_prefix => "cutruntools2_",
        pbs        => {
          "nodes"    => "1:ppn=1",
          "walltime" => "24",
          "mem"      => "10gb"
        },
      };
      push @$summary_ref, ( $genomecov_task );

      $config->{genomecov_treatments} = {
        class => "CQS::GroupPickTask",
        source_ref => $genomecov_task,
        groups => getValue($def, "treatments"),
      };

      $config->{genomecov_controls} = {
        class => "CQS::GroupPickTask",
        source_ref => $genomecov_task,
        groups => getValue($def, "controls"),
      };

      my $genomecov_seacr_task = "${genomecov_task}_seacr";
      my $seacr_path = getValue($def, "seacr_path");
      $config->{$genomecov_seacr_task} = {
        class => "CQS::ProgramWrapperOneToOne",
        target_dir => "${target_dir}/$genomecov_seacr_task",
        interpretor => "",
        program => "",
        check_program => 0,
        docker_prefix => "cutruntools2_",
        option => "
  echo seacr_relaxed=`date` 
  bash $seacr_path __FILE__ __FILE2__ norm relaxed __NAME__

  echo sort-bed_relaxed=`date` 
  sort-bed __NAME__.relaxed.bed > __NAME__.relaxed.sort.bed
  rm __NAME__.relaxed.bed

  echo get_summits_seacr_relaxed=`date` 
  python $extratoolsbin/get_summits_seacr.py __NAME__.relaxed.sort.bed | sort-bed - > __NAME__.relaxed.sort.summits.bed

  echo seacr_stringent=`date` 
  bash $seacr_path __FILE__ __FILE2__ norm stringent __NAME__

  echo sort-bed_stringent=`date` 
  sort-bed __NAME__.stringent.bed > __NAME__.stringent.sort.bed
  rm __NAME__.stringent.bed

  echo get_summits_seacr_stringent=`date` 
  python $extratoolsbin/get_summits_seacr.py __NAME__.stringent.sort.bed | sort-bed - > __NAME__.stringent.sort.summits.bed

  ",
        source_ref => "genomecov_treatments",
        parameterSampleFile2_ref => "genomecov_controls",
        output_file_ext => ".relaxed.sort.bed,.relaxed.sort.summits.bed,,.stringent.sort.bed,.stringent.sort.summits.bed",
        output_to_same_folder => 0,
        can_result_be_empty_file => 0,
        no_output => 1,
        sh_direct => 0,
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      },
      push @$summary_ref, ($genomecov_seacr_task);

      add_peak_count($config, $def, $summary_ref, $target_dir, $genomecov_seacr_task . "_count", [$genomecov_seacr_task, ".relaxed.sort.bed"]);

      $seacr_task = $genomecov_seacr_task;
      $peak_calling_task = $seacr_task;
    }
  }

  if ( getValue( $def, "perform_homer_find_peaks" ) ) {
    my $homer_makeTagDirectory = "homer_01_makeTagDirectory";
    my $homer_makeTagDirectory_option = getValue($def, "homer_makeTagDirectory_option", "");
    $config->{ $homer_makeTagDirectory } = {
      class                 => "CQS::ProgramWrapperOneToOne",
      perform               => 1,
      target_dir            => "${target_dir}/" . getNextFolderIndex($def) . "$homer_makeTagDirectory",
      interpretor => "",
      program => "",
      check_program => 0,
      option     => "
rm -rf __NAME__.failed __NAME__.makeTagDirectory.failed __NAME__.makeTagDirectory.succeed

makeTagDirectory __NAME__ $homer_makeTagDirectory_option __FILE__

status=\$?
if [[ \$status -eq 0 && -s __NAME__/tagAutocorrelation.txt ]]; then
  touch __NAME__.makeTagDirectory.succeed
else
  echo \$status > __NAME__.makeTagDirectory.failed
  mv __NAME__ __NAME__.failed
fi
",
      source_ref => $bam_ref,
      sh_direct  => 0,
      no_output => 1,
      output_file_ext => "__NAME__",
      pbs        => {
        "nodes"    => "1:ppn=1",
        "walltime" => "24",
        "mem"      => "10gb"
      },
    };
    push @$summary_ref, ( $homer_makeTagDirectory );

    $config->{tag_treatments} = {
      class => "CQS::GroupPickTask",
      source_ref => $homer_makeTagDirectory,
      groups => getValue($def, "treatments"),
    };

    my $control_def = undef;
    my $control_option = "";
    if(defined $def->{controls}){
      $config->{tag_controls} = {
        class => "CQS::GroupPickTask",
        source_ref => $homer_makeTagDirectory,
        groups => getValue($def, "controls"),
      };
      $control_option = "-i __FILE2__";
      $control_def = "tag_controls";
    }

    my $homer_genome = getValue($def, "homer_genome");
    my $rename_2_py = dirname(__FILE__) . "/../Homer/rename_annotate_peaks_path.py";

    my $homer_findPeaks = "homer_02_findPeaks";
    my $default_findPeaks_option = $cutrun_type eq "TF" ? "-style factor -center -size 200 -tbp 0" : "-style histone -tbp 0 -F 2 -P 0.01";
    my $homer_findPeaks_option = getValue($def, "homer_findPeaks_option", $default_findPeaks_option);
    $config->{$homer_findPeaks} = {
      class => "CQS::ProgramWrapperOneToOne",
      target_dir => "${target_dir}/$homer_findPeaks",
      interpretor => "",
      program => "",
      check_program => 0,
      option => "
if [[ ! -s __FILE__ ]]; then
  touch __NAME__.findPeaks.failed
  rm -f __NAME__.findPeaks.succeed
  exit 1
fi

findPeaks __FILE__ $homer_findPeaks_option $control_option -o tmp.__NAME__.peaks.txt

status=\$?
if [[ \$status -eq 0 ]]; then
  rm -f __NAME__.findPeaks.failed
  touch __NAME__.findPeaks.succeed
  mv tmp.__NAME__.peaks.txt __NAME__.peaks.txt

  grep -v '^#' __NAME__.peaks.txt | awk -v OFS='\\t' -F'\\t' '{ print \$2,\$3,\$4,\$1,\$6,\$5}'  > __NAME__.peaks.bed

  echo annotatePeaks_raw=`date`
  annotatePeaks.pl __NAME__.peaks.txt \\
    $homer_genome -raw -d \\
    __FILE__ > __NAME__.peaks.raw.txt

  echo annotatePeaks_rpkm=`date`
  annotatePeaks.pl __NAME__.peaks.txt \\
    $homer_genome -fpkm -d \\
    __FILE__ > __NAME__.peaks.fpkm.txt

  python $rename_2_py __NAME__.peaks.raw.txt __NAME__.peaks.fpkm.txt  
else
  echo \$status > __NAME__.findPeaks.failed
  rm -f tmp.__NAME__.peaks.txt __NAME__.findPeaks.succeed
fi

",
      source_ref => "tag_treatments",
      parameterSampleFile2_ref => $control_def,
      output_file_ext => ".peaks.txt,.peaks.bed,.peaks.raw.txt,.peaks.fpkm.txt",
      output_to_same_folder => 1,
      can_result_be_empty_file => 0,
      no_output => 1,
      sh_direct => 1,
      pbs => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    },
    push @$summary_ref, ($homer_findPeaks);

    my $call_task;
    my $count_task;
    my $homer_motif;

    if(getValue($def, "perform_homer_mergePeaks", 1)){
      my $rename_1_py = dirname(__FILE__) . "/../Homer/rename_find_peaks_path.py";
      my $homer_mergePeaks = "homer_03_mergePeaks";
      my $homer_mergePeaks_option = getValue($def, "homer_mergePeaks_option", "-d given");

      my $homer_mergePeaks_groups;
      if($def->{homer_mergePeaks_groups_pattern}){
        my $gpattern = $def->{homer_mergePeaks_groups_pattern};
        $homer_mergePeaks_groups = {};
        my $treatments = $def->{treatments};
        for my $samplename (keys %$treatments){
          my $groupname = $samplename;
          if($samplename =~ /$gpattern/){
            $groupname = $1;
            if(defined $2){
              $groupname = $groupname . $2;
            }
            if(defined $3){
              $groupname = $groupname . $3;
            }
            #print($groupname . " : " . $samplename . "\n");
            if (not defined $homer_mergePeaks_groups->{$groupname}){
              $homer_mergePeaks_groups->{$groupname} = [$samplename];
            }
            else{
              my $samples = $homer_mergePeaks_groups->{$groupname};
              push (@$samples, $samplename);
            }
          }
        }
      }elsif($def->{homer_mergePeaks_groups}){
        $homer_mergePeaks_groups = $def->{homer_mergePeaks_groups};
      }else{
        $homer_mergePeaks_groups = {};
        my $treatments = $def->{treatments};
        $homer_mergePeaks_groups->{$def->{task_name}} = [sort keys %$treatments];
      }
      $config->{homer_merge_groups} = {
        class => "CQS::GroupPickTask",
        source_ref => [$homer_findPeaks, ".peaks.txt"],
        groups => $homer_mergePeaks_groups,
        return_all => 1,
      };

      my $homer_merge_peaks_tag_groups = {};
      my $treatments = $def->{treatments};
      for my $group_name (keys %$homer_mergePeaks_groups){
        my $treatment_names = $homer_mergePeaks_groups->{$group_name};
        my $sample_names = [];
        for my $treatment_name (@$treatment_names){
          my $sample_name = $treatments->{$treatment_name}->[0];
          push @$sample_names, $sample_name;
        }
        $homer_merge_peaks_tag_groups->{$group_name} = $sample_names;
      }
      $config->{homer_merge_tag_groups} = {
        class => "CQS::GroupPickTask",
        source_ref => $homer_makeTagDirectory,
        groups => $homer_merge_peaks_tag_groups,
        return_all => 1,
      };

      $config->{$homer_mergePeaks} = {
        class => "CQS::ProgramWrapperOneToOne",
        target_dir => "${target_dir}/$homer_mergePeaks",
        interpretor => "",
        program => "",
        check_program => 0,
        option => "
  echo mergePeaks=`date`
  mergePeaks $homer_mergePeaks_option \\
    __FILE__ \\
    -matrix __NAME__ \\
    -venn __NAME__.venn.txt > tmp.__NAME__.all_dGiven.peaks.txt

  if [[ \$(wc -l <tmp.__NAME__.all_dGiven.peaks.txt) -ge 2 ]]; then
    export status=0
  else
    export status=1
  fi

  if [[ \$status -ne 0 ]]; then
    echo \$status > __NAME__.mergePeaks.failed
    rm -f tmp.__NAME__.all_dGiven.peaks.txt __NAME__.mergePeaks.succeed
  else
    touch __NAME__.mergePeaks.succeed
    rm -f __NAME__.mergePeaks.failed
    mv tmp.__NAME__.all_dGiven.peaks.txt __NAME__.all_dGiven.peaks.txt

    echo rename_names=`date`
    python $rename_1_py \\
      __NAME__.all_dGiven.peaks.txt \\
      __NAME__.venn.txt \\
      __NAME__.count.matrix.txt \\
      __NAME__.logPvalue.matrix.txt \\
      __NAME__.logRatio.matrix.txt

    awk -v OFS='\\t' -F'\\t' '{ print \$2,\$3,\$4,\$1,\$6,\$5}' __NAME__.all_dGiven.peaks.txt > __NAME__.all_dGiven.peaks.bed

    echo annotatePeaks_raw=`date`
    annotatePeaks.pl __NAME__.all_dGiven.peaks.txt \\
      $homer_genome -raw -d \\
      __FILE2__ > __NAME__.all_dGiven.peaks.raw.txt

    echo annotatePeaks_rpkm=`date`
    annotatePeaks.pl __NAME__.all_dGiven.peaks.txt \\
      $homer_genome -fpkm -d \\
      __FILE2__ > __NAME__.all_dGiven.peaks.fpkm.txt

    python $rename_2_py __NAME__.all_dGiven.peaks.raw.txt __NAME__.all_dGiven.peaks.fpkm.txt  
  fi
  ",
        source_ref => "homer_merge_groups",
        source_join_delimiter => " \\\n  ",
        source_type => "array",
        parameterSampleFile2_ref => "homer_merge_tag_groups",
        parameterSampleFile2_join_delimiter => " \\\n  ",
        parameterSampleFile2_type => "array",
        output_arg => "",
        output_file_ext => ".all_dGiven.peaks.txt",
        output_other_ext => ".all_dGiven.peaks.bed,.all_dGiven.peaks.fpkm.txt",
        output_to_same_folder => 1,
        can_result_be_empty_file => 0,
        no_output => 1,
        sh_direct => 1,
        pbs => {
          "nodes"     => "1:ppn=1",
          "walltime"  => "1",
          "mem"       => "10gb"
        },
      },
      push @$summary_ref, ($homer_mergePeaks);

      $count_task = $homer_mergePeaks . "_count";
      add_peak_count($config, $def, $summary_ref, $target_dir, $count_task, $homer_mergePeaks);
      
      $homer_motif = "homer_04_motif";
      $call_task = $homer_mergePeaks;
    }else{
      $count_task = $homer_findPeaks . "_count";
      add_peak_count($config, $def, $summary_ref, $target_dir, $count_task, $homer_findPeaks);

      $homer_motif = "homer_03_motif";
      $call_task = $homer_findPeaks;
    }

    addHomerAnnotation( $config, $def, $summary_ref, $target_dir, $call_task, ".bed", $homer_motif );
  }

  if(defined $peak_calling_task){
    if ( getValue( $def, "perform_homer" ) ) {
      my $homer = addHomerAnnotation( $config, $def, $summary_ref, $target_dir, $peak_calling_task, ".bed", $peak_calling_task . "_homer_motif" );
    }

    if($def->{annotate_nearest_gene}){
      my $gene_bed = getValue($def, "gene_bed");
      add_nearest_gene($config, $def, $summary_ref, $target_dir, $peak_calling_task, ".bed", $gene_bed);
    }
  }

  if ($def->{perform_bamplot}){
    defined $def->{dataset_name} or die "Define dataset_name for bamplot first!";
    if ( not defined $def->{bamplot_gff} ) {
      defined $def->{gene_names} or die "Define gene_names for bamplot first, seperate by blank space!";
      defined $def->{add_chr}    or die "Define add_chr for bamplot first, check your genome sequence!";
    }
    add_bamplot($config, $def, $summary_ref, $target_dir, $bam_ref);
  }

  if ( getValue( $def, "perform_report" ) ) {
    my $task_dic = {};

    if($config->{fastqc_raw_summary}){
      $task_dic->{fastqc_raw_summary} = $config->{fastqc_raw_summary}{target_dir};
    }

    if($config->{trimmomatic_fastqc_summary}){
      $task_dic->{trimmomatic_fastqc_summary} = $config->{trimmomatic_fastqc_summary}{target_dir};
    }

    if($config->{bowtie2_summary}){
      $task_dic->{bowtie2_summary} = $config->{bowtie2_summary}{target_dir};
    }

    if($config->{macs2_broad_count}){
      $task_dic->{macs2_broad_count} = $config->{macs2_broad_count}{target_dir};
    }

    if($config->{macs2_broad_homer_motif}){
      $task_dic->{macs2_broad_homer_motif} = $config->{macs2_broad_homer_motif}{target_dir};
    }

    if($config->{macs2_narrow_count}){
      $task_dic->{macs2_narrow_count} = $config->{macs2_narrow_count}{target_dir};
    }

    if($config->{macs2_narrow_homer_motif}){
      $task_dic->{macs2_narrow_homer_motif} = $config->{macs2_narrow_homer_motif}{target_dir};
    }

    if($config->{homer_02_findPeaks_count}){
      $task_dic->{homer_02_findPeaks_count} = $config->{homer_02_findPeaks_count}{target_dir};
    }

    if($config->{homer_03_mergePeaks_count}){
      $task_dic->{homer_03_mergePeaks_count} = $config->{homer_03_mergePeaks_count}{target_dir};
    }

    if($config->{homer_03_motif}){
      $task_dic->{homer_03_motif} = $config->{homer_03_motif}{target_dir};
    }

    if($config->{homer_04_motif}){
      $task_dic->{homer_04_motif} = $config->{homer_04_motif}{target_dir};
    }

    if($config->{bamplot}){
      $task_dic->{bamplot} = $config->{bamplot}{target_dir};
    }

    $config->{report} = {
      class                      => "CQS::BuildReport",
      perform                    => 1,
      target_dir                 => $target_dir . "/" . getNextFolderIndex($def) . "report",
      report_rmd_file            => "../Pipeline/Cutrun.rmd",
      additional_rmd_files       => "reportFunctions.R;../Pipeline/Pipeline.R",
      docker_prefix              => "report_",
      parameterSampleFile1 => {
        task_name => getValue($def, "task_name"),
        email => getValue($def, "email"),
        affiliation => getValue($def, "affiliation", ""),
        cutrun_type => $cutrun_type,
        frag_120bp => $frag_120bp,
      },
      parameterSampleFile2       => $task_dic,
      parameterSampleFile4_ref   => [keys %$task_dic],
      parameterSampleFile5_ref   => ["bowtie2", '.png$'],
      sh_direct                  => 1,
      pbs                        => {
        "nodes"     => "1:ppn=1",
        "walltime"  => "1",
        "mem"       => "10gb"
      },
    };
    push( @$summary_ref, "report" );
  }

  $config->{"sequencetask"} = {
    class      => getSequenceTaskClassname($cluster),
    perform    => 1,
    target_dir => "${target_dir}/sequencetask",
    option     => "",
    source     => {
      tasks => $summary_ref,
    },
    sh_direct => 0,
    pbs       => {
      "nodes"    => "1:ppn=8",
      "walltime" => "72",
      "mem"      => "40gb"
    },
  };

  return ($config);
}

sub performCutrun {
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

1;
