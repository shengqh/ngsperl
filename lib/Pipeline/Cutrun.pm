#!/usr/bin/perl
package Pipeline::Cutrun;

#this is a wrapper for cutruntools2, except for peak calling, which supports using igg or input as control.
#https://github.com/fl-yu/CUT-RUNTools-2.0/blob/cf72ca5d0801ccab46ca7cfdb6810628797f4c9b/src/bulk/bulk-pipeline.sh

use strict;
use warnings;
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
  initDefaultValue( $def, "mark_duplicates", 0 );

  initDefaultValue( $def, "min_read_length", 25);

  initDefaultValue( $def, "bowtie2_option", "--dovetail --phred33" );

  initDefaultValue( $def, "perform_macs2_broad", 1 );
  initDefaultValue( $def, "perform_macs2_narrow", 1 );
  initDefaultValue( $def, "perform_seacr", 1 );

  initDefaultValue( $def, "macs2_broad_option", "-f BAMPE --broad --broad-cutoff 0.1 -B --SPMR --keep-dup all");
  initDefaultValue( $def, "macs2_narrow_option", "-f BAMPE -q 0.01 -B --SPMR --keep-dup all");

  initDefaultValue( $def, "annotate_nearest_gene", 1 );

  if(defined $def->{treatments}){
    $def->{treatments_auto} = 0;
  }else{
    initDefaultValue( $def, "treatments_auto",     1 );
  }
  
  return $def;
}

sub getConfig {
  my ($def) = @_;
  $def->{VERSION} = $VERSION;

  my $target_dir = $def->{target_dir};
  create_directory_or_die($target_dir);

  $def = initializeDefaultOptions($def);

  my ( $config, $individual_ref, $summary_ref, $source_ref, $preprocessing_dir, $untrimed_ref, $cluster ) = getPreprocessionConfig($def);
  #merge summary and individual 
  push @$individual_ref, @$summary_ref;
  $summary_ref = $individual_ref;

  #https://github.com/fl-yu/CUT-RUNTools-2.0/blob/master/src/bulk/bulk-pipeline.sh

  my $genome_sequence = getValue($def, "fasta_file");

  my $cutruntools2_path = getValue($def, "cutruntools2_path");
  my $extratoolsbin = $cutruntools2_path . "/install";
  my $adapterpath = $cutruntools2_path . "/adapters";

  my $trimmomaticjarfile = "$extratoolsbin/trimmomatic-0.36.jar";
  my $kseq_test = "$extratoolsbin/kseq_test";
  my $len = getValue($def, "fastq_sequence_length");
  my $adapter_type = getValue($def, "adaptor_type");

  my $adapter_file = $adapter_type eq "Nextera" ? "NexteraPE-PE.fa": "Truseq3.PE.fa";

  my $trimmomatic_task = "trimmomatic";
  $config->{ $trimmomatic_task } = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "${target_dir}/" . getNextFolderIndex($def) . "$trimmomatic_task",
    interpretor => "",
    program => "",
    check_program => 0,
    option => "
rm -f *.failed *.succeed

echo trimmomatic=`date`
java -jar $trimmomaticjarfile PE -threads 8 -phred33 \\
  __FILE__ \\
  __NAME__.clipped.1.fastq.gz __NAME__.1.unpaired.fastq.gz \\
  __NAME__.clipped.2.fastq.gz __NAME__.2.unpaired.fastq.gz \\
  ILLUMINACLIP:$adapterpath/$adapter_file:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25 > __NAME__.trimmomatic.log 2>&1    

status=\$?
if [[ \$status -eq 0 ]]; then
  rm -f __NAME__.1.unpaired.fastq.gz __NAME__.2.unpaired.fastq.gz
  touch __NAME__.trimmomatic.succeed
else
  rm -f __NAME__.1.unpaired.fastq.gz __NAME__.2.unpaired.fastq.gz __NAME__.clipped.1.fastq.gz __NAME__.clipped.2.fastq.gz
  echo \$status > __NAME__.trimmomatic.failed
  exit \$status
fi

# echo kseq_test=`date`
# It doesn't make sense to use kseq_test to trim the reads after trimmomatic
# $kseq_test __NAME__.1.paired.fastq.gz $len __NAME__.clipped.1.fastq.gz
# $kseq_test __NAME__.2.paired.fastq.gz $len __NAME__.clipped.2.fastq.gz
# rm -f __NAME__.1.paired.fastq.gz __NAME__.2.paired.fastq.gz

echo v0.36 > __NAME__.trimmomatic.version
",
    source_ref            => $untrimed_ref,
    source_join_delimiter => " ",
    output_to_same_folder => 0,
    output_file_ext => ".clipped.1.fastq.gz,.clipped.2.fastq.gz,.trimmomatic.version",
    sh_direct             => 0,
    no_output             => 1,
    docker_prefix => "cutruntools2_",
    pbs                   => {
      "nodes"    => "1:ppn=8",
      "walltime" => getValue($def, "trimmomatic_walltime", "24"),
      "mem"      => getValue($def, "trimmomatic_mem", "40gb")
    },
  };
  push @$summary_ref, ( $trimmomatic_task );

  if ( $def->{perform_fastqc} ) {
    addFastQC( $config, $def, $summary_ref, $summary_ref, "${trimmomatic_task}_fastqc", [ $trimmomatic_task, ".fastq.gz" ], $preprocessing_dir );
  }

  my $bowtie2_option = getValue( $def, "bowtie2_option" );
  my $bowtie2_index = getValue( $def, "bowtie2_index" );
  my $awk1 = "$extratoolsbin/filter_below.awk";

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

echo filter_120bp=`date`
samtools view -h __NAME__.unsorted.bam | LC_ALL=C awk -f $awk1 | samtools view -Sb - > __NAME__.unsorted.120bp.bam

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.120bp.succeed
  rm __NAME__.unsorted.bam
else
  echo \$status > __NAME__.120bp.failed
  rm -f __NAME__.unsorted.120bp.bam
  exit \$status
fi

echo sort=`date`
samtools sort -@ 8 -o __NAME__.120bp.bam -T __NAME__ __NAME__.unsorted.120bp.bam

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.sort.succeed
  rm __NAME__.unsorted.120bp.bam
else
  echo \$status > __NAME__.sort.failed
  rm -f __NAME__.120bp.bam
  exit \$status
fi

echo index=`date`
samtools index __NAME__.120bp.bam

echo idxstats=`date`
samtools idxstats __NAME__.120bp.bam > __NAME__.120bp.bam.chromosome.count

echo flagstat=`date`
samtools flagstat __NAME__.120bp.bam > __NAME__.120bp.bam.stat 

bowtie2 --version | grep bowtie2 | grep version | cut -d ' ' -f3 | awk '{print \"bowtie2,v\"\$1}' > __NAME__.bam.version

",
    source_ref            => [$trimmomatic_task, ".fastq.gz\$"],
    source_join_delimiter => " -2 ",
    output_to_same_folder => 0,
    no_output => 1,
    output_file_ext => ".120bp.bam,.120bp.bam.chromosome.count,.120bp.bam.stat,.log,.bam.version",
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

  add_alignment_summary($config, $def, $summary_ref, $target_dir, "${bowtie2_task}_summary", "../Alignment/AlignmentUtils.r;../Alignment/Bowtie2Summary.r", ".reads.csv;.reads.png;.chromosome.csv;.chromosome.png", [ "bowtie2", ".log" ], ["bowtie2", ".chromosome.count"] );

  my $bam_files = get_result_file($config, $bowtie2_task, ".bam\$");
  my $treatment_samples = do_get_group_samplefile_map(getValue($def, "treatments"), $bam_files);
  my $control_samples = do_get_group_samplefile_map(getValue($def, "controls"), $bam_files);
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

  my $seacr_task = "seacr";
  $config->{ $seacr_task } = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "${target_dir}/" . getNextFolderIndex($def) . "$seacr_task",
    interpretor => "",
    program => "",
    check_program => 0,
    option => "
rm -f *.failed *.succeed

echo change.bdg=`date` 
python $extratoolsbin/change.bdg.py __FILE__ > __NAME___treat_integer.bdg

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.change.bdg.succeed
else
  echo \$status > __NAME__.change.bdg.failed
  exit \$status
fi

executable_path=\$(command -v -- \"Rscript\")
Rscriptbin=\$(dirname -- \"\$executable_path\")

echo SEACR=`date` 
bash $extratoolsbin/SEACR_1.1.sh __NAME___treat_integer.bdg 0.01 non stringent __NAME___treat \$Rscriptbin

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.SEACR.succeed
else
  echo \$status > __NAME__.SEACR.failed
  exit \$status
fi

echo sort-bed=`date` 
sort-bed __NAME___treat.stringent.bed > __NAME___treat.stringent.sort.bed

echo get_summits_seacr=`date` 
python $extratoolsbin/get_summits_seacr.py __NAME___treat.stringent.bed | sort-bed - > __NAME___treat.stringent.sort.summits.bed

#__OUTPUT__
",
    source_ref =>  [$macs2_narrow_task, "_treat_pileup.bdg\$"],
    output_to_same_folder => 0,
    output_file_ext => "_treat.stringent.sort.bed,_treat.stringent.sort.summits.bed",
    sh_direct             => 0,
    docker_prefix => "cutruntools2_",
    pbs                   => {
      "nodes"    => "1:ppn=1",
      "walltime" => getValue($def, "seacr_walltime", "24"),
      "mem"      => getValue($def, "seacr_mem", "40gb")
    },
  };
  push @$summary_ref, ( $seacr_task );

  if($def->{annotate_nearest_gene}){
    my $gene_bed = getValue($def, "gene_bed");
    add_nearest_gene($config, $def, $summary_ref, $target_dir, $seacr_task, "_treat.stringent.sort.bed", $gene_bed);
  }

  my $peak_file_ref = [ $seacr_task, "_treat.stringent.sort.bed" ];

  my $summit_ref = [ $seacr_task, "_treat.stringent.sort.summits.bed" ];
  my $suffix="_treat.stringent.sort.bed";
  my $summit_suffix="_treat.stringent.sort.summits.bed";
  my $summit_padded_suffix="_treat.stringent.sort.summits_padded.fa";

  my $blacklist;
  if($macs2_genome eq "hs"){
    $blacklist = "$cutruntools2_path/blacklist/hg38.blacklist.bed";
  }elsif($macs2_genome eq "mm"){
    $blacklist = "$cutruntools2_path/blacklist/mm10.blacklist.bed";
  }else{
    die "Unknown macs2_genome $macs2_genome";
  }

  my $num_bp_from_summit = 100;
  my $num_peaks = 1000;
  my $total_peaks = 2000;
  my $motif_scanning_pval = 0.0005;
  my $num_motifs = 10;

  my $tmp_dir="tmp";
  my $blacklist_filtered_dir="$tmp_dir/blacklist_filtered";
  my $msummit_dir="$tmp_dir/summits";
  my $mpadded_dir="$tmp_dir/padded";
  my $mpaddedfa_dir="$tmp_dir/padded.fa";

  my $summitfa="__NAME__$summit_padded_suffix";

  my $peak="__NAME__.peak.bed";
  my $summit="__NAME__.summit.bed";

  my $motif_task = "seacr_motif";
  $config->{ $motif_task } = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "${target_dir}/" . getNextFolderIndex($def) . "$motif_task",
    interpretor => "",
    program => "",
    check_program => 0,
    option => "
rm -f *.failed *.succeed

for d in $tmp_dir $blacklist_filtered_dir $msummit_dir $mpadded_dir $mpaddedfa_dir; do
  if [ ! -d \$d ]; then
    mkdir -p \$d
  fi
done

echo motif=`date` 
cat __FILE__ | grep -v -e \"chrM\" | sort-bed - | bedops -n 1 - $blacklist > $blacklist_filtered_dir/$peak
cat __FILE2__ | grep -v -e \"chrM\" | sort-bed - | bedops -n 1 - $blacklist > $blacklist_filtered_dir/$summit

echo \"[info] Get randomized [$num_peaks] peaks from the top [$total_peaks] peaks...\"
echo \"[info] Filtering the blacklist regions for the selected peak files\"

cat $blacklist_filtered_dir/$peak | sort -t\$'\\t' -g -k5 -r | head -n $total_peaks | shuf | head -n $num_peaks | sort-bed - > $tmp_dir/$peak

bedops -e 1 $blacklist_filtered_dir/$summit $tmp_dir/$peak > $msummit_dir/$summit
bedops --range $num_bp_from_summit -u $msummit_dir/$summit > $mpadded_dir/$summit

echo \"[info] Getting Fasta sequences\"
bedtools getfasta -fi $genome_sequence -bed $mpadded_dir/$summit -fo $mpaddedfa_dir/$summitfa

echo \"[info] Start MEME analysis for de novo motif finding ...\"
echo \"[info] Up to $num_motifs will be output ...\"
meme-chip -oc . -dreme-m $num_motifs -meme-nmotifs $num_motifs $mpaddedfa_dir/$summitfa

echo \"[info] Saving the De Novo motifs ...\"
python $extratoolsbin/read.meme.py .

rm -rf $tmp_dir

echo \"[info] De Novo motifs can be found: `pwd` ...\"

",
    source_ref =>  $peak_file_ref,
    parameterSampleFile2_ref => $summit_ref,
    output_to_same_folder => 0,
    samplename_in_result => 0,
    output_file_ext => "motifs,meme-chip.html",
    sh_direct => 0,
    no_output => 1,
    docker_prefix => "cutruntools2_",
    pbs                   => {
      "nodes"    => "1:ppn=1",
      "walltime" => getValue($def, "motif_walltime", "24"),
      "mem"      => getValue($def, "motif_mem", "40gb")
    },
  };
  push @$summary_ref, ( $motif_task );


  my $fa_dir="$tmp_dir/blacklist_filtered.fa";
  my $fimo_task = "seacr_motif_fimo";
  $config->{ $fimo_task } = {
    class                 => "CQS::ProgramWrapperOneToOne",
    perform               => 1,
    target_dir            => "${target_dir}/" . getNextFolderIndex($def) . "$fimo_task",
    interpretor => "",
    program => "",
    check_program => 0,
    option => "
rm -f *.failed *.succeed

p=$motif_scanning_pval
echo \"[info] The signficance cutoff of Fimo scaning is \${p}...\"

mbase='__NAME__'
motif_dir='__FILE2__' # a directory containing a list of *.meme files
outbam='__FILE3__'

echo \"[info] Motif files can be found: \$motif_dir\"
if [ ! -d $blacklist_filtered_dir ]; then
  mkdir -p $blacklist_filtered_dir
fi

if [ ! -d $fa_dir ]; then
  mkdir -p $fa_dir
fi

echo \"[info] Filtering the blacklist regions for the selected peak files\"
cat '__FILE__' | grep -v -e \"chrM\" | sort-bed - | bedops -n 1 - '$blacklist' > '$blacklist_filtered_dir/$peak'

echo \"[info] Getting Fasta sequences\"
bedtools getfasta -fi $genome_sequence -bed $blacklist_filtered_dir/$peak -fo $fa_dir/\${mbase}.fa
python $extratoolsbin/fix_sequence.py $fa_dir/\${mbase}.fa

echo \"[info] Scaning the De Novo motifs for each peak\"
for m in `ls -1 \$motif_dir`; do
  motif=`basename \$m .meme`
  fimo_d=fimo2.\$motif
  if [ ! -d \$fimo_d ]; then
    mkdir \$fimo_d
  fi
  fimo --thresh \$p --parse-genomic-coord -oc \$fimo_d \$motif_dir/\${motif}.meme $fa_dir/\${mbase}.fa
  gff2bed < \$fimo_d/fimo.gff | awk 'BEGIN {IFS=\"\\t\"; OFS=\"\\t\";} {print \$1,\$2,\$3,\$4,\$5,\$6}' > \$fimo_d/fimo.bed

  if [[ -s \$fimo_d/fimo.bed ]]; then
    tmp=`echo \$motif | cut -d '.' -f2 | wc -c`
    mlen=\$(( tmp - 1 ))
    make_cut_matrix -v -b '(25-150 1)' -d -o 0 -r 100 -p 1 -f 3 -F 4 -F 8 -q 0 \$outbam \$fimo_d/fimo.bed > \$fimo_d/fimo.cuts.freq.txt
    Rscript $extratoolsbin/run_centipede_parker.R \$fimo_d/fimo.cuts.freq.txt \$fimo_d/fimo.bed \$fimo_d/fimo.pdf \$mlen
  else
    touch \$fimo_d/fimo.bed.is_empty
    echo \$fimo_d/fimo.bed is empty, ignored
  fi
done

echo \"[info] Output can be found: `pwd`\"
",
    source_ref =>  $peak_file_ref,
    parameterSampleFile2_ref => [ $motif_task, "motifs"],
    parameterSampleFile3 => $treatment_samples,
    output_to_same_folder => 0,
    samplename_in_result => 0,
    output_file_ext => "fimo.html",
    sh_direct => 0,
    no_output => 1,
    docker_prefix => "cutruntools2_",
    pbs                   => {
      "nodes"    => "1:ppn=1",
      "walltime" => getValue($def, "motif_walltime", "24"),
      "mem"      => getValue($def, "motif_mem", "40gb")
    },
  };
  push @$summary_ref, ( $fimo_task );

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
