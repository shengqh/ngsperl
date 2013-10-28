#!/usr/bin/perl
package CQS::MuTect;

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
use CQS::Fasta;
use CQS::AbstractSomaticMutation;

our @ISA = qw(CQS::AbstractSomaticMutation);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name} = "MuTect";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $muTect_jar = get_param_file( $config->{$section}{muTect_jar},  "muTect_jar",  1 );
  my $faFile     = get_param_file( $config->{$section}{fasta_file},  "fasta_file",  1 );
  my $cosmicfile = get_param_file( $config->{$section}{cosmic_file}, "cosmic_file", 1 );
  my $dbsnpfile  = get_param_file( $config->{$section}{dbsnp_file},  "dbsnp_file",  1 );

  my $bychromosome = $config->{$section}{bychromosome};
  if ( !defined $bychromosome ) {
    $bychromosome = 0;
  }

  my @chromosomes = ();
  if ($bychromosome) {
    @chromosomes = get_sequence_names($faFile);
  }
  else {
    @chromosomes = ("");
  }

  my $java_option = $config->{$section}{java_option};
  if ( !defined $java_option ) {
    $java_option = "";
  }

  my %group_sample_map = %{$self->get_group_sample_map ($config, $section)};

  my $shfile = $pbsDir . "/${task_name}_mt.submit";
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";

  for my $groupName ( sort keys %group_sample_map ) {
    my @sampleFiles = @{ $group_sample_map{$groupName} };
    my $sampleCount = scalar(@sampleFiles);
    my $curDir      = create_directory_or_die( $resultDir . "/$groupName" );

    if ( $sampleCount != 2 ) {
      die "SampleFile should be normal,tumor paired.";
    }

    my $normal = $sampleFiles[0][1];
    my $tumor  = $sampleFiles[1][1];
    
    for my $chr (@chromosomes) {
      my $chrstr;
      my $chroption;
      if ( $chr eq "" ) {
        $chrstr    = "";
        $chroption = "";
      }
      else {
        $chrstr    = "_${chr}";
        $chroption = "-L $chr";
      }

      my $pbsName = "${groupName}${chrstr}_mt.pbs";
      my $pbsFile = "${pbsDir}/$pbsName";
      print SH "\$MYCMD ./$pbsName \n";

      my $out     = "${groupName}${chrstr}.somatic.out";
      my $vcf     = "${groupName}${chrstr}.somatic.vcf";
      my $passvcf = "${groupName}${chrstr}.somatic.pass.vcf";
      my $log     = "${logDir}/${groupName}${chrstr}_mt.log";

      open( OUT, ">$pbsFile" ) or die $!;
      print OUT "$pbsDesc
#PBS -o $log
#PBS -j oe

$path_file

echo muTect=`date` 

cd $curDir

if [ ! -s ${normal}.bai ]; then
  samtools index ${normal}
fi

if [ ! -s ${tumor}.bai ]; then
  samtools index ${tumor}
fi

if [ ! -s $vcf ]; then
  java $java_option -jar $muTect_jar $option $chroption --analysis_type MuTect --reference_sequence $faFile --cosmic $cosmicfile --dbsnp $dbsnpfile --input_file:normal $normal --input_file:tumor $tumor -o $out --vcf $vcf
fi 

if [[ -s $vcf && ! -s $passvcf ]]; then
  awk '{if (\$1 ~ /^##INFO/) {exit;} else {print;}}' $vcf > $passvcf
  echo \"##INFO=<ID=LOD,Number=1,Type=Float,Description=\\\"Log of (likelihood tumor event is real / likelihood event is sequencing error )\\\">\" >> $passvcf
  awk 'BEGIN{info=0}{if (\$1 ~ /^##INFO/) {info=1;} if(info) {print;}}' $vcf | grep \"^#\" >> $passvcf
  grep -v \"^##\" $out | awk '{print \";LOD=\"\$17}' > ${out}.lod
  grep -v \"^##\" $vcf | cut -f1-8 > ${vcf}.first
  paste -d \"\" ${vcf}.first ${out}.lod > ${vcf}.first_lod
  grep -v \"^##\" $vcf | cut -f9- > ${vcf}.second
  paste ${vcf}.first_lod ${vcf}.second | grep -v \"^#CHROM\" | grep -v REJECT >> $passvcf
  rm ${out}.lod ${vcf}.first ${vcf}.first_lod ${vcf}.second
fi

echo finished=`date` \n";

      close OUT;

      print "$pbsFile created \n";
    }
  }

  close(SH);

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all MuTect tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $groupName ( keys %{$groups} ) {
    my @resultFiles = ();
    my $curDir      = $resultDir . "/$groupName";
    push( @resultFiles, "$curDir/${groupName}.somatic.pass.vcf" );
    $result->{$groupName} = filter_array( \@resultFiles, $pattern );
  }
  return $result;
}

1;
