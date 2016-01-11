#!/usr/bin/perl
package GATK::MuTect;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::Fasta;
use CQS::GroupTask;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = "GATK::MuTect";
  $self->{_suffix} = "_mt";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbsDesc, $target_dir, $logDir, $pbsDir, $resultDir, $option, $sh_direct, $cluster, $thread ) = get_parameter( $config, $section );

  my $muTect_jar = get_param_file( $config->{$section}{muTect_jar},  "muTect_jar",  1 );
  my $faFile     = get_param_file( $config->{$section}{fasta_file},  "fasta_file",  1 );
  my $dbsnpfile  = get_param_file( $config->{$section}{dbsnp_file},  "dbsnp_file",  1 );
  
  #mouse genome has no corresponding cosmic database
    my $cosmic_param = "";
  my $cosmicfile = get_param_file( $config->{$section}{cosmic_file}, "cosmic_file", 0 );
  if(defined $cosmicfile){
    $cosmic_param = "--cosmic $cosmicfile";
  }
  
  my $java = get_java($config, $section);

  my $java_option = get_option( $config, $section, "java_option", "" );

  my %group_sample_map = %{ get_group_sample_map( $config, $section ) };

  my $shfile = $self->taskfile( $pbsDir, $task_name );
  open( SH, ">$shfile" ) or die "Cannot create $shfile";
  print SH get_run_command($sh_direct) . "\n";
  print SH "cd $pbsDir\n";

  for my $groupName ( sort keys %group_sample_map ) {
    my @sampleFiles = @{ $group_sample_map{$groupName} };
    my $sampleCount = scalar(@sampleFiles);
    my $curDir      = create_directory_or_die( $resultDir . "/$groupName" );

    if ( $sampleCount != 2 ) {
      die "SampleFile should be normal,tumor paired.";
    }

    my $normal = $sampleFiles[0][1];
    my $tumor  = $sampleFiles[1][1];

    my $pbsFile = $self->pbsfile( $pbsDir, $groupName );
    my $pbsName = basename($pbsFile);
    my $log     = $self->logfile( $logDir, $groupName );

    print SH "\$MYCMD ./$pbsName \n";

    my $out     = "${groupName}.somatic.out";
    my $vcf     = "${groupName}.somatic.vcf";
    my $passvcf = "${groupName}.somatic.pass.vcf";

    my $log_desc = $cluster->get_log_desc($log);

    open( OUT, ">$pbsFile" ) or die $!;
    print OUT "$pbsDesc
$log_desc

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
  $java $java_option -jar $muTect_jar $option --analysis_type MuTect --reference_sequence $faFile $cosmic_param --dbsnp $dbsnpfile --input_file:normal $normal --input_file:tumor $tumor -o $out --vcf $vcf --enable_extended_output
fi 

if [[ -s $vcf && ! -s $passvcf ]]; then
  if [ -s $out ]; then
    awk '{if (\$1 ~ /^##INFO/) {exit;} else {print;}}' $vcf > $passvcf
    echo \"##INFO=<ID=LOD,Number=1,Type=Float,Description=\\\"Log of (likelihood tumor event is real / likelihood event is sequencing error )\\\">\" >> $passvcf
    awk 'BEGIN{info=0}{if (\$1 ~ /^##INFO/) {info=1;} if(info) {print;}}' $vcf | grep \"^#\" >> $passvcf
    grep -v \"^##\" $out | awk -vCOLM=t_lod_fstar 'NR==1 { for (i = 1; i <= NF; i++) { if (\$i == COLM) { cidx = i; } } } NR > 1 {print \";LOD=\"\$cidx}' > ${out}.lod
    grep -v \"^#\" $vcf | cut -f1-8 > ${vcf}.first
    paste -d \"\" ${vcf}.first ${out}.lod > ${vcf}.first_lod
    grep -v \"^#\" $vcf | cut -f9- > ${vcf}.second
    paste ${vcf}.first_lod ${vcf}.second | grep -v \"^#CHROM\" | grep -v REJECT >> $passvcf
    rm ${out}.lod ${vcf}.first ${vcf}.first_lod ${vcf}.second
    grep -v \"^#\" $passvcf | cut -f1 | uniq -c | awk '{print \$2\"\\t\"\$1}' > ${passvcf}.chromosome
    
  else
    grep -v REJECT $vcf > $passvcf
  fi
fi

echo finished=`date` \n";

    close OUT;

    print "$pbsFile created \n";
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
