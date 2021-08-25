package Alignment::AlignmentUtils;

use strict;
use warnings;
use CQS::ConfigUtils;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw( getSortByCoordinate getChromosomeFilterCommand getAddRgCommand hasMarkDuplicate)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub getSortByCoordinate {
  my ( $config, $section ) = @_;

  my $result;
  if ( defined $config->{$section}{output_sort_by_coordinate} ) {
    $result = get_option( $config, $section, "output_sort_by_coordinate" );
  }
  else {
    $result = get_option( $config, $section, "sort_by_coordinate", 1 );
  }

  return $result;
}

sub getChromosomeFilterCommand {
  my ( $bam_file, $chromosome_grep_pattern ) = @_;

  my $result     = "";
  my $final_file = $bam_file;
  if ( $chromosome_grep_pattern ne "" ) {
    my $tmp_file = $bam_file . ".filtered.bam";
    $result = "
    echo filtering bam by chromosome pattern $chromosome_grep_pattern
    samtools idxstats $bam_file | cut -f 1 | grep $chromosome_grep_pattern | xargs samtools view -b $bam_file > $tmp_file
    samtools flagstat $bam_file > ${bam_file}.raw.stat
    rm $bam_file
    rm ${bam_file}.bai
    mv $tmp_file $bam_file
    samtools index $bam_file
";
  }
  return $result;
}

sub getAddRgCommand {
  my ( $picard, $add_RG_to_read, $bam_file, $sample_name ) = @_;
  my $result = "";
  if ($add_RG_to_read) {
    my $tmp_file = $sample_name . ".rg.bam";
    $result = "
    echo add_RG_to_read_command by picard
    java -jar $picard AddOrReplaceReadGroups I=$bam_file O=$tmp_file ID=$sample_name LB=$sample_name SM=$sample_name PL=illumina PU=$sample_name CREATE_INDEX=False
    rm $bam_file
    rm ${bam_file}.bai
    mv $tmp_file $bam_file
    samtools index $bam_file
";
  }
  return $result;
}

sub hasMarkDuplicate {
  my ($section, $defaultValue) = @_;
  if(not defined $defaultValue){
    $defaultValue = 0;
  }
  return get_option_value( $section->{"mark_duplicates"}, $defaultValue );
}

1;
