package Alignment::AlignmentUtils;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw( getChromosomeFilterCommand)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

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
