package Gcloud::GcloudUtils;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw( replaceSampleNameBWACommand )] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub replaceSampleNameBWACommand {
  my ($bwa, $sample_name) = @_;
  #bwa mem -K 100000000 -v 3 -t 16 -R '@RG\tID:TP0097_Norecur_Base\tPU:illumina\tLB:TP0097_Norecur_Base\tSM:TP0097_Norecur_Base\tPL:illumina' -Y $bash_ref_fasta
  $bwa =~ s/ID:.+?\\/ID:$sample_name\\/;
  $bwa =~ s/LB:.+?\\/LB:$sample_name\\/;
  $bwa =~ s/SM:.+?\\/SM:$sample_name\\/;
  return ($bwa);
}


1;
