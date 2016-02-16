use strict;

use Bio::SeqIO;

my $file = "H:/shengquanhu/projects/database/ucsc/GtRNAdb2.20160216.fa";
my $seqio = Bio::SeqIO->new( -file => $file, -format => 'fasta' );

my $seqnames = {};
while ( my $seq = $seqio->next_seq ) {
  print $seq->id, "\n";
  exit;
}

1;
