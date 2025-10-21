open IN,"$ARGV[0]";
open OUT,">$ARGV[1]";

$header_found = 0;
$gene_col_index = -1;

while(<IN>){
  chomp;
  
  # Find header line and get Gene.refGene column index
  if(!$header_found && /^Chr|^chr\s+/) {
    @header = split/\t/,$_;
    for my $i (0..$#header) {
      if($header[$i] eq "Gene.refGene") {
        $gene_col_index = $i;
        last;
      }
    }
    $header_found = 1;
    next;
  }
  
  next if(/^Chr/);
  next if(/^chr\s+/);
  next if($gene_col_index == -1);
  
  @a = split/\t/,$_;
  if($a[$gene_col_index] =~/[,;]+/){
    @b = split/[,;]+/,$a[$gene_col_index];
    foreach $key(@b){
      $genes{$key}++;
    }
  }
  else{
    $genes{$a[$gene_col_index]}++;
  }
}

foreach $key(sort keys%genes){
  print OUT "$key\n";
}

