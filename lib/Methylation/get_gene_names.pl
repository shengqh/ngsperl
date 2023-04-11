open IN,"$ARGV[0]";
open OUT,">$ARGV[1]";

while(<IN>){
  chomp;
  next if(/^Chr/);
  @a = split/\t/,$_;
  if($a[6] =~/\,/){
    @b = split/\,/,$a[6];
    foreach $key(@b){
      $genes{$key}++;
    }
  }
  else{
    $genes{$a[6]}++;
  }
}

foreach $key(sort keys%genes){
  print OUT "$key\n";
}

