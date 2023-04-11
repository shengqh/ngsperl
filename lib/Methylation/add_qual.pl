open IN,"$ARGV[0]";
open OUT,">$ARGV[1]";

while(<IN>){
	chomp;
	if(/^\@/){
		print OUT "$_\n";
	}
	else{
		@a=split/\t/,$_;
		$rl=length($a[9]);
		$qual="F" x $rl;
		print OUT "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\t$qual\t$a[11]\t$a[12]\n";
	}
}

