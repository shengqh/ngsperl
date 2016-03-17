use strict;
use warnings;

use Getopt::Long;

my $usage = "

Synopsis:

perl GTEx.pl -i snp_file -d GTEx_directory -o result_file

Options:
  -i|--input {file}         Input snp file with first three columns : chromosome, name and position
  -d|--gtex_dir {directory} Input GTEx directory contains all snp-tissue databases
  -o|--output {file}        Output file
  -h|--help                 This page.
";

Getopt::Long::Configure('bundling');

my $input_file;
my $gtex_dir;
my $output_file;
my $help;

GetOptions(
  'i|input=s'    => \$input_file,
  'd|gtex_dir=s' => \$gtex_dir,
  'o|output=s'   => \$output_file,
  'h|help'       => \$help,
);

if ( defined $help ) {
  print $usage;
  exit(1);
}

if ( !defined $input_file || !-e $input_file ) {
  die "Input file is not valid!";
}

if ( !defined $gtex_dir || !-e $gtex_dir ) {
  die "GTEx directory is not valid!";
}

if ( !defined $output_file ) {
  die "Output file is required!";
}

die "$output_file already exists" if ( -e $output_file );

my %res;
open( my $input, "<$input_file" ) or die "Cannot open $input_file";
while (<$input>) {
  chomp;
  my @parts = split( "\t", $_ );
  if ( scalar(@parts) < 3 ) {
    next;
  }

  my $chrom    = $parts[0];
  my $position = $parts[2];

  if ( $position =~ /^\d+$/ ) {
    my $key = "^" . $chrom . "_" . $position . "_";
    my $cmd = "grep \"$key\" ${gtex_dir}/*.snpgenes |";
    print $cmd, "\n";
    open( my $find, $cmd ) or die "Cannot execute grep command $cmd ";
    my @target = ();
    while(<$find>){
      chomp;
      print $_ . "\n";
      my @parts = split( "\t", $_ );
      push(@target, $parts[0] . ":" . $parts[26]);
    }
    close($find);
    if(scalar(@target) > 0){
      print $key, "\t", join("/", @target), "\n";
    } 
  }
}
close($input);
