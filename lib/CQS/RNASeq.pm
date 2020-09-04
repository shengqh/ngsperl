#!/usr/bin/perl
package CQS::RNASeq;

use strict;
use warnings;
use File::Copy;
use File::Basename;
use List::Compare;
use CQS::PBS;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::NGSCommon;
use CQS::ClassFactory;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
  'all' => [
    qw(tophat2 call_tophat2 tophat2_by_pbs get_tophat2_result cufflinks_by_pbs cuffmerge_by_pbs cuffdiff_by_pbs read_cufflinks_fpkm read_cuffdiff_significant_genes copy_and_rename_cuffdiff_file compare_cuffdiff novoalign shrimp2)
  ]
);

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

use Cwd;

sub get_sorted_bam_prefix {
  my ($oldbam) = @_;
  my ( $filename, $dirs, $suffix ) = fileparse( $oldbam, qr/\.[^.]*/ );
  return ( $filename . "_sorted" );
}

sub tophat2_by_pbs {
  my ( $config, $section ) = @_;
  my $obj = instantiate("Tophat2");
  $obj->perform( $config, $section );
}

sub tophat2 {
  tophat2_by_pbs(@_);
}

sub call_tophat2 {
  tophat2_by_pbs(@_);
}

#get expected tophat2 result based on tophat2 definition
sub get_tophat2_map {
  my ( $config, $section ) = @_;

  my ( $result, $issource ) = get_raw_files2( $config, $section );
  if ($issource) {
    retun $result;
  }

  my $tophatsection = $config->{$section}{source_ref};
  my $tophat_dir = $config->{$tophatsection}{target_dir} or die "${tophatsection}::target_dir not defined.";
  my ( $log_dir, $pbs_dir, $result_dir ) = init_dir( $tophat_dir, 0 );
  my %fqFiles = %{$result};

  my $tpresult = {};
  for my $sample_name ( keys %fqFiles ) {
    my $bam = "${result_dir}/${sample_name}/accepted_hits.bam";
    $tpresult->{$sample_name} = $bam;
  }
  return $tpresult;
}

sub cufflinks_by_pbs {
  my ( $config, $section ) = @_;
  my $obj = instantiate("Cufflinks");
  $obj->perform( $config, $section );
}

sub cuffmerge_by_pbs {
  my ( $config, $section ) = @_;
  my $obj = instantiate("Cuffmerge");
  $obj->perform( $config, $section );
}

sub cuffdiff_by_pbs {
  my ( $config, $section ) = @_;
  my $obj = instantiate("Cuffdiff");
  $obj->perform( $config, $section );
}

sub check_is_single() {
  my ( $sample_nameCount, @sample_files ) = @_;
  my $sampleFileCount = scalar(@sample_files);
  my $isSingle        = 1;
  if ( $sample_nameCount == $sampleFileCount ) {
    $isSingle = 1;
  }
  elsif ( $sample_nameCount * 2 == $sampleFileCount ) {
    $isSingle = 0;
  }
  else {
    die "Count of sample_name should be equals to count/half count of sample_files";
  }

  return ($isSingle);
}

sub read_cufflinks_fpkm {
  my $genefile = shift;
  open GF, "<$genefile" or die "Cannot open file $genefile";
  my $header = <GF>;
  my @headers = split( /\t/, $header );
  my ($geneindex) = grep { $headers[$_] eq "gene_id" } 0 .. $#headers;
  my ($fpkmindex) = grep { $headers[$_] eq "FPKM" } 0 .. $#headers;

  my $result = {};
  while ( my $line = <GF> ) {
    my @values = split( /\t/, $line );
    $result->{ $values[$geneindex] } = $values[$fpkmindex];
  }
  close(GF);

  return $result;
}

sub read_cuffdiff_significant_genes {
  my $file = shift;

  my $result = {};
  open IN, "<$file" or die "Cannot open file $file";
  my $header = <IN>;
  chomp($header);
  while ( my $line = <IN> ) {
    chomp($line);
    my @part = split( /\t/, $line );
    if ( $part[13] eq "yes" ) {
      $result->{ $part[2] } = $line;
    }
  }
  close IN;

  return ( $result, $header );
}

#use CQS::RNASeq;
#my $root = "/home/xxx/cuffdiff/result";
#my $targetdir = create_directory_or_die( $root . "/comparison" );
#my $config    = {
#    "RenameDiff" => {
#        target_dir => $targetdir,
#        root_dir   => $root
#    }
#}
#copy_and_rename_cuffdiff_file($config, "RenameDiff");
sub copy_and_rename_cuffdiff_file {
  my ( $config, $section ) = @_;
  my $dir = $config->{$section}{"root_dir"} or die "define ${section}::root_dir first";
  if ( !-d $dir ) {
    die "directory $dir is not exists";
  }

  my $targetdir = $config->{$section}{"target_dir"} or die "define ${section}::target_dir first";
  if ( !-d $targetdir ) {
    create_directory_or_die($targetdir);
  }

  my $gene_only = $config->{$section}{"gene_only"};
  if ( !defined($gene_only) ) {
    $gene_only = 0;
  }

  my @subdirs = list_directories($dir);
  if ( 0 == scalar(@subdirs) ) {
    die "$dir has no sub CuffDiff directories";
  }

  my @filenames;
  if ($gene_only) {
    @filenames = ("gene_exp.diff");
  }
  else {
    @filenames = ( "gene_exp.diff", "splicing.diff" );
  }

  for my $subdir (@subdirs) {
    foreach my $filename (@filenames) {
      my $file = "${dir}/${subdir}/${filename}";
      if ( -s $file ) {
        print "Dealing " . $file . "\n";
        open IN, "<$file" or die "Cannot open file $file";
        my $line = <IN>;
        $line = <IN>;
        close(IN);

        if ( defined($line) ) {
          my @parts      = split( /\t/, $line );
          my $partcount  = scalar(@parts);
          my $targetname = "${targetdir}/${subdir}.${filename}";

          copy( $file, $targetname ) or die "copy failed : $!";

          my $target_sign_name = $targetname . ".sig";
          my $cmd;
          if ($gene_only) {
            $cmd = "cat $targetname | awk '(\$3 != \"-\") && (\$$partcount==\"yes\" || \$$partcount==\"significant\")' > $target_sign_name";
          }
          else {
            $cmd = "cat $targetname | awk '\$$partcount==\"yes\" || \$$partcount==\"significant\"' > $target_sign_name";
          }

          print $cmd . "\n";
          `$cmd`;
        }
      }
    }
  }
}

sub compare_cuffdiff {
  my ( $config, $section ) = @_;

  my $info = $config->{$section};

  my ( $file1, $file2 ) = @{ $info->{"files"} };

  my ( $data1, $header )  = read_cuffdiff_significant_genes($file1);
  my ( $data2, $header2 ) = read_cuffdiff_significant_genes($file2);

  my @genes1 = keys %{$data1};
  my @genes2 = keys %{$data2};

  my $lc = List::Compare->new( \@genes1, \@genes2 );

  my @resultgenes = ();
  if ( $info->{operation} eq "minus" ) {
    @resultgenes = $lc->get_Lonly();
  }
  elsif ( $info->{operation} eq "intersect" ) {
    @resultgenes = $lc->get_intersection();
  }
  else {
    die "Only minus or intersect is supported.";
  }

  my $result_fileName = $info->{"target_file"};

  output_compare_cuffdiff_file( $data1, $data2, $result_fileName, $header, @resultgenes );
}

sub novoalign {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option ) = $self->init_parameter( $config, $section );

  my $novoindex = get_param_file( $config->{$section}{novoindex}, "novoindex", 1 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbs_dir . "/${task_name}.sh";
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };

    my $sampleFile = $sample_files[0];

    my $sam_file         = $sample_name . ".sam";
    my $bam_file         = $sample_name . ".bam";
    my $sortedBamPrefix = $sample_name . "_sorted";
    my $sortedbam_file   = $sortedBamPrefix . ".bam";

    my $pbs_name = "${sample_name}_nalign.pbs";
    my $pbs_file = "${pbs_dir}/$pbs_name";

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log    = "${log_dir}/${sample_name}_nalign.log";
    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    open( my $out, ">$pbs_file" ) or die $!;
    print $out "$pbs_desc
#PBS -o $log
#PBS -j oe

$path_file

echo novoalign=`date`

cd $cur_dir

novoalign -d $novoindex -f $sampleFile $option -o SAM > $sam_file

samtools view -b -S $sam_file -o $bam_file 

samtools sort $bam_file $sortedBamPrefix 

samtools index $sortedbam_file 

samtools flagstat $sortedbam_file > ${sortedbam_file}.stat 

echo finished=`date` 
";

    close $out;

    print "$pbs_file created \n";
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all bwa tasks.\n";

  #`qsub $pbs_file`;
}

sub shrimp2 {
  my ( $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section );

  my $genome_index = $config->{$section}{genome_index} or die "define ${section}::genome_index first";
  die "genome index ${genome_index}.genome not exist" if ( !-e "${genome_index}.genome" );
  my $is_mirna   = $config->{$section}{is_mirna}   or die "define ${section}::is_mirna first";
  my $output_bam = $config->{$section}{output_bam} or die "define ${section}::output_bam first";
  my $mirna = "-M mirna" if $is_mirna or "";

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $pbs_dir . "/${task_name}.sh";
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  if ($sh_direct) {
    print $sh "export MYCMD=\"bash\" \n";
  }
  else {
    print $sh "type -P qsub &>/dev/null && export MYCMD=\"qsub\" || export MYCMD=\"bash\" \n";
  }

  for my $sample_name ( sort keys %raw_files ) {
    my $sampleFile = $raw_files{$sample_name};

    my $shrimpFile      = $sample_name . ".shrimp";
    my $sam_file         = $sample_name . ".sam";
    my $bam_file         = $sample_name . ".bam";
    my $sortedBamPrefix = $sample_name . "_sorted";
    my $sortedbam_file   = $sortedBamPrefix . ".bam";

    my $pbs_name = "${sample_name}_shrimp2.pbs";
    my $pbs_file = "${pbs_dir}/$pbs_name";

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log    = "${log_dir}/${sample_name}_shrimp2.log";
    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    open( my $out, ">$pbs_file" ) or die $!;

    print $out "$pbs_desc
#PBS -o $log
#PBS -j oe

$path_file

echo shrimp2=`date`

cd $cur_dir
";

    if ($output_bam) {
      print $out "gmapper -L $genome_index $sampleFile $mirna $option --extra-sam-fields >$sam_file

if [ -s $sam_file ]; then
  samtools view -b -S $sam_file -o $bam_file 
  samtools sort $bam_file $sortedBamPrefix 
  samtools index $sortedbam_file 
  samtools flagstat $sortedbam_file > ${sortedbam_file}.stat 
fi

echo finished=`date` 
";
    }
    else {
      print $out "gmapper -L $genome_index $sampleFile $mirna $option --pretty >$shrimpFile
      
echo finished=`date` 
";
    }

    close $out;

    print "$pbs_file created \n";
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all shrimp2 tasks.\n";

  #`qsub $pbs_file`;
}

1;
