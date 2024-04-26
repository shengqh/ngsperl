#!/usr/bin/perl
package Annotation::Annovar;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_ann";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $buildver = $config->{$section}{buildver} or die "buildver is not defined in $section";
  $option = "-buildver $buildver $option";

  my $clean_folder = get_option( $config, $section, "clean_folder", 1 );
  my $perform_splicing = get_option( $config, $section, "perform_splicing", 1 );
  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder", 0 );

  my $annovar_gzipped = get_option( $config, $section, "annovar_gzipped", 0 );

  my $annovarDB = $config->{$section}{annovar_db} or die "annovar_db is not defined in $section";
  my $isvcf = $config->{$section}{isvcf};
  if ( !defined $isvcf ) {
    $isvcf = 0;
  }
  my $isBed = $config->{$section}{isBed};
  if ( !defined $isBed ) {
    $isBed = 0;
  }

  my $splicing_threshold;
  if($perform_splicing){
    $splicing_threshold = get_option( $config, $section, "splicing_threshold", 0 );
    if ( $splicing_threshold > 0 ) {
      $option = $option + " -splicing_threshold $splicing_threshold";
    }
  }

  my $refineSplicing = $perform_splicing && ($option =~ 'refGene');
  my $pythonSplicing = dirname(__FILE__) . "/annovarSplicing.py";
  if ( $refineSplicing & !-e $pythonSplicing ) {
    die "File not found : " . $pythonSplicing;
  }

  my $toExcel = get_option( $config, $section, "to_excel", 0 );
  my $affyFile = get_param_file( $config->{$section}{affy_file}, "affy_file", 0 );

  my $raw_files = get_raw_files( $config, $section );

  my $sampleCount = scalar keys %$raw_files;

  my $shfile;
  my $sh;
  if ( $sampleCount > 1 ) {
    $shfile = $self->get_task_filename( $pbs_dir, $task_name );
    open( $sh, ">$shfile" ) or die "Cannot create $shfile";
    print $sh get_run_command($sh_direct);
  }

  my $listfile = $self->get_file( $result_dir, $task_name, ".list", 0 );
  open( my $lt, ">$listfile" ) or die "Cannot create $listfile";

  for my $sample_name ( sort keys %{$raw_files} ) {
    my @sample_files = @{ $raw_files->{$sample_name} };

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log_file = $self->get_log_filename( $log_dir, $sample_name );

    my $cur_dir = create_directory_or_die($output_to_same_folder ? $result_dir : $result_dir . "/$sample_name");

    my $log_desc = $cluster->get_log_description($log_file);

    my $final_file = $self->get_final_file($config, $section, $cur_dir, $sample_name);
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );

    for my $sampleFile (@sample_files) {
      my ( $filename, $dir ) = fileparse($sampleFile);

      if ( $dir eq $cur_dir ) {
        $sampleFile = $filename;
      }

      my $cat = ( $filename =~ /.gz/ ) ? "zcat" : "cat";

      $filename =~ s/.gz$//g;
      $filename =~ s/.vcf$//g;

      my $annovar       = $filename . ".annovar";
      my $result        = "${annovar}.${buildver}_multianno.txt";
      my $refine_result = "${annovar}.splicing.${buildver}_multianno.txt";
      my $final_tsv         = $annovar . ".final.tsv";
      my $final = $annovar_gzipped ? $final_tsv . ".gz" : $final_tsv;
      my $excel         = $final . ".xls";

      my $runcmd;
      my $passinput;
      if ($isvcf) {
        $passinput = $filename . ".avinput";
        $runcmd    = "convert2annovar.pl -format vcf4old ${sampleFile} | cut -f1-7 | awk '{gsub(\",\\\\*\", \"\", \$0); print}'> $passinput 
  if [ -s $passinput ]; then
    table_annovar.pl $passinput $annovarDB $option --thread $thread --outfile $annovar --remove
  fi";
      }
      elsif ($isBed) {
        $passinput = $filename . ".avinput";
        $runcmd    = "perl -lane 'my \$fileColNum=scalar(\@F);my \$fileColPart=join(\"\t\",\@F[3..(\$fileColNum-1)]);print \"\$F[0]\t\$F[1]\t\$F[2]\t0\t-\t\$fileColPart\"' $sampleFile > $passinput 
  if [ -s $passinput ]; then
    table_annovar.pl $passinput $annovarDB $option --thread $thread --outfile $annovar --remove
  fi";
      }
      else {
        $passinput = $sampleFile;
        $runcmd    = "table_annovar.pl $passinput $annovarDB $option --thread $thread --outfile $annovar --remove";
      }

      print $pbs "
if [[ ! -s $result && ! -s $final ]]; then 
  $runcmd
fi
";

      if ($refineSplicing) {
        my $splicing_threshold_option = $splicing_threshold > 0 ? " -s $splicing_threshold" : "";
        print $pbs "
echo find_protein_position_for_splicing=`date`
if [[ -s $result && ! -s $refine_result ]]; then 
  python3 $pythonSplicing -i $result -d $annovarDB -o $refine_result -b $buildver $splicing_threshold_option
fi
";
        $result = $refine_result;
      }

      print $pbs "
if [[ -s $result && ! -s $final_tsv ]]; then
  sed -n '/^[^#]/q;p' ${sampleFile} | sed '\$ d' > ${final_tsv}.header
  $cat ${sampleFile} | grep -v \"^##\" | cut -f7- > ${filename}.clean
  grep -v \"^##\" ${result} > ${result}.clean
  paste ${result}.clean ${filename}.clean > ${final_tsv}.data
  cat ${final_tsv}.header ${final_tsv}.data > $final_tsv
  rm ${filename}.clean ${result}.clean ${final_tsv}.header ${final_tsv}.data
fi
";

      if ( $toExcel ) {
        my $affyoption = defined($affyFile) ? "-a $affyFile" : "";
        print $pbs "
if [[ -s $final_tsv && ! -s $excel ]]; then
  cqstools annovar_refine -i $final_tsv $affyoption -o $excel
fi
";
      }

      if($annovar_gzipped){
        print $pbs "
if [[ -s $final_tsv && ! -s $final ]]; then
  gzip $final_tsv
fi
";
      }

      print $lt "${cur_dir}/${result}\n";

      if($clean_folder){
        print $pbs "\nrm -f $passinput $result $refine_result \n";
      }
    }
    $self->close_pbs( $pbs, $pbs_file );

    if ( $sampleCount > 1 ) {
      print $sh "\$MYCMD ./$pbs_name \n";
    }
  }
  close $lt;

  if ( $sampleCount > 1 ) {
    print $sh "exit 0\n";
    close $sh;

    if ( is_linux() ) {
      chmod 0755, $shfile;
    }
    print "!!!shell file $shfile created, you can run this shell file to submit Annovar tasks.\n";
  }
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $buildver = $config->{$section}{buildver} or die "buildver is not defined in $section";
  my $raw_files = get_raw_files( $config, $section );
  my $toExcel = get_option( $config, $section, "to_excel", 0 );
  my $output_to_same_folder = get_option( $config, $section, "output_to_same_folder", 0 );
  my $annovar_gzipped = get_option( $config, $section, "annovar_gzipped", 0 );

  my $result = {};
  for my $sample_name ( sort keys %{$raw_files} ) {
    my @sample_files = @{ $raw_files->{$sample_name} };
    my $cur_dir      = $output_to_same_folder ? $result_dir : $result_dir . "/$sample_name";
    my @result_files = ();
    for my $sampleFile (@sample_files) {
      my ( $filename, $dir ) = fileparse($sampleFile);

      my $cat = ( $filename =~ /.gz/ ) ? "zcat" : "cat";

      $filename =~ s/.gz$//g;
      $filename =~ s/.vcf$//g;

      my $annovar = $filename . ".annovar";
      my $final   = $annovar . ".final.tsv";
      if ($annovar_gzipped){
        $final = $final . ".gz";
      }
      if ( $toExcel ) {
        my $excel = $final . ".xls";
        push( @result_files, $cur_dir . "/$excel" );
      }
      push( @result_files, $cur_dir . "/$final" );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
