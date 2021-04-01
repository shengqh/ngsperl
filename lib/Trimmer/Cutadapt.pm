#!/usr/bin/perl
package Trimmer::Cutadapt;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_cut";
  bless $self, $class;
  return $self;
}

sub get_final_files {
  my ( $self, $ispairend, $sample_name, $extension, $fastqextension ) = @_;
  if ($ispairend) {
    my $finalName_1      = $sample_name . $extension . ".1" . $fastqextension;
    my $finalName_2      = $sample_name . $extension . ".2" . $fastqextension;
    my $final_file_1     = "${finalName_1}.gz";
    my $final_file_2     = "${finalName_2}.gz";
    my $finalShortFile_1 = $finalName_1 . ".short.gz";
    my $finalShortFile_2 = $finalName_2 . ".short.gz";
    my $finalLongFile_1  = $finalName_1 . ".long.gz";
    my $finalLongFile_2  = $finalName_2 . ".long.gz";
    return ( $final_file_1, $final_file_2, $finalShortFile_1, $finalShortFile_2, $finalLongFile_1, $finalLongFile_2 );
  }
  else {
    my $finalName      = $sample_name . $extension . $fastqextension;
    my $final_file     = "${finalName}.gz";
    my $finalShortFile = $finalName . ".short.gz";
    my $finalLongFile  = $finalName . ".long.gz";
    return ( $final_file, $finalShortFile, $finalLongFile );
  }
}

sub get_extension {
  my ( $self, $config, $section ) = @_;

  my $curSection = get_config_section( $config, $section );

  my $extension = $curSection->{extension} or die "define ${section}::extension first";
  if ( $extension =~ /\.gz$/ ) {
    $extension =~ s/\.gz$//g;
  }

  my $fastqextension = ".fastq";
  if ( $extension =~ /\.fastq$/ ) {
    $extension =~ s/\.fastq$//g;
  }

  if ( $extension =~ /\.fq$/ ) {
    $fastqextension = ".fq";
    $extension =~ s/\.fq$//g;
  }

  return ( $extension, $fastqextension );
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = $self->init_parameter( $config, $section );

  my $curSection = get_config_section( $config, $section );

  my $thread_option = "";
  if ($thread > 1) {
    $thread_option = "-j $thread";
  }

  my $hard_trim = get_option( $config, $section, "hard_trim", 0 );

  my $remove_bases_option = "";
  my $random_bases_remove_after_trim = get_option( $config, $section, "random_bases_remove_after_trim", 0 );
  if ( $random_bases_remove_after_trim > 0 ) {
    $remove_bases_option = "-u $random_bases_remove_after_trim -u -$random_bases_remove_after_trim";
  }
  else {
    my $random_bases_remove_after_trim_5 = get_option( $config, $section, "random_bases_remove_after_trim_5", 0 );
    if ( $random_bases_remove_after_trim_5 > 0 ) {
      $remove_bases_option = "-u $random_bases_remove_after_trim_5";
    }
    my $random_bases_remove_after_trim_3 = get_option( $config, $section, "random_bases_remove_after_trim_3", 0 );
    if ( $random_bases_remove_after_trim_3 > 0 ) {
      $remove_bases_option = $remove_bases_option . "-u -$random_bases_remove_after_trim_3";
    }
  }

  my $trim_base_quality_after_adapter_trim = get_option( $config, $section, "trim_base_quality_after_adapter_trim", 0 );
  if($trim_base_quality_after_adapter_trim > 0 && $remove_bases_option eq ""){
    $remove_bases_option = "-q " . $trim_base_quality_after_adapter_trim;
  }

  my $ispairend = get_is_paired_end_option( $config, $section );

  my $adapter_option = "";
  if ( $option !~ /-a/ ) {
    if ( defined $curSection->{adapter} && length( $curSection->{adapter} ) > 0 ) {
      my @adapters = split(',', $curSection->{adapter});
      print(@adapters);
      if ($ispairend) {
        $adapter_option = " -a " . join(' -a ', @adapters) . " -A " . join(' -A ', @adapters);
      }
      else {
        $adapter_option = " -a " . join(' -a ', @adapters);
      }
    }

    if ( defined $curSection->{adapter_3} && length( $curSection->{adapter_3} ) > 0 ) {
      my @adapters = split /,/, $curSection->{adapter_3};
      if ($ispairend) {
        $adapter_option = " -a " . join(' -a ', @adapters) . " -A " . join(' -A ', @adapters);
      }
      else {
        $adapter_option = " -a " . join(' -a ', @adapters);
      }
    }
  }

  if ( $option !~ /-g/ ) {
    if ( defined $curSection->{adapter_5} && length( $curSection->{adapter_5} ) > 0 ) {
      my @adapters = split /,/, $curSection->{adapter_5};
      if ($ispairend) {
        $adapter_option = $adapter_option . " -g " . join(' -g ', @adapters) . " -G " . join(' -G ', @adapters);
      }
      else {
        $adapter_option = $adapter_option . " -g " .  join(' -g ', @adapters);
      }
    }
  }

  my $trim_poly_atgc = get_option( $config, $section, "trim_poly_atgc", 1 );
  if ($trim_poly_atgc) {
    if ( $adapter_option =~ /-a/ ) {
      $adapter_option = $adapter_option . " -a \"A{50}\" -a \"T{50}\" -a \"G{50}\" -a \"C{50}\"";
    }

    if ( $adapter_option =~ /-g/ ) {
      $adapter_option = $adapter_option . " -g \"A{50}\" -g \"T{50}\" -g \"G{50}\" -g \"C{50}\"";
    }

    if ($ispairend) {
      if ( $adapter_option =~ /-A/ ) {
        $adapter_option = $adapter_option . " -A \"A{50}\" -A \"T{50}\" -A \"G{50}\" -A \"C{50}\"";
      }

      if ( $adapter_option =~ /-G/ ) {
        $adapter_option = $adapter_option . " -G \"A{50}\" -G \"T{50}\" -G \"G{50}\" -G \"C{50}\"";
      }
    }
  }

  my ( $extension, $fastqextension ) = $self->get_extension( $config, $section );

  my $optionOnlyLimited   = '';
  my $optionRemoveLimited = $option;
  my $shortLimited        = $option =~ /(-m\s+\d+\s*)/;
  if ($shortLimited) {
    $shortLimited      = $1;
    $optionOnlyLimited = $optionOnlyLimited . " " . $shortLimited;
    $optionRemoveLimited =~ s/$shortLimited//;
  }
  my $longLimited = $option =~ /(-M\s+\d+\s*)/;
  if ($longLimited) {
    $longLimited       = $1;
    $optionOnlyLimited = $optionOnlyLimited . " " . $longLimited;
    $optionRemoveLimited =~ s/$longLimited//;
  }

  if ( index( $option, "--trim-n" ) == -1 ) {
    $option = $option . " --trim-n";
  }

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $sample_name ( sort keys %raw_files ) {

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my @sample_files = @{ $raw_files{$sample_name} };

    my @final_files = $self->get_final_files( $ispairend, $sample_name, $extension, $fastqextension );
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_files[0] );

    my @rmlist = ();
    for my $sample_file (@sample_files) {
      print $pbs "
if [ ! -s $sample_file ]; then
  echo input file not exist: $sample_file
  exit 0
fi
";
    }

    if ($hard_trim > 0){
      for my $sample_file (@sample_files) {
        my $temp_file = basename($sample_file) . ".hardtrim.fastq";
        print $pbs "
cutadapt $thread_option -l $hard_trim -o $temp_file $sample_file
    ";
        push @rmlist, $temp_file;
      }
      @sample_files = @rmlist;
    }

    if ($ispairend) {
      die "should be pair-end data but not! " .  $sample_name . ": " . "; ".join(@sample_files) if ( scalar(@sample_files) != 2 );

      #pair-end data
      my $read1file = $sample_files[0];
      my $read2file = $sample_files[1];

      my ( $read1name, $read2name, $finalShortFile_1, $finalShortFile_2, $finalLongFile_1, $finalLongFile_2 ) = $self->get_final_files( $ispairend, $sample_name, $extension, $fastqextension );

      my $limit_file_options = "";
      if ($shortLimited) {
        $limit_file_options = " --too-short-output=$finalShortFile_1 --too-short-paired-output=$finalShortFile_2";
      }
      if ($longLimited) {
        $limit_file_options = $limit_file_options . " --too-long-output=$finalLongFile_1 --too-long-paired-output=$finalLongFile_2";
      }

      if ($remove_bases_option ne "") {    # remove top random bases
        my $temp1_file = $read1name . ".cutAdapter.fastq";
        my $temp2_file = $read2name . ".cutAdapter.fastq";
        print $pbs "
cutadapt $thread_option $option $adapter_option -o $temp1_file -p $temp2_file $read1file $read2file
cutadapt $thread_option $remove_bases_option -o $read1name -p $read2name $limit_file_options $temp1_file $temp2_file 
rm $temp1_file $temp2_file
";
      }
      else {                         # NOT remove top random bases
        print $pbs "
cutadapt $thread_option $option $adapter_option -o $read1name -p $read2name $limit_file_options $read1file $read2file
";
      }
    }
    else {
      my ( $final_file, $finalShortFile, $finalLongFile ) = $self->get_final_files( $ispairend, $sample_name, $extension, $fastqextension );

      my $limit_file_options = "";
      if ($shortLimited) {
        $limit_file_options = " --too-short-output=$finalShortFile";
      }
      if ($longLimited) {
        $limit_file_options = $limit_file_options . " --too-long-output=$finalLongFile";
      }

      if ($remove_bases_option ne "") {    #remove top random bases
        my $temp_file = $final_file . ".cutAdapter.fastq";
        if ( scalar(@sample_files) == 1 ) {
          print $pbs "
cutadapt $thread_option $optionRemoveLimited $adapter_option -o $temp_file $sample_files[0]
cutadapt $thread_option $optionOnlyLimited $limit_file_options $remove_bases_option -o $final_file $temp_file
rm $temp_file
";
        }
        else {
          print $pbs "
if [ -s $temp_file ]; then
  delete $temp_file
fi
";
          for my $sample_file (@sample_files) {
            print $pbs "
cutadapt $thread_option $optionRemoveLimited $adapter_option -o temp.fastq $sample_file
cat temp.fastq >> $temp_file
rm temp.fastq
";
          }

          print $pbs "
cutadapt $thread_option $optionOnlyLimited $limit_file_options $remove_bases_option -o $final_file $temp_file 
rm $temp_file
";
        }
      }
      else {    #NOT remove top random bases
        if ( scalar(@sample_files) == 1 ) {
          print $pbs "
cutadapt $thread_option $option $adapter_option $limit_file_options -o $final_file $sample_files[0]
";
        }
        else {
          my $temp_file = $final_file . ".cutAdapter.fastq";
          print $pbs "
if [ -s $temp_file ]; then
  delete $temp_file
fi
";
          if ( length($limit_file_options) > 0 ) {
            for my $sample_file (@sample_files) {
              print $pbs "
cutadapt $thread_option $optionRemoveLimited $adapter_option -o temp.fastq $sample_file
cat temp.fastq >> $temp_file
rm temp.fastq
";
            }

            print $pbs "
cutadapt $thread_option $optionOnlyLimited $limit_file_options -o $final_file $temp_file
rm $temp_file
";
          }
          else {
            for my $sample_file (@sample_files) {
              print $pbs "
cutadapt $thread_option $option $adapter_option -o temp.fastq $sample_file
cat temp.fastq >> $temp_file
rm temp.fastq
";
            }
          }
          print $pbs "
gzip $temp_file
mv ${temp_file}.gz $final_file
";
        }
      }
    }

    if(scalar(@rmlist) > 0){
      my $rmstr = join(" ",  @rmlist);
      print $pbs "
rm $rmstr
";
    }

    print $pbs "
cutadapt --version 2>&1 | awk '{print \"Cutadapt,v\"\$1}' > ${sample_name}.version
";
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $ispairend = get_is_paired_end_option( $config, $section );
  my ( $extension, $fastqextension ) = $self->get_extension( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};

  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    if ($ispairend) {
      my ( $read1name, $read2name ) = $self->get_final_files( $ispairend, $sample_name, $extension, $fastqextension );

      push( @result_files, $result_dir . "/" . $read1name );
      push( @result_files, $result_dir . "/" . $read2name );
    }
    else {
      my ( $final_file, $finalShortFile, $finalLongFile ) = $self->get_final_files( $ispairend, $sample_name, $extension, $fastqextension );

      push( @result_files, $result_dir . "/" . $final_file );

      my $shortLimited = $option =~ /-m\s+\d+/;
      my $longLimited  = $option =~ /-M\s+\d+/;
      if ($shortLimited) {
        push( @result_files, $result_dir . "/" . $finalShortFile );
      }
      if ($longLimited) {
        push( @result_files, $result_dir . "/" . $finalLongFile );
      }
    }
    push( @result_files, $result_dir . "/${sample_name}.version" );
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
