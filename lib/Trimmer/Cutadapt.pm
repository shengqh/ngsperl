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
    return ( $sample_name . $extension . ".1" . $fastqextension . ".gz", $sample_name . $extension . ".2" . $fastqextension . ".gz" );
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $curSection = get_config_section( $config, $section );

  my $random_bases_option = "";
  my $random_bases_remove_after_trim = get_option( $config, $section, "random_bases_remove_after_trim", 0 );
  if ( $random_bases_remove_after_trim > 0 ) {
    $random_bases_option = "-u $random_bases_remove_after_trim -u -$random_bases_remove_after_trim";
  }
  else {
    my $random_bases_remove_after_trim_5 = get_option( $config, $section, "random_bases_remove_after_trim_5", 0 );
    if ( $random_bases_remove_after_trim_5 > 0 ) {
      $random_bases_option = "-u $random_bases_remove_after_trim_5";
    }
    my $random_bases_remove_after_trim_3 = get_option( $config, $section, "random_bases_remove_after_trim_3", 0 );
    if ( $random_bases_remove_after_trim_3 > 0 ) {
      $random_bases_option = $random_bases_option . "-u -$random_bases_remove_after_trim_3";
    }
  }

  my $ispairend = get_is_paired_end_option( $config, $section );

  my $adapter_option = "";
  if ( $option !~ /-a/ ) {
    if ( defined $curSection->{adapter} && length( $curSection->{adapter} ) > 0 ) {
      if ($ispairend) {
        $adapter_option = " -a " . $curSection->{adapter} . " -A " . $curSection->{adapter};
      }
      else {
        $adapter_option = " -a " . $curSection->{adapter};
      }
    }

    if ( defined $curSection->{adapter_3} && length( $curSection->{adapter_3} ) > 0 ) {
      if ($ispairend) {
        $adapter_option = " -a " . $curSection->{adapter_3} . " -A " . $curSection->{adapter_3};
      }
      else {
        $adapter_option = " -a " . $curSection->{adapter_3};
      }
    }
  }

  if ( $option !~ /-g/ ) {
    if ( defined $curSection->{adapter_5} && length( $curSection->{adapter_5} ) > 0 ) {
      if ($ispairend) {
        $adapter_option = $adapter_option . " -g " . $curSection->{adapter_5} . " -G " . $curSection->{adapter_5};
      }
      else {
        $adapter_option = $adapter_option . " -g " . $curSection->{adapter_5};
      }
    }
  }
 
  my $trim_poly_atgc = get_option( $config, $section, "trim_poly_atgc", 0 );
  if ($trim_poly_atgc) {
    if($adapter_option =~ /-a/){
      $adapter_option = $adapter_option . " -a \"A{50}\" -a \"T{50}\" -a \"G{50}\" -a \"C{50}\"";
    }

    if($adapter_option =~ /-g/){
      $adapter_option = $adapter_option . " -g \"A{50}\" -g \"T{50}\" -g \"G{50}\" -g \"C{50}\"";
    }

    if ($ispairend) {
      if($adapter_option =~ /-A/){
        $adapter_option = $adapter_option . " -A \"A{50}\" -A \"T{50}\" -A \"G{50}\" -A \"C{50}\"";
      }

      if($adapter_option =~ /-G/){
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

    for my $sample_file (@sample_files) {
      print $pbs "
if [ ! -s $sample_file ]; then
  echo input file not exits: $sample_file
  exit 0
fi
";

    }
    if ($ispairend) {
      die "should be pair-end data but not!" if ( scalar(@sample_files) != 2 );

      #pair-end data
      my $read1file = $sample_files[0];
      my $read2file = $sample_files[1];

      my ( $read1name, $read2name ) = $self->get_final_files( $ispairend, $sample_name, $extension, $fastqextension );

      if ($random_bases_option) {    # remove top random bases
        my $temp1_file = $read1name . ".cutAdapter.fastq";
        my $temp2_file = $read2name . ".cutAdapter.fastq";
        print $pbs "
cutadapt $option $adapter_option -o $temp1_file -p $temp2_file $read1file $read2file
cutadapt $random_bases_option -o $read1name -p $read2name $temp1_file $temp2_file
rm $temp1_file $temp2_file
";
      }
      else {                         # NOT remove top random bases
        print $pbs "
cutadapt $option $adapter_option -o $read1name -p $read2name $read1file $read2file
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

      if ($random_bases_option) {    #remove top random bases
        my $temp_file = $final_file . ".cutAdapter.fastq";
        if ( scalar(@sample_files) == 1 ) {
          print $pbs "
cutadapt $optionRemoveLimited $adapter_option -o $temp_file $sample_files[0]
cutadapt $optionOnlyLimited $limit_file_options $random_bases_option -o $final_file $temp_file
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
cutadapt $optionRemoveLimited $adapter_option -o temp.fastq $sample_file
cat temp.fastq >> $temp_file
rm temp.fastq
";
          }

          print $pbs "
cutadapt $optionOnlyLimited $limit_file_options $random_bases_option -o $final_file $temp_file 
rm $temp_file
";
        }
      }
      else {    #NOT remove top random bases
        if ( scalar(@sample_files) == 1 ) {
          print $pbs "
cutadapt $option $adapter_option $limit_file_options -o $final_file $sample_files[0]
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
cutadapt $optionRemoveLimited $adapter_option -o temp.fastq $sample_file
cat temp.fastq >> $temp_file
rm temp.fastq
";
            }

            print $pbs "
cutadapt $optionOnlyLimited $limit_file_options -o $final_file $temp_file
rm $temp_file
";
          }
          else {
            for my $sample_file (@sample_files) {
              print $pbs "
cutadapt $option $adapter_option -o temp.fastq $sample_file
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
    print $pbs "
cutadapt --version | awk '{print \"Cutadapt,v\"\$1}' > ${sample_name}.version
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

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
