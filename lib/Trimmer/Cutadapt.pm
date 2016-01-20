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

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $random_bases_remove_after_trim = get_option($config, $section, "random_bases_remove_after_trim", 0);
  
  if ( defined $config->{$section}{adapter} ) {
    $option = $option . " -a " . $config->{$section}{adapter};
  }

  my $extension = get_option( $config, $section, "extension" );
  my $gzipped = get_option( $config, $section, "gzipped", 1 );

  if ( $gzipped && $extension =~ /\.gz$/ ) {
    $extension =~ s/\.gz$//g;
  }

  my $optionOnlyLimited='';
  my $optionRemoveLimited=$option;
  my $shortLimited = $option =~ /(-m\s+\d+\s+)/;
  if ($shortLimited) {
  	$shortLimited=$1;
  	$optionOnlyLimited=$optionOnlyLimited.$shortLimited;
  	$optionRemoveLimited=~s/$shortLimited//;
  }
  my $longLimited  = $option =~ /(-M\s+\d+\s+)/;
  if ($longLimited) {
  	$longLimited=$1;
  	$optionOnlyLimited=$optionOnlyLimited.$longLimited;
  	$optionRemoveLimited=~s/$longLimited//;
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

    my $final_file = scalar(@sample_files) == 1 ? $sample_name . $extension : $sample_name . ".1.fastq";
    if ($gzipped) {
      $final_file = $final_file . ".gz";
    }
    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    if ( scalar(@sample_files) == 1 ) { # single reads
      my $finalName      = $sample_name . $extension;
      my $finalShortName = $finalName . ".short";
      my $finalLongName  = $finalName . ".long";

      my $final_file     = $gzipped ? "${finalName}.gz"      : $finalName;
      my $finalShortFile = $gzipped ? "${finalShortName}.gz" : $finalShortName;
      my $finalLongFile  = $gzipped ? "${finalLongName}.gz"  : $finalLongName;

      if ($random_bases_remove_after_trim) { #remove top random bases
      	my $temp_file = $final_file . ".cutAdapter.fastq";
      	print $pbs "cutadapt $optionRemoveLimited -o $temp_file $sample_files[0]\n";
      	
      	print $pbs "cutadapt $optionOnlyLimited -u $random_bases_remove_after_trim -u -$random_bases_remove_after_trim -o $final_file ";
      	if ($shortLimited) {
           print $pbs " --too-short-output=$finalShortFile";
        }
        if ($longLimited) {
           print $pbs " --too-long-output=$finalLongFile";
        }
        print $pbs " $temp_file \n";
      } else { #NOT remove top random bases
      	 print $pbs "cutadapt $option -o $final_file ";
         if ($shortLimited) {
           print $pbs " --too-short-output=$finalShortFile";
         }
         if ($longLimited) {
           print $pbs " --too-long-output=$finalLongFile";
         }
         print $pbs " $sample_files[0] \n";
      }
    }
    else {
      #pair-end data
      my $read1file = $sample_files[0];
      my $read2file = $sample_files[1];
      my $read1name = $sample_name . ".1.fastq.gz";
      my $read2name = $sample_name . ".2.fastq.gz";

      if ( $shortLimited || $longLimited ) {
        my $temp1name = $sample_name . ".1.tmp.fastq";
        my $temp2name = $sample_name . ".2.tmp.fastq";
      	my $temp1_file = $read1name . ".cutAdapter.fastq";
      	my $temp2_file = $read2name . ".cutAdapter.fastq";
      	    
        #https://cutadapt.readthedocs.org/en/stable/guide.html#illumina-truseq
        if ($random_bases_remove_after_trim) { # remove top random bases
#           print $pbs "cutadapt $option -o $temp1name -p $temp2name $read1file $read2file \n";
#       	print $pbs "cutadapt $option -o $temp2_file -p $temp1_file $temp2name $temp1name \n";
#        	print $pbs "rm $temp2name $temp1name \n";

        	print $pbs "cutadapt $optionRemoveLimited -o $temp1_file $read1file \n";
        	print $pbs "cutadapt $optionRemoveLimited -o $temp2_file $read2file \n";
        	
        	print $pbs "cutadapt $optionOnlyLimited -u $random_bases_remove_after_trim -u -$random_bases_remove_after_trim -o $temp1name -p $temp2name $temp1_file $temp2_file \n";
        	print $pbs "cutadapt $optionOnlyLimited -u $random_bases_remove_after_trim -u -$random_bases_remove_after_trim -o $read2name -p $read1name $temp2name $temp1name \n";
        	print $pbs "rm $temp2name $temp1name \n";
        	
        } else { # NOT remove top random bases
            print $pbs "cutadapt $option -o $temp1name -p $temp2name $read1file $read2file \n";
        	print $pbs "cutadapt $option -o $read2name -p $read1name $temp2name $temp1name \n";
        	print $pbs "rm $temp2name $temp1name \n";
        }
     }
      else { #no short or long limited
      	if ($random_bases_remove_after_trim) { # remove top random bases
      	    my $temp1_file = $read1name . ".cutAdapter.fastq";
      	    my $temp2_file = $read2name . ".cutAdapter.fastq";
      	    print $pbs "cutadapt $option -o $temp1_file $read1file \n";
        	print $pbs "cutadapt $option -o $temp2_file $read2file \n";	
      	    print $pbs "cutadapt -u $random_bases_remove_after_trim -u -$random_bases_remove_after_trim -o $read1name $temp1_file \n";
      	    print $pbs "cutadapt -u $random_bases_remove_after_trim -u -$random_bases_remove_after_trim -o $read2name $temp2_file \n";
      	} else { # NOT remove top random bases
      		print $pbs "cutadapt $option -o $read1name $read1file \n";
        	print $pbs "cutadapt $option -o $read2name $read2file \n";	
      	}
      }
    }
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

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my $extension = $config->{$section}{extension} or die "define ${section}::extension first";
  my $gzipped = get_option( $config, $section, "gzipped", 1 );

  if ( $gzipped && $extension =~ /\.gz$/ ) {
    $extension =~ s/\.gz$//g;
  }

  my $shortLimited = $option =~ /-m\s+\d+/;
  my $longLimited  = $option =~ /-M\s+\d+/;

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};

  for my $sample_name ( keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };

    my @result_files = ();
    if ( scalar(@sample_files) == 1 ) {
      my $finalName      = $sample_name . $extension;
      my $finalShortName = $finalName . ".short";
      my $finalLongName  = $finalName . ".long";

      my $final_file     = $gzipped ? "${finalName}.gz"      : $finalName;
      my $finalShortFile = $gzipped ? "${finalShortName}.gz" : $finalShortName;
      my $finalLongFile  = $gzipped ? "${finalLongName}.gz"  : $finalLongName;

      push( @result_files, $result_dir . "/" . $final_file );

      if ($shortLimited) {
        push( @result_files, $result_dir . "/" . $finalShortFile );
      }
      if ($longLimited) {
        push( @result_files, $result_dir . "/" . $finalLongFile );
      }
    }
    else {

      #pair-end data
      my $read1name = $sample_name . ".1.fastq.gz";
      my $read2name = $sample_name . ".2.fastq.gz";

      push( @result_files, $result_dir . "/" . $read1name );
      push( @result_files, $result_dir . "/" . $read2name );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
