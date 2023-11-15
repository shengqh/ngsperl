#!/usr/bin/perl
package Format::MergeFastqGzipped;

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
  $self->{_suffix} = "_mf";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $ispaired    = get_is_paired_end_option( $config, $section );
  my $is_collated = get_option( $config, $section, "is_collated", 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );

    if ($ispaired) {
      my $final_1_file  = $sample_name . ".1.fastq.gz";
      my $final_1_file_tmp  = $sample_name . ".1.fastq.tmp.gz";
      my $final_2_file  = $sample_name . ".2.fastq.gz";
      my $final_2_file_tmp  = $sample_name . ".2.fastq.tmp.gz";
      my $log_desc      = $cluster->get_log_description($log);
      my $pbs           = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_1_file );

    print $sh "
if [[ ! -s $result_dir/$final_1_file ]]; then
  \$MYCMD ./$pbs_name 
fi
";

      my $file_count = scalar(@sample_files);
      die "file count $file_count is not even for sample $sample_name @sample_files " if $file_count % 2 != 0;

      if ( scalar(@sample_files) == 2 ) {
        print $pbs "ln -s $sample_files[0] $final_1_file \n";
        print $pbs "ln -s $sample_files[1] $final_2_file \n";

        #print $pbs "cp $sample_files[0] $final_file \n";
      }
      else {
        print $pbs "
rm -f $final_1_file $final_1_file_tmp $final_2_file $final_2_file_tmp $sample_name.failed $sample_name.succeed
";
        if ($is_collated) {
          for ( my $sample_index = 0 ; $sample_index < $file_count / 2 ; $sample_index++ ) {
            my $curSample = $sample_files[$sample_index];
            print $pbs "
echo merging $sample_files[$sample_index] ...
cat $sample_files[$sample_index] >> $final_1_file_tmp 
";
          }
          for ( my $sample_index = $file_count / 2 ; $sample_index < $file_count ; $sample_index++ ) {
            my $curSample = $sample_files[$sample_index];
            print $pbs "
echo merging $sample_files[$sample_index] ...
cat $sample_files[$sample_index] >> $final_2_file_tmp
";
          }
        }
        else {
          for ( my $sample_index = 0 ; $sample_index < $file_count ; $sample_index += 2 ) {
            my $curSample = $sample_files[$sample_index];
            print $pbs "
echo merging $sample_files[$sample_index] ...
cat $sample_files[$sample_index] >> $final_1_file_tmp 

echo merging $sample_files[$sample_index+1] ...
cat $sample_files[$sample_index+1] >> $final_2_file_tmp 
";
          }
        }

        print $pbs "
status=\$?
if [[ \$status -ne 0 ]]; then
  rm -f $final_1_file $final_1_file_tmp $final_2_file $final_2_file_tmp
  touch $sample_name.failed
else
  mv $final_1_file_tmp $final_1_file 
  mv $final_2_file_tmp $final_2_file
  touch $sample_name.succeed
fi
";

      }
      $self->close_pbs( $pbs, $pbs_file );

    }
    else {
      my $final_file  = $sample_name . ".fastq.gz";
      my $log_desc    = $cluster->get_log_description($log);
      my $pbs         = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    print $sh "
if [[ ! -s $result_dir/$final_file ]]; then
  \$MYCMD ./$pbs_name 
fi
";

      if ( scalar(@sample_files) == 1 ) {
        print $pbs "ln -s $sample_files[0] $final_file \n";

        #print $pbs "cp $sample_files[0] $final_file \n";
      }
      else {
        print $pbs "if [ -s $final_file ]; then
  rm $final_file
fi
";
        for my $sample_file (@sample_files) {
          print $pbs "
echo merging $sample_file ...
cat $sample_file >> $final_file 
";
        }
      }
      $self->close_pbs( $pbs, $pbs_file );
    }
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all Bam2Fastq tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $ispaired    = get_is_paired_end_option( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    if ($ispaired) {
      my $final_1_file = $sample_name . ".1.fastq.gz";
      my $final_2_file = $sample_name . ".2.fastq.gz";
      push( @result_files, $result_dir . "/" . $final_1_file );
      push( @result_files, $result_dir . "/" . $final_2_file );
    }
    else {
      my $final_file = $sample_name . ".fastq.gz";
      push( @result_files, $result_dir . "/" . $final_file );
    }
    $result->{$sample_name} = filter_array( \@result_files, $pattern );

  }
  return $result;
}

1;
