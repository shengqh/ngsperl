#!/usr/bin/perl
package Format::MergeTagAlign;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::GroupTask;
use CQS::StringUtils;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_mt";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %groups = %{ $self->get_grouped_raw_files( $config, $section, "groups" ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $group_name ( sort keys %groups ) {
    my @sample_files = @{ $groups{$group_name} };

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $group_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $final_tagalign  = $group_name . ".tagalign";
    my $final_file  = $group_name . ".tagalign.gz";
    my $log_desc    = $cluster->get_log_description($log);
    my $pbs         = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

    if ( scalar(@sample_files) == 1 ) {
      print $pbs "cp $sample_files[0] $final_file \n";
    }
    else {
      print $pbs "if [ -s $final_tagalign ]; then
  rm $final_tagalign
fi
";
      for my $sample_file (@sample_files) {
        print $pbs "
echo merging $sample_file ...
zcat $sample_file >> $final_tagalign 
";
      }
      print $pbs "
echo gzipping $final_tagalign ...
gzip $final_tagalign \n";
    }
    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %groups = %{ $self->get_grouped_raw_files( $config, $section, "groups" ) };

  my $result = {};
  for my $group_name ( keys %groups ) {
    my $final_file  = $group_name . ".tagAlign.gz";

    my @result_files = ();
    push( @result_files, $result_dir . "/" . $final_file );
    $result->{$group_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
