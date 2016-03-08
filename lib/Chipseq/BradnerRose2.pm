#!/usr/bin/perl
package Chipseq::BradnerRose2;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::StringUtils;
use CQS::GroupTask;
use CQS::NGSCommon;
use Data::Dumper;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_br";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $genome = get_option( $config, $section, "genome" );
  my $pipeline_dir = get_directory( $config, $section, "pipeline_dir", 1 );

  if ( $option eq "" ) {
    $option = "-s 12500 -t 2500";
  }

  my %raw_files = %{ $self->get_grouped_raw_files( $config, $section, "groups" ) };

  #print Dumper(%raw_files) . "\n";

  my %control_files = %{ $self->get_grouped_raw_files( $config, $section, "controls" ) };

  #print Dumper(%control_files) . "\n";

  my %binding_site_beds = %{ get_raw_files( $config, $section, "binding_site_bed" ) };
  my $binding_site_filter = $config->{$section}{"binding_site_filter"};

  #print Dumper(%binding_site_beds) . "\n";

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $group_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$group_name} };
    my $treatment    = "-r " . $sample_files[0];

    my @controls = @{ $control_files{$group_name} };
    my $control  = "-c " . $controls[0];

    my @binding_files = @{ $binding_site_beds{$group_name} };

    my $cur_dir = create_directory_or_die( $result_dir . "/$group_name" );

    my $cpRawFile  = "${group_name}.copynumber";
    my $cpCallFile = "${group_name}.call";
    my $cpSegFile  = "${group_name}.segment";

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $group_name );

    my $log_desc = $cluster->get_log_description($log);
    if ( scalar(@binding_files) > 1 ) {
      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir );
      print $pbs "
cd $pipeline_dir
";

      for my $binding_file (@binding_files) {
        my ( $filename, $dirs, $suffix ) = fileparse( $binding_file, ".bed\$" );
        $filename =~ s/\./_/g;

        my $final_file = "${cur_dir}/${filename}_peaks_AllEnhancers.table.txt";
        print $pbs "if [ ! -s $final_file ]; then \n";

        my $newfile = "${cur_dir}/${filename}${suffix}";
        if ( defined $binding_site_filter ) {
          print $pbs "  grep \"$binding_site_filter\" $binding_file > $newfile \n";
        }
        else {
          print $pbs "  cp -f $binding_file $newfile \n";
        }
        print $pbs "  python ROSE2_main.py -g $genome -i $newfile $treatment $control -o $cur_dir $option 
fi

";
      }
      $self->close_pbs( $pbs, $pbs_file );
    }
    else {
      my ( $filename, $dirs, $suffix ) = fileparse( $binding_files[0], ".bed\$" );
      $filename =~ s/\./_/g;
      my $final_file = "${cur_dir}/${filename}_peaks_AllEnhancers.table.txt";

      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir, $final_file );
      print $pbs "
cd $pipeline_dir
";
      my $binding_file = $binding_files[0];

      my $newfile = "${cur_dir}/${filename}${suffix}";
      if ( defined $binding_site_filter ) {
        print $pbs "grep \"$binding_site_filter\" $binding_file > $newfile \n";
      }
      else {
        print $pbs "cp -f $binding_file $newfile \n";
      }
      print $pbs "python ROSE2_main.py -g $genome -i $newfile $treatment $control -o $cur_dir $option \n"; 
      $self->close_pbs( $pbs, $pbs_file );
    }

    print $sh "\$MYCMD ./$pbs_name \n";
  }

  print $sh "exit 0\n";
  close $sh;

  print "!!!shell file $shfile created, you can run this shell file to submit tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section );

  my %raw_files = %{ $self->get_grouped_raw_files( $config, $section, "groups" ) };
  my %binding_site_beds = %{ get_raw_files( $config, $section, "binding_site_bed" ) };

  my $result = {};
  for my $group_name ( sort keys %raw_files ) {
    my $cur_dir       = $result_dir . "/$group_name";
    my @result_files  = ();
    my @binding_files = @{ $binding_site_beds{$group_name} };
    for my $binding_file (@binding_files) {
      my ( $filename, $dirs, $suffix ) = fileparse( $binding_file, ".bed\$" );
      $filename =~ s/\./_/g;
      my $final_file = "${filename}_peaks_AllEnhancers.table.txt";
      push( @result_files, $cur_dir . "/${filename}_peaks_AllEnhancers.table.txt" );
    }
    $result->{$group_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
