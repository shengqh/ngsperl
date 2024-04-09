#!/usr/bin/perl
package Visualization::Bamplot;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use File::Copy;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::UniqueTask;
use CQS::StringUtils;

our @ISA = qw(CQS::UniqueTask);

my $directory;

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_bp";
  $self->{_docker_prefix} = "bamplot_";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = $self->init_parameter( $config, $section );

  my $bam_files = get_raw_files( $config, $section );
  my $draw_individual = get_option( $config, $section, "is_draw_individual" );
  my $rainbow_color   = get_option( $config, $section, "is_rainbow_color" );
  my $default_color   = get_option( $config, $section, "default_color", "0,0,0" );
  my $gff_file = parse_param_file( $config, $section, "gff_file", 1 );

  my $draw_by_r        = get_option( $config, $section, "draw_by_r",        1 );
  my $draw_by_r_width  = get_option( $config, $section, "draw_by_r_width",  0 );
  my $draw_by_r_height = get_option( $config, $section, "draw_by_r_height", 0 );

  if ($option !~ /-y/){
    $option = $option . " -y uniform";
  }

  my $colors;
  if ( has_raw_files( $config, $section, "colors" ) ) {
    $colors = get_raw_files( $config, $section, "colors" );
  }

  my $groups;
  if ( has_raw_files( $config, $section, "groups" ) ) {
    $groups = get_raw_files( $config, $section, "groups" );
  }
  else {
    my @bamnames = keys %$bam_files;
    $groups = { $task_name => \@bamnames };
  }

  my $pbs_file = $self->get_pbs_filename( $pbs_dir, $task_name );
  my $pbs_name = basename($pbs_file);
  my $log      = $self->get_log_filename( $log_dir, $task_name );
  my $log_desc = $cluster->get_log_description($log);

  my $final_file = "${task_name}_plots.pdf";
  my $pbs      = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );

  for my $name ( sort keys %{$groups} ) {
    my @curbam_names = @{ $groups->{$name} };
    my @curbam_files = ();
    my @bam_colors   = ();
    for my $bamName (@curbam_names) {
      push( @curbam_files, $bam_files->{$bamName}->[0] );
      if ($colors) {
        if ( $colors->{$bamName} ) {
          push( @bam_colors, $colors->{$bamName}[0] );
        }
        else {
          push( @bam_colors, $default_color );
        }
      }
      else {
        push( @bam_colors, $default_color );
      }
    }

    for my $bam (@curbam_files) {
      print $pbs "if [ ! -s ${bam}.bai ]; then
  samtools index $bam
fi
";
    }

    if ($draw_individual) {
      for my $i ( 0 .. $#curbam_names ) {
        my $bamName  = $curbam_names[$i];
        my $bamFile  = $curbam_files[$i];
        my $bamColor = $bam_colors[$i];

        my $curgff = "${bamName}.gff";

        #copy( $gff_file, "${pbs_dir}/${curgff}" );

        my $colorStr = $rainbow_color ? "" : "--color $bamColor";

        print $pbs "
cp -f $gff_file $curgff
bamPlot_turbo $option -b \"$bamFile\" -n \"$bamName\" -i $curgff $colorStr -o .

";
      }
    }
    else {
      my $curgff = "${name}.gff";

      my $curbam_nameStr = join( ',', @curbam_names );
      my $curbam_fileStr = join( ',', @curbam_files );
      my $colorStr = $rainbow_color ? "" : "--color " . join( ':', @bam_colors );

      my $cur_final_pdf  = "${name}.pdf";
      my $cur_final_png  = "${name}.png";
      print $pbs "

cp -f $gff_file $curgff
bamPlot_turbo $option -b \"$curbam_fileStr\" -n \"$curbam_nameStr\" -i $curgff $colorStr -o .

";

      if ($draw_by_r) {
        my $rscript = dirname(__FILE__) . "/bamPlot_turbo.R";
        my $summaryfile = "${name}/${name}_summary.txt";
        print $pbs "
if [ -s $summaryfile ]; then
  #R --vanilla -f $rscript --args $summaryfile $cur_final_pdf UNIFORM MULTIPLE MULTIPLE_PAGE $draw_by_r_width $draw_by_r_height
  R --vanilla -f $rscript --args $summaryfile $cur_final_png UNIFORM MULTIPLE MULTIPLE_PAGE $draw_by_r_width $draw_by_r_height
fi
";
      }
    }
  }

  $self->close_pbs( $pbs, $pbs_file );

  print "!!!pbs file $pbs_file created, you can submit this file or run directly.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $bam_files = get_raw_files( $config, $section );
  my $groups;
  if ( has_raw_files( $config, $section, "groups" ) ) {
    $groups = get_raw_files( $config, $section, "groups" );
  }
  else {
    my @bamnames = keys %$bam_files;
    $groups = { $task_name => \@bamnames };
  }

  my $draw_individual = get_option( $config, $section, "is_draw_individual" );

  my $result       = {};
  my @result_files = ();
  for my $name ( sort keys %{$groups} ) {
    if ($draw_individual) {
      my @curbam_names = @{ $groups->{$name} };
      for my $bamName (@curbam_names) {
        push( @result_files, "${result_dir}/${bamName}.pdf" );
      }
    }
    else {
      push( @result_files, "${result_dir}/${name}_plots.pdf" );
    }
  }
  $result->{$task_name} = filter_array( \@result_files, $pattern );
  return $result;
}

1;
