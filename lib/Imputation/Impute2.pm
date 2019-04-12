#!/usr/bin/perl
package Imputation::Impute2;

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
use GWAS::GwasUtils;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_imp";
  bless $self, $class;
  return $self;
}

sub containPosition {
  my ( $positions, $start, $end ) = @_;
  for my $position ( @{$positions} ) {
    if ( $position >= $start && $position <= $end ) {
      return 1;
    }
  }

  return 0;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $isPhased = get_option( $config, $section, "isPhased", 0 );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  my $ranges = {};
  my $interval_bed = get_option_file( $config, $section, "range_file", "interval_bed" );
  open( my $fin, '<', $interval_bed ) or die "Could not open $interval_bed\n";
  while ( my $line = <$fin> ) {
    chomp $line;
    my @parts = split( "\t", $line );
    my $chrom = $parts[0];
    if (not defined $ranges->{$chrom}){
      $ranges->{$chrom} = [];
    }
    
    my $locus = $ranges->{$chrom};
    push @$locus, [$parts[1], $parts[2]];
  }

  my $raw_files  = get_raw_files( $config, $section );
  my $mapFiles = readGwasDataFile($config, $section, "genetic_map_file",$raw_files  );
  my $haploFiles = readGwasDataFile( $config, $section, "haplo_file", $raw_files );
  my $legendFiles = readGwasDataFile( $config, $section, "legend_file", $raw_files );
  
  for my $sample_name ( sort keys %$raw_files ) {
    my ( $chrName ) = $sample_name =~ /_chr(\S+)$/;
    if (not defined $ranges->{$chrName}){
      die "sample = $sample_name, chr = $chrName, not in $interval_bed";
    }
    my $locus = $ranges->{$chrName};
    
    my @sample_files = @{ $raw_files->{$sample_name} };
    my $sample       = $sample_files[0];

    my @mFiles = @{ $mapFiles->{$sample_name} };
    my $map    = $mFiles[0];

    my @hFiles     = @{ $haploFiles->{$sample_name} };
    my $haploFile  = $hFiles[0];

    my @lFiles     = @{ $legendFiles->{$sample_name} };
    my $legendFile = $lFiles[0];

    my $cur_dir = create_directory_or_die( $result_dir . "/$sample_name" );

    my $gen_file     = "${sample_name}.gen";
    my $gen_tmp_file = "${sample_name}.gen.tmp";

    for my $loc (@$locus) {
      my $start = $loc->[0];
      my $end   = $loc->[1];

      my $cursample = $sample_name . "_" . $start . "_" . $end;

      my $pbs_file = $self->get_pbs_filename( $pbs_dir, $cursample );
      my $pbs_name = basename($pbs_file);
      my $log      = $self->get_log_filename( $log_dir, $cursample );

      my $tmpFile = "${cursample}.tmp";

      my $log_desc = $cluster->get_log_description($log);

      my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir );

      print $pbs "

if [ -s $gen_file ]; then
  echo job has already been done. if you want to do again, delete $cur_dir/$gen_file and submit job again.
  exit 0;
fi

if [ -s $tmpFile ]; then
  echo job has already been done. if you want to do again, delete $cur_dir/$tmpFile and submit job again.
  exit 0;
fi

impute2 $option -known_haps_g $sample -m $map -int $start $end -h $haploFile -l $legendFile -o $tmpFile
";
      $self->close_pbs( $pbs, $pbs_file );
      print $sh "\$MYCMD ./$pbs_name \n";
    }

    close(INFILE);
  }
  print $sh "exit 0\n";
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my %rangeFiles = %{ get_raw_files( $config, $section, "range_file" ) };

  my $result = {};
  for my $sample_name ( sort keys %raw_files ) {
    my @result_files = ();

    my @sample_files = @{ $raw_files{$sample_name} };
    my $sample       = $sample_files[0];

    my @rFiles    = @{ $rangeFiles{$sample_name} };
    my $rangeFile = $rFiles[0];

    open( INFILE, "<", $rangeFile ) or die("Couldn't open $rangeFile for reading!\n");
    my $isfirst = 1;
    while (<INFILE>) {
      my $start = ( split( /\s+/, $_ ) )[0];
      my $end   = ( split( /\s+/, $_ ) )[1];

      my $cursample = $sample_name . "_" . $start . "_" . $end;

      my $tmpFile = "${cursample}.tmp";

      push( @result_files, "${result_dir}/${sample_name}/${tmpFile}" );
    }
    close(INFILE);

    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
