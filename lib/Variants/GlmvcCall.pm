#!/usr/bin/perl
package Variants::GlmvcCall;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::GroupTask;
use CQS::NGSCommon;
use CQS::StringUtils;

our @ISA = qw(CQS::GroupTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_gc";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread ) = $self->init_parameter( $config, $section );

  my $glmvcfile = get_param_file( $config->{$section}{execute_file}, "execute_file", 1, not $self->using_docker() );
  my $source_type = $config->{$section}{source_type} or die "source_type is not defined in $section";

  my $rnaediting_db = get_directory( $config, $section, "rnaediting_db", 0 );
  if ( defined $rnaediting_db ) {
    $option = $option . " --rnaediting_db $rnaediting_db ";
  }

  my $annovar_buildver = $config->{$section}{annovar_buildver};
  if ( defined $annovar_buildver ) {
    $option = $option . " --annovar_buildver $annovar_buildver ";

    my $annovar_protocol = $config->{$section}{annovar_protocol};
    if ( defined $annovar_protocol ) {
      $option = $option . " --annovar_protocol $annovar_protocol ";
    }

    my $annovar_operation = $config->{$section}{annovar_operation};
    if ( defined $annovar_operation ) {
      $option = $option . " --annovar_operation $annovar_operation ";
    }

    my $annovar_db = $config->{$section}{annovar_db};
    if ( defined $annovar_db ) {
      $option = $option . " --annovar_db $annovar_db ";
    }
  }

  my $distance_exon_gtf = get_param_file( $config->{$section}{distance_exon_gtf}, "distance_exon_gtf", 0 );
  if ( defined $distance_exon_gtf ) {
    $option = $option . " --distance_exon_gtf $distance_exon_gtf ";
  }

  my $anno = defined $rnaediting_db || defined $annovar_buildver || defined $annovar_buildver;

  my $raw_files = get_raw_files( $config, $section );

  my %group_sample_map = ();

  my $fafile           = "";
  my $mpileupParameter = "";
  my $isbam            = lc($source_type) eq "bam";
  if ($isbam) {
    $fafile = get_param_file( $config->{$section}{fasta_file}, "fasta_file (for mpileup)", 1 );
    $mpileupParameter = $config->{$section}{mpileup_option};
    if ( defined $mpileupParameter ) {
      if ( $mpileupParameter eq "" ) {
        undef($$mpileupParameter);
      }
    }

    my $groups = get_raw_files( $config, $section, "groups" );
    for my $group_name ( sort keys %{$groups} ) {
      my @samples = @{ $groups->{$group_name} };
      my @gfiles  = ();
      my $index   = 0;
      foreach my $sample_name (@samples) {
        my @bam_files = @{ $raw_files->{$sample_name} };
        push( @gfiles, $bam_files[0] );
      }
      $group_sample_map{$group_name} = \@gfiles;
    }
  }
  else {
    %group_sample_map = %{$raw_files};
  }

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";
  print $sh "cd $pbs_dir \n";

  for my $group_name ( sort keys %group_sample_map ) {
    my @sample_files = @{ $group_sample_map{$group_name} };
    my $sampleCount  = scalar(@sample_files);
    my $cur_dir      = create_directory_or_die( $result_dir . "/$group_name" );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $group_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $group_name );

    print $sh "\$MYCMD ./$pbs_name \n";

    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir );

    my $finalvcf = "${group_name}.vcf";

    if ($isbam) {
      if ( $sampleCount != 2 ) {
        die "SampleFile should be normal,tumor paired.";
      }

      my $normal = $sample_files[0];
      my $tumor  = $sample_files[1];
      my $final  = $anno ? "${group_name}.annotation.tsv" : "${group_name}.tsv";

      my $cmd;
      if ( defined $mpileupParameter ) {
        $cmd = "samtools mpileup -f $fafile $mpileupParameter $normal $tumor | mono-sgen $glmvcfile call -t console $option -o ${cur_dir}/${group_name}";
      }
      else {
        $cmd = "mono $glmvcfile call -c $thread -t bam -f $fafile $option --normal $normal --tumor $tumor -o ${cur_dir}/${group_name}";
      }

      print $pbs "
if [ -s $final ]; then
  echo job has already been done. if you want to do again, delete ${cur_dir}/${final} and submit job again.
  exit 0;
fi      
      
if [ ! -s ${normal}.bai ]; then
  samtools index ${normal}
fi

if [ ! -s ${tumor}.bai ]; then
  samtools index ${tumor}
fi

$cmd

";
    }
    else {
      print $pbs "mono $glmvcfile call -t mpileup -m $sample_files[0] $option -o ${cur_dir}/${group_name} \n";
    }

    print $pbs "grep -v \"^#\" $finalvcf | cut -f1 | uniq -c | awk '{print \$2\"\\t\"\$1}' > ${finalvcf}.chromosome";
    $self->close_pbs( $pbs, $pbs_file );

  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all " . $self->{_name} . " tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $rnaediting_db     = get_directory( $config, $section, "rnaediting_db", 0 );
  my $annovar_buildver  = $config->{$section}{annovar_buildver};
  my $distance_exon_gtf = get_param_file( $config->{$section}{distance_exon_gtf}, "distance_exon_gtf", 0 );
  my $anno              = defined $rnaediting_db || defined $annovar_buildver || defined $annovar_buildver;

  my $groups = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $group_name ( keys %{$groups} ) {
    my @result_files = ();
    my $cur_dir      = $result_dir . "/$group_name";
    push( @result_files, "$cur_dir/${group_name}.tsv" );
    push( @result_files, "$cur_dir/${group_name}.vcf" );
    if ($anno) {
      push( @result_files, "$cur_dir/${group_name}.annotation.tsv" );
    }
    $result->{$group_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

sub get_pbs_files {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = $self->init_parameter( $config, $section, 0 );

  my $pairs = get_raw_files( $config, $section, "groups" );

  my $result = {};
  for my $pair_name ( sort keys %{$pairs} ) {
    my $pbs_name = $self->pbs_name($pair_name);
    my $pbs_file = $pbs_dir . "/$pbs_name";
    $result->{$pair_name} = $pbs_file;
  }
  return $result;
}

1;
