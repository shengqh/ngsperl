#!/usr/bin/perl
package TCGA::SplitSNP;

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
  $self->{_suffix} = "_ss";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster, $thread, $memory ) = get_parameter( $config, $section );
  $self->{_task_prefix} = get_option( $config, $section, "prefix", "" );
  $self->{_task_suffix} = get_option( $config, $section, "suffix", "" );

  my $prefix = get_option( $config, $section, "prefix", "" );
  my $plink_prefixes = get_raw_files( $config, $section );

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  for my $sample_name ( sort keys %$plink_prefixes ) {
    my $plink_prefix = $plink_prefixes->{$sample_name}->[0];
    my $fam_file     = $plink_prefix . ".fam";

    my $file_name = $prefix . $sample_name;

    my $cur_dir  = create_directory_or_die( $result_dir . "/" . $sample_name );
    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $log_desc = $cluster->get_log_description($log);

    print $sh "\$MYCMD ./$pbs_name \n";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $cur_dir );

    print $pbs "
if [ -s ${file_name}_Primary_Solid_Tumor.bed ]; then
  echo $sample_name has been prepared. Delete $cur_dir/${file_name}_Primary_Solid_Tumor.bed to regenerate the result.
  exit 1;
fi
    
grep \"\\-01\\s\" $fam_file > Primary_Solid_Tumor.fam
plink -bfile $plink_prefix --keep Primary_Solid_Tumor.fam --out ${file_name}_Primary_Solid_Tumor --make-bed 

grep \"\\-10\\s\" $fam_file > Blood_Derived_Normal.fam
plink -bfile $plink_prefix --keep Blood_Derived_Normal.fam --out ${file_name}_Blood_Derived_Normal --make-bed 

grep \"\\-11\\s\" $fam_file > Solid_Tissue_Normal.fam
plink -bfile $plink_prefix --keep Solid_Tissue_Normal.fam --out ${file_name}_Solid_Tissue_Normal --make-bed

rm Primary_Solid_Tumor.fam Blood_Derived_Normal.fam Solid_Tissue_Normal.fam *.log *.nosex
 
";

    $self->close_pbs( $pbs, $pbs_file );
  }

  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all fastx_trimmer tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my %raw_files = %{ get_raw_files( $config, $section ) };
  my $prefix = get_option( $config, $section, "prefix", "" );

  my $result = {};
  for my $sample_name ( keys %raw_files ) {
    my @result_files = ();
    my $cur_dir      = $result_dir . "/" . $sample_name;

    my $final_name = $prefix . $sample_name;
    push( @result_files, "${cur_dir}/${final_name}_Primary_Solid_Tumor.bed" );
    push( @result_files, "${cur_dir}/${final_name}_Primary_Solid_Tumor.bim" );
    push( @result_files, "${cur_dir}/${final_name}_Primary_Solid_Tumor.fam" );
    push( @result_files, "${cur_dir}/${final_name}_Blood_Derived_Normal.bed" );
    push( @result_files, "${cur_dir}/${final_name}_Blood_Derived_Normal.bim" );
    push( @result_files, "${cur_dir}/${final_name}_Blood_Derived_Normal.fam" );
    push( @result_files, "${cur_dir}/${final_name}_Solid_Tissue_Normal.bed" );
    push( @result_files, "${cur_dir}/${final_name}_Solid_Tissue_Normal.bim" );
    push( @result_files, "${cur_dir}/${final_name}_Solid_Tissue_Normal.fam" );
    
    $result->{$sample_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
