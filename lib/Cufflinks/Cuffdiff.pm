#!/usr/bin/perl
package Cufflinks::Cuffdiff;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::PairTask;
use CQS::NGSCommon;
use Data::Dumper;

our @ISA = qw(CQS::PairTask);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_cdiff";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my $transcript_gtf = parse_param_file( $config, $section, "transcript_gtf", 1 );

  my $faFile = get_param_file( $config->{$section}{fasta_file}, "fasta_file", 1 );

  my $raw_files = get_raw_files( $config, $section );

  #print Dumper($raw_files);

  my $groups = get_raw_files( $config, $section, "groups" );

  #print Dumper($groups);

  my $pairs = get_raw_files( $config, $section, "pairs" );

  #print Dumper($pairs);

  my $mapfile = $result_dir . "/${task_name}_group_sample.map";
  open( MAP, ">$mapfile" ) or die "Cannot create $mapfile";
  print MAP "GROUP_INDEX\tSAMPLE_NAME\tGROUP_SAMPLE\tGROUP\tIndex\n";

  my %tpgroups         = ();
  my %group_sample_map = ();
  for my $group_name ( sort keys %{$groups} ) {
    my @samples = @{ $groups->{$group_name} };
    my @gfiles  = ();
    my $index   = 0;
    foreach my $sample_name ( sort @samples ) {
      my @bam_files = @{ $raw_files->{$sample_name} };
      push( @gfiles, $bam_files[0] );
      my $group_index = $group_name . "_" . $index;
      print MAP $group_index . "\t" . $sample_name . "\t" . $group_name . "_" . $sample_name . "\t" . $group_name . "\t" . $index . "\n";
      $group_sample_map{$group_index} = $sample_name;
      $index = $index + 1;
    }
    $tpgroups{$group_name} = join( ",", @gfiles );
  }
  close(MAP);

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct);

  for my $pair_name ( sort keys %{$pairs} ) {
    my ( $ispaired, $gNames ) = get_pair_groups( $pairs, $pair_name );

    #print Dumper($gNames);

    my @group_names = @{$gNames};
    my @bams        = ();
    foreach my $group_name (@group_names) {
      push( @bams, $tpgroups{$group_name} );
    }
    my $bamstrs = join( " ", @bams );

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $pair_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $pair_name );

    my $cur_dir = create_directory_or_die( $result_dir . "/$pair_name" );

    my $labels = join( ",", @group_names );

    my $log_desc = $cluster->get_log_description($log);

    my $final_file = "gene_exp.diff";

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir, $final_file );
    print $pbs "cuffdiff $option -o . -L $labels -b $faFile $transcript_gtf $bamstrs";
    $self->close_pbs( $pbs, $pbs_file );

    print $sh "\$MYCMD ./$pbs_name \n";
  }

  print $sh "exit 0\n";
  close $sh;

  my $sigfile = $pbs_dir . "/${task_name}_sig.pl";
  open( $sh, ">$sigfile" ) or die "Cannot create $sigfile";

  print $sh "#!/usr/bin/perl
use strict;
use warnings;

use CQS::RNASeq;

my \$config = {
  rename_diff => {
    target_dir => \"${target_dir}/result/comparison\",
    root_dir   => \"${target_dir}/result\",
    gene_only  => 0
  },
};

copy_and_rename_cuffdiff_file(\$config, \"rename_diff\");

1;

";
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }
  print "!!!shell file $shfile created, you can run this shell file to submit tasks.\n";
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $pairs = get_raw_files( $config, $section, "pairs" );

  my $result = {};
  for my $pair_name ( sort keys %{$pairs} ) {
    my $cur_dir      = $result_dir . "/$pair_name";
    my @result_files = ();
    push( @result_files, $cur_dir . "/gene_exp.diff" );
    push( @result_files, $cur_dir . "/genes.read_group_tracking" );
    push( @result_files, $cur_dir . "/splicing.diff" );
    push( @result_files, $result_dir . "/${task_name}_group_sample.map" );

    $result->{$pair_name} = filter_array( \@result_files, $pattern );
  }
  return $result;
}

1;
