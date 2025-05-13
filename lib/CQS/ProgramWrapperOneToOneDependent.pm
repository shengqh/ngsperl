#!/usr/bin/perl
package CQS::ProgramWrapperOneToOneDependent;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::NGSCommon;
use CQS::StringUtils;
use CQS::Task;
use File::Spec;
use Data::Dumper;
use CQS::ProgramWrapperOneToOne;

our @ISA = qw(CQS::ProgramWrapperOneToOne);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_o2o";
  $self->{_output_to_same_folder} = 1;
  $self->{max_jobs} = 10;
  bless $self, $class;
  return $self;
}

sub print_sh_pbs {
  my ($self, $sh, $pbs_name, $final_file, $pbs_index) = @_;

  my $pbs_index_key= "jid" . $pbs_index;
  my $max_jobs = $self->{max_jobs};

  if ($pbs_index <= $max_jobs) {
    print $sh "
if [[ ! -s $final_file ]]; then
  $pbs_index_key=\$(sbatch ./$pbs_name | awk '{print \$NF}')
else
  $pbs_index_key=1000000000
fi
";
  }else{
    my $dep_keys = [];
    my $dep_index = $pbs_index - $max_jobs;
    while($dep_index > 0) {
      my $dep_index_key = "jid" . $dep_index;
      push @$dep_keys, $dep_index_key;
      $dep_index -= $max_jobs;
    }
    my $dep_index_keys = join(':$', @$dep_keys);
    print $sh "
if [[ ! -s $final_file ]]; then
  $pbs_index_key=\$(sbatch --dependency=afterany:\$$dep_index_keys ./$pbs_name | awk '{print \$NF}')
else
  $pbs_index_key=1000000000
fi
";
  }
}

sub perform {
  my ( $self, $config, $section ) = @_;
  $self->{max_jobs} = get_option( $config, $section, "max_jobs", 10 );
  $self->CQS::ProgramWrapperOneToOne::perform( $config, $section );
}

1;
