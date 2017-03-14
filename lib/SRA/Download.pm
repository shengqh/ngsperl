#!/usr/bin/perl
package SRA::Download;

use strict;
use warnings;
use File::Basename;
use CQS::PBS;
use CQS::ConfigUtils;
use CQS::SystemUtils;
use CQS::FileUtils;
use CQS::Task;
use CQS::StringUtils;
use Utils::CollectionUtils;
use LWP::Simple;
use LWP::UserAgent;
use URI::Escape;

our @ISA = qw(CQS::Task);

sub new {
  my ($class) = @_;
  my $self = $class->SUPER::new();
  $self->{_name}   = __PACKAGE__;
  $self->{_suffix} = "_ds";
  bless $self, $class;
  return $self;
}

sub perform {
  my ( $self, $config, $section ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct, $cluster ) = get_parameter( $config, $section );

  my %raw_files = %{ get_raw_files( $config, $section ) };

  my $shfile = $self->get_task_filename( $pbs_dir, $task_name );
  open( my $sh, ">$shfile" ) or die "Cannot create $shfile";
  print $sh get_run_command($sh_direct) . "\n";

  my $ua = new LWP::UserAgent;
  $ua->agent( "AgentName/0.1 " . $ua->agent );

  for my $sample_name ( sort keys %raw_files ) {
    my @sample_files = @{ $raw_files{$sample_name} };
    my $list_file    = $sample_files[0];

    my $pbs_file = $self->get_pbs_filename( $pbs_dir, $sample_name );
    my $pbs_name = basename($pbs_file);
    my $log      = $self->get_log_filename( $log_dir, $sample_name );
    my $log_desc = $cluster->get_log_description($log);

    my $pbs = $self->open_pbs( $pbs_file, $pbs_desc, $log_desc, $path_file, $result_dir );

    print $sh "\$MYCMD ./$pbs_name \n";
    my $taMap = readDictionaryByIndex( $list_file, 19, 16, 1 );

    for my $gsm ( sort keys %$taMap ) {
      my $srr = $taMap->{$gsm};
      my $url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" . $gsm;
      my $req = HTTP::Request->new(GET => $url);

      # Pass request to the user agent and get a response back
      my $res = $ua->request($req);

      if ( $res->is_success ) {
        my $rescontent = $res->content;
        #print $rescontent;
        
        my @sraUrls = ( $rescontent =~ m/<a href=\"(ftp:\/\/ftp\-trace.ncbi.nlm.nih.gov\/sra.*?)\">.ftp/g );
        foreach my $sraUrl (@sraUrls) {
          my $fileUrl = $sraUrl . "/" . $srr . "/" .$srr . ".sra";
          print $pbs "if [ ! -s ${srr}.sra ]; then
  wget $fileUrl
fi

";
        }
      }
    }

    $self->close_pbs( $pbs, $pbs_file );
  }
  close $sh;

  if ( is_linux() ) {
    chmod 0755, $shfile;
  }

  print "!!!shell file $shfile created, you can run this shell file to submit all SRA::Download tasks.\n";

  #`qsub $pbs_file`;
}

sub result {
  my ( $self, $config, $section, $pattern ) = @_;

  my ( $task_name, $path_file, $pbs_desc, $target_dir, $log_dir, $pbs_dir, $result_dir, $option, $sh_direct ) = get_parameter( $config, $section, 0 );

  my $result = {};

  return $result;
}

1;
