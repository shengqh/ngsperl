use strict;
use File::Basename;
use LWP::Simple;
use LWP::UserAgent;
use Bio::SeqIO;
use URI::Escape;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw(getStructure readIntronFromBed buildMatureFasta)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

sub getStructure {
  my $structureFile = shift;

  if ( !-e $structureFile ) {
    my $tmpFile = $structureFile . ".tmp";

    my $ua = new LWP::UserAgent;
    $ua->agent( "AgentName/0.1 " . $ua->agent );

    # Create a request
    my $url = 'http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/';
    my $req = new HTTP::Request GET => $url;

    # Pass request to the user agent and get a response back
    my $res = $ua->request($req);

    if ( $res->is_success ) {
      open( my $tmp, ">$tmpFile" ) or die "Cannot create $tmpFile";
      my $rescontent = $res->content;

      my @categories = ( $rescontent =~ m/folder.gif" alt="\[DIR\]"> <a href="(.*?)\/"/g );
      foreach my $category (@categories) {
        print $category, "\n";

        my $categoryurl     = $url . $category;
        my $categoryreq     = new HTTP::Request GET => $categoryurl;
        my $categoryres     = $ua->request($categoryreq);
        my $categorycontent = $categoryres->content;

        my @species_array = $categorycontent =~ m/folder.gif" alt="\[DIR\]"> <a href="(.*?)\/"/g;

        foreach my $species (@species_array) {
          if ( $species =~ /_old/ ) {
            next;
          }
          my $speciesurl     = $categoryurl . "/" . $species;
          my $speciesreq     = new HTTP::Request GET => $speciesurl;
          my $speciesres     = $ua->request($speciesreq);
          my $speciescontent = $speciesres->content;

          #print $speciescontent, "\n";

          my $tarUrl;
          my $faUrl;
          my $filePrefix;
          if ( $speciescontent =~ /href="(.*?)\.tar.gz">/ ) {
            $filePrefix = $1;
            my $tarfile = $1 . ".tar.gz";
            $tarUrl = $speciesurl . "/" . uri_escape($tarfile);
          }
          else {
            print "Cannot find tar url of " . $species . " in " . $category . "; Species url = " . $speciesurl . "\n";
            next;
          }

          if ( $speciescontent =~ /href="(.*?\.fa)">/ ) {
            my $fastafile = $1;
            $faUrl = $speciesurl . "/" . uri_escape($fastafile);
          }
          else {
            print "Cannot find fa url of " . $species . " in " . $category . "; Species url = " . $speciesurl . "\n";
            next;
          }

          print $tmp $species . "\t" . $category . "\t" . $filePrefix . "\t" . $tarUrl . "\t" . $faUrl . "\n";
        }
      }

      if ( -s $structureFile ) {
        unlink($structureFile);
      }
      rename( $tmpFile, $structureFile );
    }
  }

  my $result = [];
  open( my $sr, $structureFile ) or die "Could not open file '$structureFile' $!";
  while ( my $row = <$sr> ) {
    chomp $row;
    my @parts = split( "\t", $row );
    push( @$result, \@parts );
  }

  close($sr);

  return $result;
}

sub readIntronFromBed {
  my $bedFile = shift;

  open( my $bedio, $bedFile ) or die "Could not open file '$bedFile' $!";
  my $result = {};
  while ( my $row = <$bedio> ) {
    chomp $row;
    my @parts     = split( "\t", $row );
    my $name      = $parts[3];
    my @sizes     = map( int($_), split( ",", $parts[10] ) );
    my @positions = map( int($_), split( ",", $parts[11] ) );

    $result->{$name} = { sizes => \@sizes, positions => \@positions };
  }

  return $result;
}

sub buildMatureFasta {
  my ( $bedFile, $fastaFile, $matureFile ) = @_;

  my $beds = readIntronFromBed($bedFile);

  my $seqio = Bio::SeqIO->new( -file => $fastaFile, -format => 'fasta' );

  my $tmpFile = $matureFile . ".tmp";
  open( my $output, ">$tmpFile" ) or die "Cannot create file $tmpFile";
  while ( my $seq = $seqio->next_seq ) {
    my $id       = $seq->id;
    my $sequence = $seq->seq;
    my $desc     = $seq->desc;

    my $values = undef;
    for my $key ( keys %$beds ) {
      my $keylen = length($key);
      if ( $key eq substr( $id, -$keylen ) ) {
        $values = $beds->{$key};
        last;
      }
    }

    $id = "$id $desc";

    if ( !defined $values ) {
      print STDERR "Cannot find " . $id . " in bed file of " . $bedFile . "\n";
      print $output ">$id\n$sequence\n";
      next;
    }

    my $sizes     = $values->{sizes};
    my $positions = $values->{positions};

    my $numberOfExon = scalar(@$sizes);
    if ( $numberOfExon > 1 ) {
      my $tmpSeq = "";
      for my $i ( 0 .. ( $numberOfExon - 1 ) ) {
        my $start  = @$positions[$i];
        my $length = @$sizes[$i];
        $tmpSeq = $tmpSeq . substr( $sequence, $start, $length );
      }
      $sequence = $tmpSeq;
      $id       = $id . " intron_removed";
    }
    print $output ">$id\n$sequence\n";
  }
  close($output);

  if ( -s $matureFile ) {
    unlink($matureFile);
  }
  rename( $tmpFile, $matureFile );
}
