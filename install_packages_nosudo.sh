#!/bin/bash

host=$(echo $HOSTNAME | cut -d'.' -f1)
lib_path="/nobackup/h_cqs/perl5/perl5_lib_${host}"
#echo $lib_path

if( [ ! -d $lib_path ] ) 
then
  mkdir -p $lib_path
fi

echo "Install perl modules to $lib_path"

curl -L http://cpanmin.us | perl - -l $lib_path Data::Dumper
curl -L http://cpanmin.us | perl - -l $lib_path Data::Table
curl -L http://cpanmin.us | perl - -l $lib_path File::Basename
curl -L http://cpanmin.us | perl - -l $lib_path File::Copy
curl -L http://cpanmin.us | perl - -l $lib_path File::Slurp
curl -L http://cpanmin.us | perl - -l $lib_path File::Spec
curl -L http://cpanmin.us | perl - -l $lib_path File::Temp
curl -L http://cpanmin.us | perl - -l $lib_path Getopt::Long
curl -L http://cpanmin.us | perl - -l $lib_path Hash::Merge
curl -L http://cpanmin.us | perl - -l $lib_path HTML::Entities
curl -L http://cpanmin.us | perl - -l $lib_path HTML::Parser
curl -L http://cpanmin.us | perl - -l $lib_path IO::Uncompress::Gunzip
curl -L http://cpanmin.us | perl - -l $lib_path JSON
curl -L http://cpanmin.us | perl - -l $lib_path List::Compare
curl -L http://cpanmin.us | perl - -l $lib_path List::Util
curl -L http://cpanmin.us | perl - -l $lib_path LWP::Simple
curl -L http://cpanmin.us | perl - -l $lib_path LWP::UserAgent
curl -L http://cpanmin.us | perl - -l $lib_path Net::FTP
curl -L http://cpanmin.us | perl - -l $lib_path Scalar::Util
curl -L http://cpanmin.us | perl - -l $lib_path Statistics::Descriptive
curl -L http://cpanmin.us | perl - -l $lib_path String::Util
curl -L http://cpanmin.us | perl - -l $lib_path Test::Deep
curl -L http://cpanmin.us | perl - -l $lib_path Test::Simple
curl -L http://cpanmin.us | perl - -l $lib_path WWW::Mechanize
curl -L http://cpanmin.us | perl - -l $lib_path Algorithm::Loops
curl -L http://cpanmin.us | perl - -l $lib_path JSON::Streaming::Reader
curl -L http://cpanmin.us | perl - -l $lib_path Tie::IxHash
curl -L http://cpanmin.us | perl - -l $lib_path Text::CSV

