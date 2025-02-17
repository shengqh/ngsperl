#!/bin/bash

if [[ $HOSTNAME == 'cqs1.vampire' || $HOSTNAME == 'cqs3.vampire' ]]; then
  export TARGET_FOLDER=/nobackup/h_vangard_1/shengq2/perl5_cqs13lib
elif [[ $HOSTNAME == 'cqs4.vampire' ]]; then
  export TARGET_FOLDER=/nobackup/h_cqs/softwares/perl5_lib_cqs4
else
  export TARGET_FOLDER="~/perl5"
fi

curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER Data::Dumper
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER Data::Table
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER File::Basename
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER File::Copy
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER File::Slurp
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER File::Spec
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER File::Temp
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER Getopt::Long
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER Hash::Merge
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER HTML::Entities
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER HTML::Parser
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER IO::Uncompress::Gunzip
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER JSON
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER List::Compare
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER List::Util
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER LWP::Simple
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER LWP::UserAgent
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER Net::FTP
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER Scalar::Util
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER Statistics::Descriptive
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER String::Util
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER Test::Deep
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER Test::Simple
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER WWW::Mechanize
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER Algorithm::Loops
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER JSON::Streaming::Reader
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER Tie::IxHash
curl -L http://cpanmin.us | perl - -l $TARGET_FOLDER Text::CSV

