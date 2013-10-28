package CQS::HashUtils;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [ qw( count_unique 
                                    return_unique
                                    return_unique_count) ] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub count_unique {
    my @array = @_;
    my %count;
    map { $count{$_}++ } @array;

    #print them out:
    #map {print "$_ = ${count{$_}}\n"} sort keys(%count);
    #or just return the hash:
    
    return %count;
}

sub return_unique {
    my @array = @_;
    my %count;
    map {$count{$_} = 1} @array;
    return sort keys(%count);
}

sub return_unique_count {
    my @array = @_;
    my %count;
    map {$count{$_} = 1} @array;
    
    my $result = keys(%count);
    return $result;
}