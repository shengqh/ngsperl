package CQS::SystemUtils;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw( get_run_now is_linux)] );

our @EXPORT = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';

sub get_run_now {
	my $result = "";
	if ( $#ARGV >= 0 ) {
		$result = $ARGV[0] eq "y";
	}
	return ($result);
}

sub is_linux {
    my $os = $^O;
		#print("os=". $os . "\n");
		my $result = $os eq "linux" ? 1 : 0;
		#print("result=". $result . "\n");
    return ( $result );
}

1;
