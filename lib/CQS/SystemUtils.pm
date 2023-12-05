package CQS::SystemUtils;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [qw( get_run_now is_linux add_copy_data)] );

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

sub add_copy_data {
  my ( $config, $def, $source_ref ) = @_;
  
  $config->{"copy_data"} = {
    class => 'CQS::ProgramWrapperOneToOne',
    target_dir => $def->{target_dir} . '/copy_data',
    option => "

rm -f __NAME__.done __NAME__.failed

cp __FILE__ .

status=\$?
if [[ \$status -eq 0 ]]; then
  touch __NAME__.done
else
  touch __NAME__.failed
fi

# __OUTPUT__

",
    interpretor => '',
    program => '',
    check_program => 0,
    source_arg => '',
    source_ref => $source_ref,
    source_join_delimiter => ' ',
    output_by_file => 1,
    output_ext => '.done',
    output_to_same_folder => 0,
    use_tmp_folder => 0,
    pbs => {
      'email' => $def->{email},
      'emailType' => $def->{emailType},
      'nodes' => '1:ppn=1',
      'walltime' => '24',
      'mem' => '10gb'
    },
  };
  return($config);
}


1;
