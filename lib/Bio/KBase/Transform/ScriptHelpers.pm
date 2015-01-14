use strict;
use warnings;
package Bio::KBase::Transform::ScriptHelpers;
use Data::Dumper;
use Config::Simple;
use Getopt::Long::Descriptive;
use Text::Table;
use JSON::XS;
use Bio::KBase::Auth;
use Exporter;
use File::Stream;
use File::Slurp;
use parent qw(Exporter);
our @EXPORT_OK = qw( parse_input_table get_input_fh get_output_fh load_input write_output write_text_output );

sub parse_input_table {
	my $filename = shift;
	my $columns = shift;#[name,required?(0/1),default,delimiter]
	if (!-e $filename) {
		print "Could not find input file:".$filename."!\n";
		exit();
	}
	if($filename !~ /\.([ct]sv|txt)$/){
    	die("$filename does not have correct suffix (.txt or .csv or .tsv)");
	}
	open(my $fh, "<", $filename) || return;
	my $headingline = <$fh>;
	chomp($headingline);
	my $delim = undef;
	if ($headingline =~ m/\t/) {
		$delim = "\\t";
	} elsif ($headingline =~ m/,/) {
		$delim = ",";
	}
	if (!defined($delim)) {
		die("$filename either does not use commas or tabs as a separator!");
	}
	my $headings = [split(/$delim/,$headingline)];
	my $data = [];
	while (my $line = <$fh>) {
		chomp($line);
		push(@{$data},[split(/$delim/,$line)]);
	}
	close($fh);
	my $headingColums;
	for (my $i=0;$i < @{$headings}; $i++) {
		$headingColums->{$headings->[$i]} = $i;
	}
	my $error = 0;
	for (my $j=0;$j < @{$columns}; $j++) {
		if (!defined($headingColums->{$columns->[$j]->[0]}) && defined($columns->[$j]->[1]) && $columns->[$j]->[1] == 1) {
			$error = 1;
			print "Model file missing required column '".$columns->[$j]->[0]."'!\n";
		}
	}
	if ($error == 1) {
		exit();
	}
	my $objects = [];
	foreach my $item (@{$data}) {
		my $object = [];
		for (my $j=0;$j < @{$columns}; $j++) {
			$object->[$j] = undef;
			if (defined($columns->[$j]->[2])) {
				$object->[$j] = $columns->[$j]->[2];
			}
			if (defined($headingColums->{$columns->[$j]->[0]}) && defined($item->[$headingColums->{$columns->[$j]->[0]}])) {
				$object->[$j] = $item->[$headingColums->{$columns->[$j]->[0]}];
			}
			if (defined($columns->[$j]->[3])) {
				if (defined($object->[$j]) && length($object->[$j]) > 0) {
					my $d = $columns->[$j]->[3];
					$object->[$j] = [split(/$d/,$object->[$j])];
				} else {
					$object->[$j] = [];
				}
			}
		}
		push(@{$objects},$object);
	}
	return $objects;
}

sub get_input_fh
{
    my($opts) = @_;

    my $fh;
    if ($opts->{input})
    {
	open($fh, "<", $opts->{input}) or die "Cannot open input file $opts->{input} :$!";
    }
    else
    {
	$fh = \*STDIN;
     }
    return $fh;
}	

sub get_output_fh
{
    my($opts) = @_;

    my $fh;
    if ($opts->{output})
    {
	open($fh, ">", $opts->{output}) or die "Cannot open input file $opts->{output} :$!";
    }
    else
    {
	$fh = \*STDOUT;
    }
    return $fh;
}	

sub load_input
{
    my($opts) = @_;

    my $fh;
    if ($opts->{input})
    {
	open($fh, "<", $opts->{input}) or die "Cannot open input file $opts->{input} :$!";
    }
    else
    {
	$fh = \*STDIN;
    }

    my $text = read_file($fh);
    undef $fh;
    my $obj = decode_json($text);
    return $obj;
}

sub write_output
{
    my($genome, $opts) = @_;

    my $coder = JSON::XS->new->pretty;
    my $text = $coder->encode($genome);
    if ($opts->{output})
    {
	write_file($opts->{output}, \$text);
    }
    else
    {
	write_file(\*STDOUT, \$text);
    }
}

sub write_text_output
{
    my($text, $opts) = @_;

    if ($opts->{output})
    {
	write_file($opts->{output}, \$text);
    }
    else
    {
	write_file(\*STDOUT, \$text);
    }
}



1;
