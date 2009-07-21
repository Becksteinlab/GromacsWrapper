#!/usr/bin/perl

#
# convert_to_gnuplot.pl - convert the single-column matrix data output from GridMAT-MD
# into a format that can be read by Gnuplot
#
# For use with GridMAT-MD, version 1.0.2 and later
#
# GridMAT-MD was developed by W. J. Allen and J. A. Lemkul in the lab of D. R. Bevan
# If you use GridMAT-MD, please cite the paper quoted in the output of that script
#
# Visit the GridMAT-MD website at: http://www.bevanlab.biochem.vt.edu/GridMAT-MD
#

use strict;

# open the input from the command line

unless (@ARGV) {
	die "Usage: $0 <input.dat> <x-grid points> <y-grid points>\n";
}

open(IN, $ARGV[0]) or die "Cannot open input: $!\n";
chomp(my @input = <IN>);
close(IN);

my $name = $ARGV[0];
my @name_array = split('\.', $name);
my $short_name = $name_array[0];

# get grid point information
my $x_grid = $ARGV[1];
my $y_grid = $ARGV[2];

# Data format: x y Z-value
# i.e., 1 1 3.295; 2 1 3.320; etc.
# nested counters necessary

# strip out unnecessary header line(s) from @input
my @input_clean;

foreach $_ (@input) {
	if ($_ =~ /^#/) {
		# do nothing
	} else {
		push(@input_clean, $_);
	}
}

# gnuplot matrix format:
# z11 z12 z13 ...
# z21 z22 z23 ...
# ...
# (each Z-value for an (x,y) pair

my $out_name = $short_name."_gnuplot.dat";

print "\nOutput file: $out_name\n\n";

open(OUT_RAW, ">>raw_gnuplot.dat") or die "Cannot open raw output: $!\n";
print OUT_RAW "# Matrix data: Z-values for each (x,y) pair\n";

for (my $y = 1; $y < $y_grid; $y++) {

	for (my $x = 1; $x < $x_grid; $x++) {
		
		for (my $z_value = 1; $z_value <= $x_grid; $z_value++) {
			
			# remove each line sequentially = don't worry about shifting indices?
			my $z = shift(@input_clean);
			printf OUT_RAW "%.3f\t", $z;
		
		}

		# newline at end of each x-row
		printf OUT_RAW "\n";

	}

}

close(OUT_RAW);

# Something weird is going on with all those loops...printing rows of zeroes at the end
# This next part is a bit of a hack to clean up the output

open(CLEAN_UP, "<raw_gnuplot.dat") or die "Cannot clean up output: $!\n";
my @clean = <CLEAN_UP>;
my @final;

foreach $_ (@clean) {

	if ($_ =~ /0.000/) {
		# do nothing
	} else {
		push(@final, $_);
	}
}

# write final output

open(OUT_FINAL, ">>${out_name}") or die "Cannot open output: $!\n";
foreach $_ (@final) {
	print OUT_FINAL $_;
}
close(OUT_FINAL);

unlink "raw_gnuplot.dat";

exit;
