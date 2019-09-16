#!/usr/bin/perl

#############################################################
#Utility script to generate restriction fragment maps from fasta sequences
#Copyright (C) (2019) Tom Sexton & Yousra Ben Zouari
#
#This program is free software; you can distribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3
#
#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
##############################################################


use strict;
use warnings;

my(%rc) = ("A", "T", "C", "G", "G", "C", "T", "A", "N", "N");

# returns reverse-complement of the argument
sub get_rc {
	my($str, $len) = @_;
	my($rcstr);
	for (my($i) = $len - 1; $i >= 0; $i--) {
		$rcstr .= $rc{substr($str, $i, 1)};
	}
	return($rcstr);
}

die "Usage: $0 FIRST_CUTTER TRACKDB_SEQ_DIR FRAG_LEN OUT_FILE\n" if $#ARGV < 3;

my($site) = $ARGV[0];
my($site_len) = length($site);
my($frag_len) = $ARGV[2];
my($outfile) = $ARGV[3];

my(@seq_fns) = <$ARGV[1]/*.fa>;

my(@hits);

# verify that the site is reverse-symmetric (i.e. the reverse-complement of the site is equal to the site)
$site eq get_rc($site, $site_len) || die "Site must be reverse-symmetric";

print "Extracting $site\n";

my($frag_id) = 0;

my($ser_id) = 0;
my($seq_fn);
foreach $seq_fn (@seq_fns) {
	open(SEQ, $seq_fn) || die "cannot open fasta $seq_fn\n";

	print "will read $seq_fn\n";
	my($chrom) = $seq_fn=~/chr(.+).fa/;

	my $seq = "";
	while (my $line = <SEQ>) {
		chomp $line;
		next unless ($line);
		$seq .= uc($line);
	}
	close(SEQ);

	print "done reading $seq_fn, got ".length($seq)." bps\n";
	my($max_i) = length($seq)-$site_len;
	my($start_coord) = index($seq, $site);
	my($mid_coord) = index($seq, $site, $start_coord + $site_len);
	my($end_coord) = index($seq, $site, $mid_coord + $site_len);

	while ($mid_coord != -1 && $end_coord != -1) {
		# strand +
		my($frag) = substr($seq, $mid_coord + $site_len, $frag_len);
		if(length($frag) == $frag_len) {
			my($hit) = ($frag_id+1)."\tchr$chrom\t".($mid_coord-$site_len-1)."\t";
			$hit .= "\t" . ($end_coord - $mid_coord);  # fragment length
			push(@hits, $hit);
			$ser_id++;
		}

		$start_coord = $mid_coord;
		$mid_coord = $end_coord;
		$end_coord = index($seq, $site, $mid_coord + $site_len);
		$frag_id++;
	}
	$frag_id++;
	
	print "counted $ser_id at chrom $chrom\n";
}

open(TAB, ">$outfile");

# write the header
print TAB "FRAG_ID\tCHROM\tCOORD\tFRAG_LEN\n";

# write hits
for(my($i) = 0; $i <= $#hits; $i++) {
	my(@hit) = split("\t", $hits[$i]);
	print TAB "$hit[0]";
	for (my($i) = 1; $i <= $#hit; $i++) {
		print TAB "\t$hit[$i]";
	}
	print TAB "\n";
}

close(TAB);
