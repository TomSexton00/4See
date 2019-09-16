#!/usr/bin/perl

#############################################################
#Utility script to convert mapped sequences from coordinate space to restriction fragment space
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

die "Usage: $0 <frag_fn> <map_fn> <chr_col> <coord_col> <strand_col> <cis_fn> <chrom> <vp.pos> <name> <read length>" if ($#ARGV != 9);
my $frag_fn = $ARGV[0];
my $read_fn = $ARGV[1];
my $chr_col = $ARGV[2]-1;
my $coord_col = $ARGV[3]-1;
my $strand_col = $ARGV[4]-1;
my $cis_file = $ARGV[5];
my $cis_chrom = $ARGV[6];
my $vp = $ARGV[7];
my $pre = $ARGV[8];
my $rlen = $ARGV[9];

#Data structures
my %frags;
#Lookup for quick searching between coord and frag
my %coord2index;

#Reads mapping to chromosome extremities
my $discard = 0;

##############################################################################
# Reads fends file
##############################################################################

open IN,$frag_fn or die $!;
print "Reading frags file $frag_fn...\n";
my $header = <IN>;
my %h = parse_header($header);
my $frag_count = 0;
while (my $line = <IN>) {
	++$frag_count;
	print "$frag_count\n" if ($frag_count % 1000000 == 0);
	chomp $line;
	my @f = split(/\t/,$line);
	my $chr = $f[$h{CHROM}];
	my $frag = $f[$h{FRAG_ID}];
	my $coord = $f[$h{COORD}];
	my $frag_len = $f[$h{FRAG_LEN}];
	my $start = $coord;
	my $end = $coord+$frag_len;
	my $mid = $coord+int($frag_len/2);
	$frags{$frag}={start=>$start,mid=>$mid,end=>$end,count=>0};
	$coord2index{$chr}->{$coord}=$frag;
}
close IN or die $!;
print "$frag_count frags read in\n";

foreach my $chr(keys %coord2index) {
	my @sorted = sort {$a <=> $b} keys %{$coord2index{$chr}};
	$coord2index{$chr}->{sorted}=\@sorted;
}

##############################################################################
# Parse reads file
##############################################################################

my $read_count = 0;
my $kept;
open IN,$read_fn or die $!;
print "Reading in $read_fn...\n";
LOOP: while (my $line = <IN>) {
	++$read_count;
	print "$read_count\n" if ($read_count % 1000000 == 0);
	chomp $line;
	my @f = split(/\t/,$line);
	my $chr = $f[$chr_col];
	my $coord = $f[$coord_col];
	my $strand = $f[$strand_col];
	
	next if (!defined($coord2index{$chr}));

	#Trim coord to first nt after expected cutter site, expecting 0-base coordinates
	if ($strand eq "+" or $strand eq "1") {
		$coord+=4;
	}
	elsif ($strand eq "-" or $strand eq "-1") {
		$coord+=$rlen-5;
	}
	else {
		die "Unrecognised strand in read file\n";
	}

	my $frag = coord_to_frag($chr,$coord);
	if ($frag==-1) {
		$frag = coord_to_frag($chr,$coord-1);
		if ($frag==-1) {
			$frag = coord_to_frag($chr,$coord+1);
			if ($frag==-1) {
				++$discard;
				next LOOP;
			}
		}
	}

	++$kept;
	my $score = $frags{$frag}->{count};
	++$score;
	$frags{$frag}->{count}=$score;
}
close IN or die $!;

# Define 2% threshold for self-lig/PCR-dominant reads
my $threshold = $kept * 0.02;

#############################################################################
# Output
#############################################################################

open CIS,">",$cis_file or die $!;
print CIS $pre,"\tchr",$cis_chrom,"\t",$vp,"\n";
print "Writing $cis_file...\n";

my @chrs = sort {$a cmp $b} keys %coord2index;

foreach my $chr(@chrs) {
	foreach my $fragcoord(@{$coord2index{$chr}->{sorted}}) {
		my $frag = $coord2index{$chr}->{$fragcoord};
		if ($frags{$frag}->{count}>$threshold) {
			print "Ignored frag: ",$chr,": ",$frags{$frag}->{mid}," - ",(($frags{$frag}->{count})/$kept)*100," %\n";
		}
		else {
			if ($chr eq $cis_chrom) {
				print CIS $frags{$frag}->{mid},"\t",$frags{$frag}->{count},"\n";
			}
		}
	}
}
close CIS or die $!;

print "$discard discarded reads out of $read_count\n";

##############################################################################
# Sub-routines
##############################################################################

 sub parse_header {
	my ($header) = @_;
	chomp ($header);
	my @f = split(/\t/,$header);
	my %result;
	for (my $i = 0; $i <= $#f; $i++) {
		$result{$f[$i]} = $i;
	}
	return %result;
}

sub binary_search {
	my $arr = shift;
	my $value = shift;
	my $left = 0;
	my $right = $#$arr;
	while ($left<=$right) {
		my $mid = ($right+$left) >> 1;
		my $c = $arr->[$mid] <=> $value;
		return $mid if ($c==0);
		if ($c>0) {
			$right = $mid - 1;
		}
		else {
			$left = $mid + 1;
		}
	}
	$left = -1 if ($left>$#$arr);
	$right = -1 if ($right < 0);
	return $left-1;
}

sub coord_to_frag {
	my ($chr,$coord) = (@_);
	my $index_p = binary_search($coord2index{$chr}->{sorted},$coord);
	return (-1) if ($index_p == -1);
	my $fcoord_p = $coord2index{$chr}->{sorted}[$index_p];
	defined ($coord2index{$chr}->{$fcoord_p}) and defined ($fcoord_p) or die $!;
	my $frag = $coord2index{$chr}->{$fcoord_p};
	return ($frag);
}
