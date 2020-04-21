#!/usr/bin/perl

use Bio::SeqIO;

my $qfile = $ARGV[0] or die "This script will remove entries
from a specified fasta file (arg1) that are shorter than minimum length (arg2)\n";
my $minlen=$ARGV[1] or die "This script will remove entries annotated containing words
from a specified fasta file (arg1) that are shorter than minimum length (arg2)\n";

my $in= Bio::SeqIO->new(-file => $ARGV[0],
                           -format => 'Fasta');
                           
my $out=Bio::SeqIO->new(-file => ">noshorts_".$ARGV[0],
                           -format => 'Fasta',
                           -alphabet => 'dna');
my $shit=0;
my $good=0;
while ( my $seq = $in->next_seq ) {
	$len=$seq->length;
	if ($len>=$minlen) {
		$out->write_seq($seq);
		$good++;
	}
	else { $shit++; }
}
print "noshorts:\nretained:\t$good\ndiscarded:\t$shit\n";