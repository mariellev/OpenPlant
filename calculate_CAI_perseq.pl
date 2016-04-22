#!/usr/bin/perl

# script to calculate the Codon Adaptation Index for each sequence in a dataset
# use the relative adaptiveness of each codon calculated in calculate_codon_usage.pl
# output is a list with 2 columns: column 1: sequence ID, column2: CAI value for that sequence

use strict;
use warnings;
use Bio::SeqIO;

if ($ARGV[0] eq '-h' || $ARGV[0] eq '-help') {
    &help();
    exit;
}

# input file of CDS sequences
my $fasta_file = shift;
# file of relative adaptiveness (w) calculated in previous step
# reference file for whole transcriptome can be typically used
my $rel_adaptiveness_values = shift;

my $in = Bio::SeqIO->new(
    -file   => "$fasta_file",
    -format => 'fasta'
);

my %CDS;
while ( my $obj = $in->next_seq ) {
    my $CDSid    = $obj->id;
    my $sequence = $obj->seq;
    $sequence =~ tr/U/T/;
    $sequence =~ tr/u/t/;
    $sequence = uc($sequence);
    $CDS{$CDSid} = $sequence;
}

# extract relative adaptiveness values
my %rel_adapt;
open( WJ, $rel_adaptiveness_values ) or die 'cannot open $rel_adaptiveness $!';

while (<WJ>) {
    my $line = $_;
    chomp $line;
    ( my $aa, my $codon, my $wj ) = split( / : /, $line );
    $rel_adapt{$codon} = $wj;
}
close WJ;

# checking it works
foreach my $accession ( sort keys %CDS ) {
    my $length = length( $CDS{$accession} );
    if ( $length % 3 != 0 ) {
        print "Nucleotides not a multiple of 3\n";
    }
    my $j;
    my $log_wj        = 0;
    my $sum_log_wj    = 0;
    my $number_codons = 0;
    for ( $j = 0 ; $j < $length ; $j = $j + 3 ) {
        my $codon = substr( $CDS{$accession}, $j, 3 );
        if ( exists $rel_adapt{$codon} ) {
            $log_wj     = log( $rel_adapt{$codon} );
            $sum_log_wj = $log_wj + $sum_log_wj;
            $number_codons++;
        }
    }
    my $CAI = exp( ( 1 / $number_codons ) * $sum_log_wj );
    #print "$accession\t$number_codons\t$sum_log_wj\t$CAI\n";
    print "$accession\t$CAI\n";
    if ( $CAI == 1 ) {
        print "PERFECT SEQUENCE: $accession\n";
    }
}

sub help {
    print "USAGE:\n";
    print "./calculate_CAI_perseq.pl [fasta_filename] [file_of_wj_values] > [outputfile]\n";
}



