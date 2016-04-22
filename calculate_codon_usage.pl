#!/usr/bin/perl/

# script to calculate the codon usage for set of a sequences
# concatenate multiple CDS of each sequence to come up with
# count for each codon
# highlight preferred codon for each amino acid
# calculate relative adaptiveness wj for each codon of kind j
# requires genetic_code.txt

use strict;
use warnings;
use Bio::SeqIO;

if ($ARGV[0] eq '-h' || $ARGV[0] eq '-help') {
    &help();
    exit;
}

#input  file of CDS sequences
my $fasta_file       = shift;
# name of output file
my $rel_adaptiveness = shift;
my $gen_code = 'genetic_code.txt';
my @allCDS       = "";
my $allSequences = "";
my @counts;
# subroutines variable
my %nucleot = ();
my @codons;

open( GCODE, $gen_code ) or die "cannot open $gen_code $!";
my @gen_code = <GCODE>;
close GCODE;

my $in = Bio::SeqIO->new(
    -file   => "$fasta_file",
    -format => 'fasta'
);

while ( my $obj = $in->next_seq ) {
    my $seqid    = $obj->id;
    my $sequence = $obj->seq;
    $sequence =~ s/U/T/;
    $sequence = uc($sequence);
    push @allCDS, $sequence;
    my $length = length($sequence);
    if ( $length % 3 != 0 ) {
        print "$seqid: Nucleotides not a multiple of 3\n";
    }
}

# go through retrieved CDS and count the codons
foreach my $seq (@allCDS) {
    # concatenate all sequences 
    $allSequences = "$allSequences" . "$seq";
}

# count the codon in frame 0 (in-frame)
my %nucleotide_f0 = &count_codons( 0, $allSequences, \@gen_code );

# checking it works
foreach my $aa ( sort keys %nucleotide_f0 ) {
    for my $cod ( sort keys %{ $nucleotide_f0{$aa} } ) {
        print "$aa,$cod,$nucleotide_f0{$aa}{$cod}\n";
        push @counts, $nucleotide_f0{$aa}{$cod};
    }
}

my %Xmax_f0 = &retrieve_most_frequent( \%nucleotide_f0 );

# calculate relative adaptiveness wj for codon of kind j
open( WJ, ">$rel_adaptiveness" ) or die "cannot create $rel_adaptiveness $!\n";
foreach my $aa ( sort keys %nucleotide_f0 ) {
    for my $cod ( sort keys %{ $nucleotide_f0{$aa} } ) {
        # print "$cod : $aa: $nucleotide_f0{$aa}{$cod}*($nucleotide_f0{$aa}{$cod}/$Xmax_f0{$aa})\n";
        my $w = $nucleotide_f0{$aa}{$cod} / $Xmax_f0{$aa};
        print WJ "$aa : $cod : $w\n";
    }
}

# retrieve each codon and keep count of each of their occurrence
sub count_codons {
    my ( $j, $seq, $refGen_code ) = @_;
    my $starting_point = $j;
    my $k;
    @codons = ();
    if ( $starting_point == 0 ) {
        $k = 0;
    }
    elsif ( $starting_point == 1 ) {
        $k = 2;
    }
    elsif ( $starting_point == 2 ) {
        $k = 1;
    }
    # Get the array from the reference
    my @genetic_code = @{$refGen_code};
    # reset the hashtable
    %nucleot = ();
    for ( $j ; $j < ( length($seq) - $k ) ; $j = $j + 3 ) {
        my $codon = substr( $seq, $j, 3 );
        for ( my $i = 0 ; $i < scalar @genetic_code ; $i++ ) {
            my @fields = split( / /, $genetic_code[$i] );
            my $aa     = $fields[1];
            my $cod    = $fields[0];
            if ( $codon eq $cod ) {
                $nucleot{$aa}{$codon}++;
                #print "$aa $codon\n";
            }
        }
    }
    return %nucleot;
}

# retrieve the most frequent codon for each amino acid
# from http://www.scriptscoop.net/t/9ca1612d592e/perl-extract-values-from-multi-dimensional-hash.html
sub retrieve_most_frequent {
    my $nucl     = shift;
    my %nuclhash = %$nucl;
    my %Xmax;
    for my $AminoAcid ( sort keys %nuclhash ) {
        my $key = (
            reverse
              sort { $nuclhash{$AminoAcid}{$a} <=> $nuclhash{$AminoAcid}{$b} }
              keys %{ $nuclhash{$AminoAcid} } )[0];
        print "$AminoAcid\t$key\t$nuclhash{$AminoAcid}{$key}\n";
        $Xmax{$AminoAcid} = $nuclhash{$AminoAcid}{$key};
    }
    return %Xmax;
}

sub help {
    print "USAGE:\n";
    print "./calculate_codon_usage [fasta_filename] [output_filename]\n";
    print "genetic_code.txt file needs to be placed in the current directory\n";
}
