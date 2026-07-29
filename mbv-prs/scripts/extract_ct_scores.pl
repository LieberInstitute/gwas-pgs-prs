#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);

my ($clumped, $lookup, $output);
GetOptions(
    'clumped=s' => \$clumped,
    'lookup=s'  => \$lookup,
    'output=s'  => \$output,
) or die usage();
die usage() unless defined $clumped && defined $lookup && defined $output;

my %lookup;
open my $lookup_fh, '<', $lookup or die "Cannot read $lookup: $!\n";
<$lookup_fh>;
while (my $line = <$lookup_fh>) {
    chomp $line;
    my ($ld_snp, @score) = split /\t/, $line, -1;
    die "Duplicate LD SNP in lookup: $ld_snp\n" if exists $lookup{$ld_snp};
    $lookup{$ld_snp} = \@score;
}
close $lookup_fh or die "Failed while reading $lookup\n";

open my $in, '<', $clumped or die "Cannot read $clumped: $!\n";
open my $out, '>', $output or die "Cannot write $output: $!\n";
print {$out} "SNP\tA1\tBETA\tP\n";

my ($snp_col, $selected) = (undef, 0);
while (my $line = <$in>) {
    next if $line =~ /^\s*$/;
    my @f = split /\s+/, trim($line);
    if (!defined $snp_col) {
        for my $i (0 .. $#f) {
            $snp_col = $i if $f[$i] eq 'SNP';
        }
        next if defined $snp_col;
        next;
    }

    my $ld_snp = $f[$snp_col];
    die "Clumped SNP missing from lookup: $ld_snp\n" unless exists $lookup{$ld_snp};
    my ($target_snp, $a1, $beta, $p) = @{ $lookup{$ld_snp} };
    print {$out} join("\t", $target_snp, $a1, $beta, $p), "\n";
    ++$selected;
}
close $in or die "Failed while reading $clumped\n";
close $out or die "Cannot close $output: $!\n";
print "$selected\n";

sub trim {
    my ($text) = @_;
    $text =~ s/^\s+//;
    $text =~ s/\s+$//;
    return $text;
}

sub usage {
    return "Usage: $0 --clumped FILE --lookup FILE --output FILE\n";
}
