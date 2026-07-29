#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use IO::Compress::Gzip qw($GzipError);
use IO::Uncompress::Gunzip qw($GunzipError);

my ($input, $output, $stats, $strict);
GetOptions(
    'input=s'  => \$input,
    'output=s' => \$output,
    'stats=s'  => \$stats,
    'strict'   => \$strict,
) or die "Usage: $0 --input FILE --output FILE.gz [--stats FILE]\n";

die "Missing --input\n" unless defined $input;
die "Missing --output\n" unless defined $output;
$stats //= "$output.validation.tsv";

my $in = open_input($input);
my $tmp = "$output.tmp.$$";
my $candidate = "$output.candidate.$$.gz";
my $out = IO::Compress::Gzip->new($candidate)
    or die "Cannot write $candidate: $GzipError\n";

my $header = <$in>;
die "Empty PRS base file: $input\n" unless defined $header;
$header =~ s/\r?\n$//;
my @header = split /\t/, $header, -1;
my @required = qw(CHR BP SNP A1 A2 BETA P);
die "Expected header: @required\n" unless join("\t", @header) eq join("\t", @required);

print {$out} join("\t", @required), "\n";
my %seen;
my %duplicates;
my %rejected;
my ($rows, $candidate_rows) = (0, 0);

while (my $line = <$in>) {
    $line =~ s/\r?\n$//;
    next if $line eq '';
    ++$rows;

    my @f = split /\t/, $line, -1;
    if (@f != 7) {
        reject_row('column_count', "Row $rows has " . scalar(@f) . " columns, expected 7");
        next;
    }
    my ($chr, $bp, $snp, $a1, $a2, $beta, $p) = @f;

    $chr =~ s/^chr//i;
    unless ($chr =~ /^\d+$/ && $chr >= 1 && $chr <= 22) {
        reject_row('non_autosomal', "Row $rows has non-autosomal chromosome: $f[0]");
        next;
    }
    unless ($bp =~ /^\d+$/ && $bp > 0) {
        reject_row('invalid_position', "Row $rows has invalid base-pair position: $bp");
        next;
    }

    $a1 = uc $a1;
    $a2 = uc $a2;
    unless ($a1 =~ /^[ACGT]+$/ && $a2 =~ /^[ACGT]+$/ && $a1 ne $a2) {
        reject_row('invalid_allele', "Row $rows has invalid allele: $a1/$a2");
        next;
    }

    my $canonical = "chr$chr:$bp:$a2:$a1";
    unless ($snp eq $canonical) {
        reject_row('canonical_allele_mismatch',
            "Row $rows SNP does not match A2/A1 alleles: $snp != $canonical");
        next;
    }
    if ($seen{$snp}++) {
        fail($candidate, "Duplicate SNP ID at row $rows: $snp") if $strict;
        $duplicates{$snp} = 1;
    }

    unless (finite_number($beta)) {
        reject_row('invalid_beta', "Row $rows has invalid BETA: $beta");
        next;
    }
    unless (finite_number($p) && $p > 0 && $p <= 1) {
        reject_row('invalid_p', "Row $rows has invalid P: $p");
        next;
    }

    ## emit numeric chromosomes for compatibility with PRSice and PLINK
    print {$out} join("\t", $chr, $bp, $snp, $a1, $a2, $beta, $p), "\n";
    ++$candidate_rows;
}

close $in or die "Failed while reading $input\n";
$out->close() or die "Failed while writing $candidate: $GzipError\n";

## remove every occurrence of a duplicate ID instead of choosing one arbitrarily
my $candidate_in = open_input($candidate);
my $final_out = IO::Compress::Gzip->new($tmp)
    or die "Cannot write $tmp: $GzipError\n";
my $candidate_header = <$candidate_in>;
print {$final_out} $candidate_header;
my ($accepted, $duplicate_rows) = (0, 0);
while (my $line = <$candidate_in>) {
    my @f = split /\t/, $line, -1;
    if ($duplicates{$f[2]}) {
        ++$duplicate_rows;
        next;
    }
    print {$final_out} $line;
    ++$accepted;
}
close $candidate_in or die "Failed while reading $candidate\n";
$final_out->close() or die "Failed while writing $tmp: $GzipError\n";
unlink $candidate or die "Cannot remove temporary file $candidate: $!\n";
rename $tmp, $output or die "Cannot rename $tmp to $output: $!\n";

open my $stats_fh, '>', $stats or die "Cannot write $stats: $!\n";
print {$stats_fh} "metric\tvalue\n";
print {$stats_fh} "input_rows\t$rows\n";
print {$stats_fh} "accepted_rows\t$accepted\n";
print {$stats_fh} "rejected_rows\t", $rows - $accepted, "\n";
print {$stats_fh} "duplicate_ids\t", scalar(keys %duplicates), "\n";
print {$stats_fh} "rejected_duplicate_rows\t$duplicate_rows\n";
for my $reason (sort keys %rejected) {
    print {$stats_fh} "rejected_$reason\t$rejected{$reason}\n";
}
close $stats_fh or die "Cannot close $stats: $!\n";

sub open_input {
    my ($path) = @_;
    if ($path =~ /\.gz$/) {
        ## BGZF is a concatenated gzip stream, so all members must be consumed
        my $fh = IO::Uncompress::Gunzip->new($path, MultiStream => 1)
            or die "Cannot read $path: $GunzipError\n";
        return $fh;
    }
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    return $fh;
}

sub finite_number {
    my ($value) = @_;
    return $value =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
}

sub reject_row {
    my ($reason, $message) = @_;
    fail($candidate, $message) if $strict;
    ++$rejected{$reason};
}

sub fail {
    my ($path, $message) = @_;
    unlink $path if -e $path;
    die "$message\n";
}
