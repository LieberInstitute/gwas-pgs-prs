#!/usr/bin/env perl
use strict;
use warnings;

use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);

my ($ld_pattern, $target_pvar, $output, $stats, $work_dir);
GetOptions(
    'ld-pattern=s'  => \$ld_pattern,
    'target-pvar=s' => \$target_pvar,
    'output=s'      => \$output,
    'stats=s'       => \$stats,
    'work-dir=s'    => \$work_dir,
) or die usage();
die usage() unless defined $ld_pattern && defined $target_pvar && defined $output;
die "LD pattern must contain {CHR}\n" unless $ld_pattern =~ /\{CHR\}/;
$stats //= "$output.stats.tsv";
$work_dir //= "$output.work";
make_path($work_dir);

# split the large target once so chromosome mapping stays memory-bounded
my %split;
for my $chr (1 .. 22) {
    open my $fh, '>', "$work_dir/target.chr$chr.tsv"
        or die "Cannot write target split for chromosome $chr: $!\n";
    $split{$chr} = $fh;
}
open my $target_in, '<', $target_pvar or die "Cannot read $target_pvar: $!\n";
while (my $line = <$target_in>) {
    next if $line =~ /^#/;
    chomp $line;
    my @f = split /\s+/, $line;
    next unless @f >= 5;
    my $chr = normalize_chr($f[0]);
    next unless defined $chr;
    print { $split{$chr} } join("\t", @f[1, 2, 3, 4]), "\n";
}
close $target_in or die "Failed while reading $target_pvar\n";
for my $fh (values %split) {
    close $fh or die "Cannot close target split: $!\n";
}

open my $out, '>', $output or die "Cannot write $output: $!\n";
my %total = (ld_rows => 0, matched => 0, missing_target => 0,
             target_allele_mismatch => 0, ambiguous_ld => 0,
             ambiguous_target => 0);
my @per_chr;

for my $chr (1 .. 22) {
    my (%target, %target_pos);
    open my $target_chr, '<', "$work_dir/target.chr$chr.tsv"
        or die "Cannot read target split for chromosome $chr: $!\n";
    while (my $line = <$target_chr>) {
        chomp $line;
        my ($pos, $id, $ref, $alt) = split /\t/, $line, -1;
        my $key = allele_key($pos, $ref, $alt);
        $target_pos{$pos} = 1;
        add_unique(\%target, $key, $id);
    }
    close $target_chr or die "Cannot close target split for chromosome $chr\n";

    my $prefix = $ld_pattern;
    $prefix =~ s/\{CHR\}/$chr/g;
    open my $bim, '<', "$prefix.bim" or die "Cannot read $prefix.bim: $!\n";
    my %count = (ld_rows => 0, matched => 0, missing_target => 0,
                 target_allele_mismatch => 0, ambiguous_ld => 0,
                 ambiguous_target => 0);
    my %ld;
    while (my $line = <$bim>) {
        chomp $line;
        my @f = split /\s+/, $line;
        next unless @f >= 6;
        ++$count{ld_rows};
        my $key = allele_key($f[3], $f[4], $f[5]);
        add_unique(\%ld, $key, $f[1]);
    }
    close $bim or die "Failed while reading $prefix.bim\n";

    for my $key (keys %ld) {
        if (!defined $ld{$key}) {
            ++$count{ambiguous_ld};
            next;
        }
        my ($pos) = split /:/, $key, 2;
        if (!exists $target{$key}) {
            ++$count{ $target_pos{$pos} ? 'target_allele_mismatch' : 'missing_target' };
            next;
        }
        if (!defined $target{$key}) {
            ++$count{ambiguous_target};
            next;
        }
        print {$out} "$ld{$key}\t$target{$key}\n";
        ++$count{matched};
    }
    push @per_chr, [$chr, map { $count{$_} }
        qw(ld_rows matched missing_target target_allele_mismatch ambiguous_ld ambiguous_target)];
    $total{$_} += $count{$_} for keys %count;
}
close $out or die "Cannot close $output: $!\n";

open my $stats_fh, '>', $stats or die "Cannot write $stats: $!\n";
print {$stats_fh} "chromosome\tld_rows\tmatched\tmissing_target\ttarget_allele_mismatch\tambiguous_ld\tambiguous_target\n";
print {$stats_fh} join("\t", @{$_}), "\n" for @per_chr;
print {$stats_fh} join("\t", 'TOTAL', map { $total{$_} }
    qw(ld_rows matched missing_target target_allele_mismatch ambiguous_ld ambiguous_target)), "\n";
close $stats_fh or die "Cannot close $stats: $!\n";

sub allele_key {
    my ($pos, $a, $b) = @_;
    ($a, $b) = (uc($a), uc($b));
    return $a lt $b ? "$pos:$a:$b" : "$pos:$b:$a";
}

sub add_unique {
    my ($mapping, $key, $value) = @_;
    if (exists $mapping->{$key} && (!defined $mapping->{$key} || $mapping->{$key} ne $value)) {
        $mapping->{$key} = undef;
    } else {
        $mapping->{$key} = $value;
    }
}

sub normalize_chr {
    my ($chr) = @_;
    $chr =~ s/^chr//i;
    return unless $chr =~ /^\d+$/ && $chr >= 1 && $chr <= 22;
    return 0 + $chr;
}

sub usage {
    return "Usage: $0 --ld-pattern PREFIX_{CHR} --target-pvar FILE " .
           "--output FILE [--stats FILE] [--work-dir DIR]\n";
}
