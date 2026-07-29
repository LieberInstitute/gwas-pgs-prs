#!/usr/bin/env perl
use strict;
use warnings;

use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);
use IO::Uncompress::Gunzip qw($GunzipError);

my ($base, $target_pvar, $ld_pattern, $out_dir, $chromosomes);
GetOptions(
    'base=s'         => \$base,
    'target-pvar=s'  => \$target_pvar,
    'ld-pattern=s'   => \$ld_pattern,
    'out-dir=s'      => \$out_dir,
    'chromosomes=s'  => \$chromosomes,
) or die usage();

for my $required ($base, $target_pvar, $ld_pattern, $out_dir) {
    die usage() unless defined $required;
}
$chromosomes //= join(',', 1 .. 22);
my @chromosomes = parse_chromosomes($chromosomes);
my %wanted = map { $_ => 1 } @chromosomes;

make_path($out_dir, "$out_dir/split");
split_target($target_pvar, "$out_dir/split", \%wanted);
split_base($base, "$out_dir/split", \%wanted);

for my $chr (@chromosomes) {
    my $ref_prefix = $ld_pattern;
    $ref_prefix =~ s/\{CHR\}/$chr/g;
    my $bim = "$ref_prefix.bim";
    die "Missing LD BIM: $bim\n" unless -s $bim;

    my (%ref, %ref_pos, %target, %target_pos);
    load_bim($bim, \%ref, \%ref_pos);
    load_target("$out_dir/split/target.chr$chr.tsv", \%target, \%target_pos);
    map_base($chr, "$out_dir/split/base.chr$chr.tsv", \%ref, \%ref_pos,
        \%target, \%target_pos, $out_dir);
}

sub split_target {
    my ($path, $dir, $selected) = @_;
    my %fh = output_handles($dir, 'target', $selected);
    open my $in, '<', $path or die "Cannot read $path: $!\n";
    while (my $line = <$in>) {
        next if $line =~ /^#/;
        chomp $line;
        my @f = split /\s+/, $line;
        next unless @f >= 5;
        my $chr = normalize_chr($f[0]);
        next unless defined $chr && $selected->{$chr};
        print { $fh{$chr} } join("\t", @f[1, 2, 3, 4]), "\n";
    }
    close $in or die "Failed while reading $path\n";
    close_handles(\%fh);
}

sub split_base {
    my ($path, $dir, $selected) = @_;
    my %fh = output_handles($dir, 'base', $selected);
    my $in = open_input($path);
    my $header = <$in>;
    die "Empty base file: $path\n" unless defined $header;
    while (my $line = <$in>) {
        chomp $line;
        my @f = split /\t/, $line, -1;
        next unless @f == 7;
        my $chr = normalize_chr($f[0]);
        next unless defined $chr && $selected->{$chr};
        print { $fh{$chr} } join("\t", @f[1 .. 6]), "\n";
    }
    close $in or die "Failed while reading $path\n";
    close_handles(\%fh);
}

sub output_handles {
    my ($dir, $label, $selected) = @_;
    my %fh;
    for my $chr (sort { $a <=> $b } keys %{$selected}) {
        open my $out, '>', "$dir/$label.chr$chr.tsv"
            or die "Cannot write $dir/$label.chr$chr.tsv: $!\n";
        $fh{$chr} = $out;
    }
    return %fh;
}

sub close_handles {
    my ($handles) = @_;
    for my $fh (values %{$handles}) {
        close $fh or die "Cannot close split output: $!\n";
    }
}

sub load_bim {
    my ($path, $mapping, $positions) = @_;
    open my $in, '<', $path or die "Cannot read $path: $!\n";
    while (my $line = <$in>) {
        chomp $line;
        my @f = split /\s+/, $line;
        next unless @f >= 6;
        my $key = allele_key($f[3], $f[4], $f[5]);
        $positions->{$f[3]} = 1;
        add_unique($mapping, $key, $f[1]);
    }
    close $in or die "Failed while reading $path\n";
}

sub load_target {
    my ($path, $mapping, $positions) = @_;
    open my $in, '<', $path or die "Cannot read $path: $!\n";
    while (my $line = <$in>) {
        chomp $line;
        my ($pos, $id, $ref, $alt) = split /\t/, $line, -1;
        my $key = allele_key($pos, $ref, $alt);
        $positions->{$pos} = 1;
        add_unique($mapping, $key, $id);
    }
    close $in or die "Failed while reading $path\n";
}

sub map_base {
    my ($chr, $path, $ref, $ref_pos, $target, $target_pos, $dir) = @_;
    open my $in, '<', $path or die "Cannot read $path: $!\n";
    open my $clump, '>', "$dir/chr$chr.clump.tsv"
        or die "Cannot write clump input: $!\n";
    open my $lookup, '>', "$dir/chr$chr.lookup.tsv"
        or die "Cannot write score lookup: $!\n";
    print {$clump} "SNP\tP\n";
    print {$lookup} "LD_SNP\tTARGET_SNP\tA1\tBETA\tP\n";

    my %count = (
        base_rows => 0, matched => 0, missing_ld => 0,
        ld_allele_mismatch => 0, ambiguous_ld => 0,
        missing_target => 0, target_allele_mismatch => 0,
        ambiguous_target => 0,
    );

    while (my $line = <$in>) {
        chomp $line;
        my ($pos, $snp, $a1, $a2, $beta, $p) = split /\t/, $line, -1;
        ++$count{base_rows};
        my $key = allele_key($pos, $a1, $a2);

        if (!exists $ref->{$key}) {
            ++$count{ $ref_pos->{$pos} ? 'ld_allele_mismatch' : 'missing_ld' };
            next;
        }
        if (!defined $ref->{$key}) {
            ++$count{ambiguous_ld};
            next;
        }
        if (!exists $target->{$key}) {
            ++$count{ $target_pos->{$pos} ? 'target_allele_mismatch' : 'missing_target' };
            next;
        }
        if (!defined $target->{$key}) {
            ++$count{ambiguous_target};
            next;
        }

        ++$count{matched};
        print {$clump} "$ref->{$key}\t$p\n";
        print {$lookup} join("\t", $ref->{$key}, $target->{$key}, $a1, $beta, $p), "\n";
    }
    close $in or die "Failed while reading $path\n";
    close $clump or die "Cannot close clump input: $!\n";
    close $lookup or die "Cannot close score lookup: $!\n";

    open my $stats, '>', "$dir/chr$chr.mapping-stats.tsv"
        or die "Cannot write mapping stats: $!\n";
    print {$stats} "metric\tvalue\n";
    for my $metric (qw(base_rows matched missing_ld ld_allele_mismatch ambiguous_ld
                       missing_target target_allele_mismatch ambiguous_target)) {
        print {$stats} "$metric\t$count{$metric}\n";
    }
    close $stats or die "Cannot close mapping stats: $!\n";
}

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

sub parse_chromosomes {
    my ($text) = @_;
    my %seen;
    my @result;
    for my $token (split /,/, $text) {
        if ($token =~ /^(\d+)-(\d+)$/) {
            push @result, $1 <= $2 ? ($1 .. $2) : reverse($2 .. $1);
        } elsif ($token =~ /^\d+$/) {
            push @result, 0 + $token;
        } else {
            die "Invalid chromosome token: $token\n";
        }
    }
    @result = grep { $_ >= 1 && $_ <= 22 && !$seen{$_}++ } @result;
    die "No autosomal chromosomes selected\n" unless @result;
    return @result;
}

sub open_input {
    my ($path) = @_;
    if ($path =~ /\.gz$/) {
        ## validated bases are BGZF-compatible concatenated gzip streams
        my $fh = IO::Uncompress::Gunzip->new($path, MultiStream => 1)
            or die "Cannot read $path: $GunzipError\n";
        return $fh;
    }
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    return $fh;
}

sub usage {
    return "Usage: $0 --base FILE.gz --target-pvar FILE --ld-pattern PREFIX_{CHR} " .
           "--out-dir DIR [--chromosomes 1-22]\n";
}
