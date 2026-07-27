BEGIN {
    FS = OFS = "\t"
    print "CHR", "BP", "SNP", "A1", "A2", "BETA", "SE", "P", "NCAS", "NCON", "NEFF"
}

function is_number(value) {
    return value ~ /^[-+]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$/
}

function is_nonnegative_number(value) {
    return is_number(value) && value >= 0
}

## validate the two pasted headers before processing 57 million rows
FNR == 1 {
    if (NF != 34 || $1 != "all.data.id" || $19 != "all.data.id" ||
        $20 != "src" || $24 != "pass") {
        print "ERROR: unexpected annotation/association headers" > "/dev/stderr"
        failed = 1
    }
    next
}

{
    total++

    ## the internal IDs must remain aligned row-for-row
    if (NF != 34 || $1 != $19 || $1 != total) {
        if (id_errors < 10)
            print "ERROR: row/ID mismatch at data row", total, $1, $19 > "/dev/stderr"
        id_errors++
        next
    }

    if ($24 != "Y") {
        failed_qc++
        next
    }

    chrom = $5
    sub(/^chr/, "", chrom)
    chrom_number = chrom + 0
    if (chrom !~ /^[0-9]+$/ || chrom_number < 1 || chrom_number > 22 || $8 != "A") {
        non_autosomal++
        next
    }

    ## supplied alleles are forward-strand, alphabetical A/B alleles
    nalleles = split($7, allele, "/")
    if ($18 != "+" || nalleles != 2 || allele[1] !~ /^[ACGT]$/ ||
        allele[2] !~ /^[ACGT]$/) {
        symbolic_or_non_snv++
        next
    }

    if (!is_number($22) || !is_number($23) || !is_number($21) ||
        $23 <= 0 || $21 <= 0 || $21 > 1) {
        invalid_statistics++
        next
    }

    ## imputed rows provide counts directly; genotyped rows provide genotype counts
    if ($20 == "I") {
        if (!is_nonnegative_number($25) || !is_nonnegative_number($27)) {
            invalid_sample_size++
            next
        }
        ncontrol = $25
        ncase = $27
    } else if ($20 == "G") {
        if (!is_nonnegative_number($29) || !is_nonnegative_number($30) ||
            !is_nonnegative_number($31) || !is_nonnegative_number($32) ||
            !is_nonnegative_number($33) || !is_nonnegative_number($34)) {
            invalid_sample_size++
            next
        }
        ncontrol = $29 + $30 + $31
        ncase = $32 + $33 + $34
    } else {
        invalid_source++
        next
    }

    if (ncase <= 0 || ncontrol <= 0) {
        invalid_sample_size++
        next
    }

    ## 23andMe effect is log odds per alphabetically higher B allele
    neff = 4 * ncase * ncontrol / (ncase + ncontrol)
    print chrom_number, $6, $4, allele[2], allele[1], $22, $23, $21,
        ncase, ncontrol, sprintf("%.10g", neff)
    kept++
}

END {
    print "23andMe preparation counts:" > "/dev/stderr"
    print "  total=" total > "/dev/stderr"
    print "  kept_autosomal_snvs=" kept > "/dev/stderr"
    print "  failed_23andme_qc=" failed_qc + 0 > "/dev/stderr"
    print "  symbolic_or_non_snv=" symbolic_or_non_snv + 0 > "/dev/stderr"
    print "  non_autosomal=" non_autosomal + 0 > "/dev/stderr"
    print "  invalid_statistics=" invalid_statistics + 0 > "/dev/stderr"
    print "  invalid_sample_size=" invalid_sample_size + 0 > "/dev/stderr"
    print "  invalid_source=" invalid_source + 0 > "/dev/stderr"
    print "  id_errors=" id_errors + 0 > "/dev/stderr"

    if (failed || id_errors)
        exit 2
}
