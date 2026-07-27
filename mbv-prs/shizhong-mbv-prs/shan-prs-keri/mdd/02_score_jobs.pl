use warnings;

open(OUT, ">score_jobs.txt");
print OUT "plink2 --pfile /dcs05/lieber/hanlab/shan/prs/keri/data/merged_maf01_snps ";
print OUT "--score ./out/profile header ignore-dup-ids ";
print OUT "--q-score-range prange ./out/profile_p header --threads 1 --out ./out/score\n";
close(OUT);