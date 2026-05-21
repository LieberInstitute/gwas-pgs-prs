# p-value informed clumping for each chr
# generate input files for PLINK2 --score

use warnings;
use Getopt::Long;
use File::Basename;
use Cwd;

GetOptions ("chr=i" => \$filechr)  
or die("Error in command line arguments\n");

# output dir
chdir $filechr;

# setup weights information
$pvalc=7; # col index for pvalue
$chrc=0; # col index for chr
$posc=1; # col index for pos-hg19 
$weic=5; # col index for effect size 
$a1c=3; # col index for allele 1, effect allele
$a2c=4;

# setup clump parameters
$clump_r2=0.1;
$clump_dist=1000;

# read snps in 1KG reference
$file="/dcs04/lieber/statsgen/shizhong/database/1KG/QC/EUR_noIndels/hg38/EUR_noIndels_chr".$filechr."_maf0.01_hg38.bim";
print "read 1KG snps...........\n";
print "$file\n";
open(IN, $file);
while(<IN>){
	chomp;
	@tokens=split(' ',$_);
	$bar=$tokens[0]."_".$tokens[3];
	# snp name used for PLINK --clump
	$snp1kg{$bar}=$tokens[1];
	if($tokens[4] lt $tokens[5]){
		$bar2=$tokens[0]."_".$tokens[3]."_".$tokens[4]."_".$tokens[5];
	}
	else{
		$bar2=$tokens[0]."_".$tokens[3]."_".$tokens[5]."_".$tokens[4];
	}
	$tag_count{$bar2}++;
}
close(IN);

# read snps in gwas data, plink2 format
$file = "/dcs05/lieber/hanlab/shan/prs/keri/data/merged_maf01_snps.pvar";
open(IN, $file);
<IN>;
while(<IN>){
	if($_ =~ /^#/ or $_ =~ /^X/){
		next;
	}
	chomp;
	@tokens=split(' ',$_);
	if($tokens[0] == $filechr){
		$bar=$tokens[0]."_".$tokens[1]; 
		# snp name used for PLINK2 --score
		$snpTargetSample{$bar}=$tokens[2];
		if($tokens[3] lt $tokens[4]){
			$bar2=$tokens[0]."_".$tokens[1]."_".$tokens[3]."_".$tokens[4];
		}
		else{
			$bar2=$tokens[0]."_".$tokens[1]."_".$tokens[4]."_".$tokens[3];
		}
		$tag_count{$bar2}++;
	}
}
close(IN);

# keep weight snps overlapping with 1KG and lieber gwas snps 
# update weight snp names to 1KG snp names
print "read weights and update snp names......\n";
$weightfile = "/dcs04/lieber/statsgen/shizhong/database/GWAS/PGC/mdd2025/gwas_hg38/pgc-mdd2025_no23andMe_eur_hg38";
if($weightfile =~ /\.gz$/){
	open(IN, "gunzip -c $weightfile |");
}
else{
	open(IN, $weightfile);
}
$out = "chr".$filechr;
open(OUT, ">$out");
print OUT "SNP\tP\n";
<IN>;
while(<IN>){
	if(/^#/){
		next;
	}
	chomp;
	@tokens=split(' ',$_);
	$bar=$tokens[$chrc]."_".$tokens[$posc];
	if($tokens[$a1c] lt $tokens[$a2c]){
		$bar2=$tokens[$chrc]."_".$tokens[$posc]."_".$tokens[$a1c]."_".$tokens[$a2c];
	}
	else{
		$bar2=$tokens[$chrc]."_".$tokens[$posc]."_".$tokens[$a2c]."_".$tokens[$a1c];
	}
	$tag_count{$bar2}++;
	if($tag_count{$bar2} == 3){
		# snp names should be 1000 genomes
	 	print OUT "$snp1kg{$bar}\t$tokens[$pvalc]\n";
		$weight{$bar} = $tokens[$weic]; # beta
		$a1{$bar} = $tokens[$a1c];
		$pval{$bar} = $tokens[$pvalc];
	}
}
close(IN);
close(OUT);

# clump
print "start clumping..........\n";
$plink="/dcl02/lieber/shan/shizhong/software/plink/plink";
$ref="/dcs04/lieber/statsgen/shizhong/database/1KG/QC/EUR_noIndels/hg38";
$command = "$plink --bfile ${ref}/EUR_noIndels_chr${filechr}_maf0.01_hg38 --clump chr$filechr --clump-p1 1 --clump-p2 1 --clump-r2 $clump_r2 --clump-kb $clump_dist --out chr$filechr";
system($command);

# output profile for clumped snps
$infile="chr".$filechr.".clumped";
$out = "profile";
$out2 = "profile_p";
open(IN, $infile);
open(OUT, ">$out");
open(OUT2, ">$out2");
<IN>;
while(<IN>){
	chomp;
	if (not $_ =~/^$/) {
		@tokens=split(' ',$_);
		$bar=$tokens[0]."_".$tokens[3];
		# snp names should be names in our gwas data
		if(defined $snpTargetSample{$bar}){
		print OUT "$snpTargetSample{$bar} $a1{$bar} $weight{$bar}\n";
		print OUT2 "$snpTargetSample{$bar}  $pval{$bar}\n"
		}
	}
}
close(IN);
close(OUT);
close(OUT2);