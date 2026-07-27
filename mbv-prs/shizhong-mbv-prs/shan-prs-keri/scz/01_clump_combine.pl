use warnings;
use Getopt::Long;
use File::Basename;
use Cwd;

mkdir "out";
open(OUT, ">./out/profile");
open(OUT2, ">./out/profile_p");
print OUT "SNP A1 BETA\n";
print OUT2 "SNP P\n";
foreach $i (1..22){
	chdir $i;
	open(IN, "profile");
	while(<IN>){
		print OUT "$_";
	}
	close(IN);
	open(IN, "profile_p");
	while(<IN>){
		print OUT2 "$_";
	}
	close(IN);
	chdir "..";
}
close(OUT);
close(OUT2);