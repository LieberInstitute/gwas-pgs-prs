use warnings;

open(OUT, ">clump_jobs.txt");
foreach $chr (1..22){
	mkdir $chr;
	print OUT "perl clump.pl --chr $chr\n";
}
close(OUT);