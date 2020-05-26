#!/usr/bin/env perl
die "perl $0 [depth] [reference] [length] [out] [vcf]\n" unless(@ARGV==5);

open VCF, "gzip -dc $ARGV[4]|" or die $!;
%vcf;
while (<VCF>){
chomp;
next if (/^\#/);
my @base=split/\t/,$_;
my $ref=length($base[3]);
my $alt=length($base[4]);
my $len=$ref-$alt;
if ($len >=1){
	for (my $i=1;$i<=$len;$i++){
		my $base1=substr($base[3],$alt-1+$i,1);##unuse
		my $pos=$base[1]+$alt-1+$i;#unuse
		$vcf{$pos}=$base1;
		#print "$pos\t$base1\n";
	}
}
}
close (VCF);
open IN,"$ARGV[0]" or die $!;
my %hash;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	if (exists $vcf{$t[1]}){next;}else{
	$hash{$t[1]} = 1 if($t[2] < $ARGV[2]);}
}
close IN;

open IN,"$ARGV[1]" or die $!;
open OU,">$ARGV[3]" or die $!;
my $ref;
while(<IN>){
	chomp;
	if(/>/){
		print OU "$_\n";
		#next;
	}else{
		$ref .= "$_";
	}

}
close IN;
my @seq = split //,$ref;
my $consensus;
foreach(my $i = 0; $i < scalar @seq; $i++){
	my $n = $i + 1;
##unuse
=cut
	if (exists $vcf{$n}){
	$consensus .="$vcf{$n}";
	}els
=cut
	if(exists $hash{$n}){
		$consensus .= "N";
		#print "$n\n";
	}else{
		$consensus .= "$seq[$i]";
	}
}
print OU "$consensus\n";
close OU;
