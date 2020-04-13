#!/usr/bin/env perl
die "perl $0 [depth] [reference] [length] [out]\n" unless(@ARGV==4);

open IN,"$ARGV[0]" or die $!;
my %hash;
while(<IN>){
	chomp;
	my @t = split /\t/,$_;
	$hash{$t[1]} = 1 if($t[2] <= $ARGV[2]);
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
	if(exists $hash{$n}){
		$consensus .= "N";
	}else{
		$consensus .= "$seq[$i]";
	}
}
print OU "$consensus\n";
close OU;
