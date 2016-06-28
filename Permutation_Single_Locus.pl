#!/usr/bin/perl   
 use strict;    
 use warnings;
 use Getopt::Long;
 use Pod::Usage qw(pod2usage);
 use List::Util qw(sum);
##########################################################
###Default Values
my $p_freq        = 0.5;
my $permutations  = 1000;
my @observed1     = ();
my @observed2     = ();
##########################################################

##########################################################
###Usage Information and Help

=head1 SYNOPSIS

    Permutation.pl -perm 1000 -observed1 10 20 30

##############################################################
##########Chi-square test of independence##################### 
##############################################################
Information from website
(http://math.hws.edu/javamath/ryan/ChiSquare.html)


=head1 OPTIONS

	-p          frequency of allele
	-perm       number of permutations to perform
	-observed1  observed genotype counts (order AA AB BB)
	-observed2  observed genotype counts (order AA AB BB)

=cut

##########################################################
###Modified Values and Calling Help Section
GetOptions(
    q(help)               => \my $help,
    q(verbose)            => \my $verbose,
    "p=s"                 => \$p_freq,
    "perm=s"              => \$permutations,
    "observed1=s{,}"      => \@observed1,
    "observed2=s{,}"      => \@observed2
) or pod2usage(q(-verbose) => 1);
pod2usage(q(-verbose) => 1) if $help;
##########################################################

##########################################################
###Get chi-square for original observation
my @observed_expected = @{&expected(\@observed1, \@observed2)};
my $observed_chi = &chi_square(\@observed1, \@observed2, \@observed_expected);
##########################################################

##########################################################
###Generate Distribution of Two Datasets (Based on observed count totals and p_allele_frequency)
###P frequency if first calculated from dataset
$p_freq = (($observed1[0] + $observed2[0]) + 0.5*($observed1[1] + $observed2[1]))/(sum(@observed1) + sum(@observed2));
my $allele_freq_thresh = $p_freq * 100;
my @chi_square_dist = ();
while ($permutations > 0){
	my %allele1 = ("AA" => 0, "AB" => 0, "BB" => 0); ##AA => count, AB => count, BB => count
	my %allele2 = ("AA" => 0, "AB" => 0, "BB" => 0); ##AA => count, AB => count, BB => count
	###Find random distribution based on p_frequency and observed count totals
	my $total_count = sum(@observed1);
	while ($total_count > 0){
		my $allele_1 = "A";
		my $allele_2 = "A";
		my $rand_1 = 1 + int(rand(100));
		my $rand_2 = 1 + int(rand(100));
		if ($rand_1 > $allele_freq_thresh){
			$allele_1 = "B";
		}
		if ($rand_2 > $allele_freq_thresh){
			$allele_2 = "B";
		}
		if ($allele_1 ne $allele_2){
			$allele1{"AB"} ++;
		}else{
			$allele1{"$allele_1$allele_2"} ++;
		}
		$total_count --;
	}
	###Find random distribution based on p_frequency and observed count totals
	$total_count = sum(@observed2);
	while ($total_count > 0){
		my $allele_1 = "A";
		my $allele_2 = "A";
		my $rand_1 = 1 + int(rand(100));
		my $rand_2 = 1 + int(rand(100));
		if ($rand_1 > $allele_freq_thresh){
			$allele_1 = "B";
		}
		if ($rand_2 > $allele_freq_thresh){
			$allele_2 = "B";
		}
		if ($allele_1 ne $allele_2){
			$allele2{"AB"} ++;
		}else{
			$allele2{"$allele_1$allele_2"} ++;
		}
		$total_count --;
	}
	my @random_count1 = ();
	my @random_count2 = ();
	foreach my $geno (sort keys %allele1){
		push (@random_count1, $allele1{$geno});
		push (@random_count2, $allele2{$geno});
	}
	my @random_expected = @{&expected(\@random_count1, \@random_count2)};
	my $chi_square = &chi_square(\@random_count1, \@random_count2, \@random_expected);
	push (@chi_square_dist, $chi_square);
	$permutations --;
}
##########################################################

##########################################################
###Find P-value
my $count_over_observed_chi = 0;
foreach (@chi_square_dist){
	if ($_ > $observed_chi){
		$count_over_observed_chi ++;
	}
}
my $p_value = $count_over_observed_chi/(scalar @chi_square_dist);
print "Frequency of Allele 1: $p_freq\n"; 
print "P-value: $count_over_observed_chi/".@chi_square_dist." = $p_value\n";


##########################################################
###Find chi-square value
sub chi_square {
	my @geno1 = @{$_[0]};
	my @geno2 = @{$_[1]};
	my @expected_combined = @{$_[2]};
	my @combined = (@geno1, @geno2);
	my $count = 0;
	my $chi_square = 0;
	foreach (@combined){
		if ($expected_combined[$count] > 0){
			$chi_square += ($_ - $expected_combined[$count])**2/$expected_combined[$count];
		}
		$count ++;
	}
	return ($chi_square);
}

#########################################################
###Find Expected Distribution
sub expected {
	my @observed_1 = @{$_[0]};
	my @observed_2 = @{$_[1]};
	my $col_sum1  = sum(@observed_1);
	my $col_sum2  = sum(@observed_2);
	my $total_sum = $col_sum1 + $col_sum2;
	my $AA_sum    = $observed_1[0] + $observed_2[0];
	my $AB_sum    = $observed_1[1] + $observed_2[1];
	my $BB_sum    = $observed_1[2] + $observed_2[2];
	my @expected_1 = ($col_sum1*$AA_sum/$total_sum, $col_sum1*$AB_sum/$total_sum, $col_sum1*$BB_sum/$total_sum);
	my @expected_2 = ($col_sum2*$AA_sum/$total_sum, $col_sum2*$AB_sum/$total_sum, $col_sum2*$BB_sum/$total_sum);
	my @expected_combined = (@expected_1, @expected_2);
	return (\@expected_combined);
}
