#!/usr/bin/perl -w
use strict;
use IO::File;
my $dir_name="S";
my $path1="e:/cattle/A2.analysis/runspace/ampeak/$dir_name";
my @threshold=(0.15, 0.2);
for my $threshold (@threshold) {
	my $fh_geno_out=IO::File->new(">e:/cattle/A2.analysis/runspace/compare/$dir_name/${dir_name}_auto_genos.$threshold.txt");
	chdir "$path1";
	for my $bianhao (2397..4141) {
		my $ampfile="$bianhao.$threshold.peakAmp.txt";
		if(-e "$path1/$ampfile"){
			my $fh_in = IO::File->new("$ampfile",'r');
			my $prefix;
			if($ampfile=~/(\w+)\.$threshold\.peakAmp\.txt/){
				$prefix=$1;
			}
			$fh_geno_out->print("$prefix");
			my $line_count=0;
			my $primary_seq;
			my $secondary_seq;
			my @sigs;
			while(<$fh_in>){
				chomp;
				my $line = $_;
				$line_count++;
				my %data;
				if($line_count>1){
					my ($A,$C,$G,$T,$ratio,$sig)=split /\s+/, $line;
					my %bases_index=(0=>"A",1=>"C",2=>"G",3=>"T");
					my @unsorted=($A,$C,$G,$T);
					my @sorted_indexes=sort {$unsorted[$b]<=>$unsorted[$a]} 0..$#unsorted;
					my @sorted_values=@unsorted[@sorted_indexes];
					my ($primary_index,$secondary_index)=@sorted_indexes[0,1];
					my $primary_base=$bases_index{$primary_index};
					$primary_seq.=$primary_base;
					my $secondary_base=$bases_index{$secondary_index};
					$secondary_seq.=$secondary_base;
					push @sigs, $sig;
				}
			}
			my @sites6=map_6sites($primary_seq);
			my $sites_count=@sites6;
			if($sites_count == 6){
				for my $site_based_on_1 (@sites6) {
					my $site_index=$site_based_on_1-1;
					my $primary_BASE=substr($primary_seq, $site_index, 1);
					my $secondary_BASE=substr($secondary_seq, $site_index, 1);
					my $sig=$sigs[$site_index];
					my $geno="$primary_BASE$primary_BASE";
					if($sig eq "TRUE"){
						$geno="$primary_BASE$secondary_BASE";
					}
					$fh_geno_out->print("\t$geno");
				}
				$fh_geno_out->print("\n");
			}
			else{
				$fh_geno_out->print("\t某位点有突变，人工检查\n");
			}
		}
		else{
			$fh_geno_out->print("$bianhao\n");
		}
	}
}

sub map_6sites {
	my $seq=shift;
	my @sites_cor_based_on_1;

	my @tags=qw(GGGCCCATCC GCCGCCTTTC GCCTGAAGTA CTCCTAAGCA TTACTGAAAG CAGCCTCTTC);
	
	for my $tag( @tags ) {
		if($seq=~/$tag/){
			my $previous=$`;
			my $matched=$&;
			my $post=$';
			my $seq_before_site="$previous$matched";
			my $coor_site=length($seq_before_site)+1;
			#print "$tag: $coor_site\n";
			push @sites_cor_based_on_1, $coor_site;
		}
		else{
			my @tag_bases=split //, $tag;
			my $match_flag=0;
			for my $base_pos (0..$#tag_bases) {
				my @new_tag_bases=@tag_bases;
				$new_tag_bases[$base_pos]="[ATGC]";
				
				my $new_tag=join("", @new_tag_bases);
				if($seq=~/$new_tag/){
					my $previous=$`;
					my $matched=$&;
					my $post=$';
					my $seq_before_site="$previous$matched";
					my $coor_site=length($seq_before_site)+1;
					#print "$new_tag: $coor_site\n";
					push @sites_cor_based_on_1, $coor_site;
					$match_flag++;
				}
			}
			if($match_flag==0){
					#print "$tag: no match\n";
			}
		}
	}
	if(@sites_cor_based_on_1 ==6){
		return @sites_cor_based_on_1;
	}
}
