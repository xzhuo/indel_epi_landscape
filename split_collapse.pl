# use strict;
use warnings;
use List::Util qw /reduce max min/;
use List::MoreUtils qw /uniq/;
BEGIN{
	my $in_f = substr $ARGV[0], 0, -4;
	my $out = $in_f.".out.txt";
	my $filterout = $in_f.".filterout.txt";
	open($fh, '>', $out) or die "Could not open file";
	open($fh2, '>', $filterout) or die "Could not open file";
}
{
	my $chr_hg = $F[0];
	my $start_hg = $F[1];
	my $end_hg = $F[2];
	my @strand_hg = split /,/, $F[3];
	my @summit_hg = grep {$_ ne "."} split /,/, $F[4];
	my @idr_hg = grep {$_ ne "."} split /,/, $F[5];
	my @chr_pt = split /,/, $F[6];
	my @start_pt = split /,/, $F[7];
	my @end_pt = split /,/, $F[8];
	my @strand_pt = split /,/, $F[9];
	my @summit_pt = grep {$_ ne "."} split /,/, $F[10];
	my @idr_pt = grep {$_ ne "."} split /,/, $F[11];
	my @uniq_chr_pt = uniq(@chr_pt);

	my $fail = 0;
	my $final_strand = "+";

	if ($#uniq_chr_pt > 0){
		$fail = 1;
	}
	else{
		if ("-" ~~ @strand_hg || "-" ~~ @strand_pt){
			$final_strand = "-";
			for my $i(0..$#strand_hg){
				if ($strand_hg[$i] eq $strand_pt[$i]){
					$fail = 1;
				}
			}
		}
	}
	unless ($fail){
		for my $i(0..$#strand_hg){
			my @starts = (@start_pt);
			my @ends = (@end_pt);
			splice @starts, $i, 1;
			splice @ends, $i, 1;
			if (@starts){
				if ($start_pt[$i] > max(@ends)+50000 || $end_pt[$i] < min(@starts)-50000){
					$fail = 1;
				}
			}
		}
	}
	if($fail){
		print $fh2 $_;
	}
	else{
		my $idr_hg = 0;
		my $fold_hg = 0;
		my $idr_pt = 0;
		my $fold_pt = 0;
		if (@idr_hg){
			my @signal_hg = map{[split /-/, $_]} @idr_hg;
			my $index_hg = reduce {$signal_hg[$a][1] > $signal_hg[$b][1]?$a:($signal_hg[$a][1] == $signal_hg[$b][1]?($signal_hg[$a][0] > $signal_hg[$b][0]?$a:$b):$b)} 0..$#signal_hg;
			$idr_hg = $signal_hg[$index_hg][1];
			$fold_hg = $signal_hg[$index_hg][0];
		}
		if (@idr_pt){
			my @signal_pt = map{[split /-/, $_]} grep {$_ =~ /-/} @idr_pt;
			my $index_pt = reduce {$signal_pt[$a][1] > $signal_pt[$b][1]?$a:($signal_pt[$a][1] == $signal_pt[$b][1]?($signal_pt[$a][0] > $signal_pt[$b][0]?$a:$b):$b)} 0..$#signal_pt;
			$idr_pt = $signal_pt[$index_pt][1];
			$fold_pt = $signal_pt[$index_pt][0];
		}
		my $strand_hg = "+";
		my $chr_pt = $uniq_chr_pt[0];
		my $start_pt = min(@start_pt);
		my $end_pt = max(@end_pt);
		my $strand_pt = $final_strand;
		print $fh "$chr_hg\t$start_hg\t$end_hg\t$strand_hg\t$fold_hg\t$idr_hg\t$chr_pt\t$start_pt\t$end_pt\t$strand_pt\t$fold_pt\t$idr_pt";
	}
}
END{
	close $fh;
	close $fh2;
}
