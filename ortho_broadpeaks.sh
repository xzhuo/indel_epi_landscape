m=$1
i=$2
j=$3
# perl -MList::Util=reduce -lane '@idr=split /,/, $F[3];@signal=split /,/, $F[4]; @peaks=split /,/, $F[5];$i=reduce{$idr[$a]>$idr[$b]?$a:($idr[$a]==$idr[$b]?($signal[$a]>$signal[$b]?$a:$b):$b)} 1..$#idr; print "$F[0]\t$F[1]\t$F[2]\t$idr[$i]\t$signal[$i]\t$peaks[$i]"' ../human_${m}_all.merge.lite.bed > human_${m}_all.merge.bestpeak.lite.bed
# perl -MList::Util=reduce -lane '@idr=split /,/, $F[3];@signal=split /,/, $F[4]; @peaks=split /,/, $F[5];$i=reduce{$idr[$a]>$idr[$b]?$a:($idr[$a]==$idr[$b]?($signal[$a]>$signal[$b]?$a:$b):$b)} 1..$#idr; print "$F[0]\t$F[1]\t$F[2]\t$idr[$i]\t$signal[$i]\t$peaks[$i]"' ../chimp_${m}_all.Pvalueall.merge.lite.bed > chimp_${m}_all.Pvalueall.merge.bestpeak.lite.bed

# sort -k4,4nr -k5,5nr chimp_${m}_all.Pvalueall.merge.bestpeak.lite.bed |head -${n}|sort -k1,1 -k2,2n -k3,3n > chimp_${m}_all.merge.bestpeak.matchingN.lite.bed
# # sort -k5,5nr chimp_${m}_all.Pvalueall.merge.bestpeak.lite.bed |head -${n}|sort -k1,1 -k2,2n -k3,3n > chimp_${m}_all.merge.bestpeak.matchingN.lite.bed

# hg_peak=human_${m}_all.merge.bestpeak.lite.bed

# pt_peak=chimp_${m}_all.merge.bestpeak.lite.bed
# pt_peak=chimp_${m}_all.merge.bestpeak.matchingbp.lite.bed
# pt_peak=chimp_${m}_all.merge.bestpeak.matchingN.lite.bed

hg_peak=$i

pt_peak=$j
# pt_peak=chimp_${m}_all.merge.lite.bed

perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[5]\t$F[4]-$F[3]\t+"' $hg_peak > human_${m}_all.merge.bestpeak.lite.strand.bed
python2 ~/github/utility/CrossMap+.py bed /bar/genomes/hg38/hg38ToPanTro5.over.chain.gz human_${m}_all.merge.bestpeak.lite.strand.bed > hg38_crossmap_panTro5.${m}.lite.strand.bed
python2 ~/github/utility/CrossMap+.py bed /bar/genomes/panTro5/panTro5ToHg38.over.chain.gz <(awk 'BEGIN{OFS="\t"}{if($7!="Fail"){print $8,$9,$10,$13,$7,$1,$2,$3,$4,$5,$6}}' hg38_crossmap_panTro5.${m}.lite.strand.bed) > hg38_crossmap_panTro5_back_hg38.${m}.lite.strand.bed
# CrossMap.py bed /bar/genomes/hg38/hg38ToPanTro5.over.chain.gz <(perl -lane 'BEGIN{$,="\t"}{print $F[0],$F[3]-1,$F[3]}' human_${m}_all.merge.bestpeak.lite.strand.bed) > hg38_crossmap_panTro5.${m}.summit.bed
perl -lane 'BEGIN{$,="\t"}{@i=split /:/, $F[4]; print $F[5],$F[6],$F[7],$F[8],$F[9],$F[10],$F[4],$F[0],$F[1],$F[2],$F[8],$F[9],$F[3] if $i[1] eq $F[12] && $i[2] == $F[13] && $i[3] == $F[14]}' hg38_crossmap_panTro5_back_hg38.${m}.lite.strand.bed > hg38_crossmap_panTro5_bothway.${m}.lite.strand.bed


# grep -v 'Fail' hg38_crossmap_panTro5.${m}.lite.strand.bed > hg38_crossmap_pass_panTro5.${m}.lite.strand.bed 
# bedtools groupby -i hg38_crossmap_pass_panTro5.${m}.lite.strand.bed -g 1,2,3,4,5,6 -c 8,9,10,13 -o distinct,collapse,collapse,distinct > hg38_crossmap_pass_panTro5.${m}.region.bed
# # perl -MList::Util=min,max -lane 'BEGIN{$,="\t"}{@s=split(/,/, $F[7]);@e=split(/,/,$F[8]); print $F[0],$F[1],$F[2],$F[3],$F[4],$F[5],$F[6],min(@s),max(@e),$F[9] if $F[6]!~/,/ && $F[9]!~/,/ && ((min(@s)==$s[0] && max(@e)==$e[-1]) or (min(@e)==$e[-1] && max(@s)==$s[0])) && abs(max(@e)-min(@s)-$F[2]+$F[1])<50000}' hg38_crossmap_pass_panTro5.${m}.region.bed > hg38_crossmap_pass_panTro5.${m}.region.filter.bed
# perl -MList::Util=min,max -lane 'BEGIN{$,="\t"}{@s=split(/,/, $F[7]);@e=split(/,/,$F[8]); print $F[0],$F[1],$F[2],$F[3],$F[4],$F[5],$F[6],min(@s),max(@e),$F[9] if $F[6]!~/,/ && $F[9]!~/,/ && issort(@s) && issort(@e) && abs(max(@e)-min(@s)-$F[2]+$F[1])<50000; sub issort{$sort=1;$order=$_[0]<$_[1]?"ord":"rev";for $i(2..$#_){if (($_[$i-1] > $_[$i] && $order eq "ord") || ($_[$i-1] < $_[$i] && $order eq "rev")){$sort=0;last}}return $sort}}' hg38_crossmap_pass_panTro5.${m}.region.bed > hg38_crossmap_pass_panTro5.${m}.region.filter.bed
# perl -MList::Util=min,max -lane 'BEGIN{$,="\t"}{@s=split(/,/, $F[7]);@e=split(/,/,$F[8]); print $F[0],$F[1],$F[2],$F[3],$F[4],$F[5] unless $F[6]!~/,/ && $F[9]!~/,/ && issort(@s) && issort(@e) && abs(max(@e)-min(@s)-$F[2]+$F[1])<50000; sub issort{$sort=1;$order=$_[0]<$_[1]?"ord":"rev";for $i(2..$#_){if (($_[$i-1] > $_[$i] && $order eq "ord") || ($_[$i-1] < $_[$i] && $order eq "rev")){$sort=0;last}}return $sort}}' hg38_crossmap_pass_panTro5.${m}.region.bed > hg38_crossmap_pass_panTro5.${m}.region.filterout.bed


python3 ~/github/utility/crossmap_parser.py -b -i hg38_crossmap_panTro5_bothway.${m}.lite.strand.bed > hg38_crossmap_panTro5.${m}.filter.bed
perl -lane 'print $_ if length($F[6]) < 6 and $F[6] ne "chrM"' hg38_crossmap_panTro5.${m}.filter.bed > hg38_crossmap_panTro5.${m}.filter.lite.bed
# bigWigAverageOverBed /bar/genomes/panTro5/panTro5.mappability.75.bigwig <(perl -lane 'print "$F[6]\t$F[7]\t$F[8]\t$F[0]:$F[1]:$F[2]:$F[3]:$F[4]:$F[5]:$F[6]:$F[7]:$F[8]:$F[9]"' hg38_crossmap_pass_panTro5.${m}.filter.bed) hg38_${m}_peak_inpanTro5.mappablepanTro5.txt -bedOut=hg38_crossmap_pass_panTro5.${m}.region.filter.mappablepanTro5.bed
# bigWigAverageOverBed /bar/genomes/hg38/hg38.mappability.75.bigwig <(perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[0]:$F[1]:$F[2]:$F[3]:$F[4]:$F[5]:$F[6]:$F[7]:$F[8]:$F[9]"' hg38_crossmap_pass_panTro5.${m}.filter.bed) hg38_${m}_peak_inpanTro5.mappablehg38.txt -bedOut=hg38_crossmap_pass_panTro5.${m}.region.filter.mappablehg38.bed
# perl -lane 'BEGIN{$,="\t"}{if(scalar @ARGV == 1){$h{$F[3]}=$F[4]}else{print $F[3],$F[4],$h{$F[3]}}}' hg38_crossmap_pass_panTro5.${m}.region.filter.mappablepanTro5.bed hg38_crossmap_pass_panTro5.${m}.region.filter.mappablehg38.bed |sort -k1,1V -k2,2n -k3,3n > hg38_crossmap_pass_panTro5.${m}.region.filter.mappable.bed
# perl -lane 'BEGIN{$,="\t"}{@i=split(/:/, $F[0]);print @i if $F[1]>0.7 && $F[2]>0.7}' hg38_crossmap_pass_panTro5.${m}.region.filter.mappable.bed > hg38_crossmap_pass_panTro5.${m}.region.filter.mappable0.7.bed


perl -lane 'BEGIN{$,="\t"}{print $F[0],$F[1],$F[2],"human_$.",$F[5],$F[3]}' hg38_crossmap_panTro5.${m}.filter.lite.bed |sort -k1,1 -k2,2n -k3,3n > hg38_peak_onhg38.${m}.lite.bed
perl -lane 'BEGIN{$,="\t"}{print $F[6],$F[7],$F[8],"human_$.",$F[5],$F[9]}' hg38_crossmap_panTro5.${m}.filter.lite.bed |sort -k1,1 -k2,2n -k3,3n > hg38_peak_onpanTro5.${m}.lite.bed


perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[5]\t$F[4]-$F[3]\t+"' $pt_peak > chimp_${m}_all.merge.bestpeak.lite.strand.bed
python2 ~/github/utility/CrossMap+.py bed /bar/genomes/panTro5/panTro5ToHg38.over.chain.gz chimp_${m}_all.merge.bestpeak.lite.strand.bed > panTro5_crossmap_hg38.${m}.lite.strand.bed
python2 ~/github/utility/CrossMap+.py bed /bar/genomes/hg38/hg38ToPanTro5.over.chain.gz <(awk 'BEGIN{OFS="\t"}{if($7!="Fail"){print $8,$9,$10,$13,$7,$1,$2,$3,$4,$5,$6}}' panTro5_crossmap_hg38.${m}.lite.strand.bed) > panTro5_crossmap_hg38_back_panTro5.${m}.lite.strand.bed
# CrossMap.py bed /bar/genomes/panTro5/panTro5ToHg38.over.chain.gz <(perl -lane 'BEGIN{$,="\t"}{print $F[0],$F[3]-1,$F[3]}' chimp_${m}_all.merge.bestpeak.lite.strand.bed) > panTro5_crossmap_hg38.${m}.summit.bed
perl -lane 'BEGIN{$,="\t"}{@i=split /:/, $F[4]; print $F[5],$F[6],$F[7],$F[8],$F[9],$F[10],$F[4],$F[0],$F[1],$F[2],$F[8],$F[9],$F[3] if $i[1] eq $F[12] && $i[2] == $F[13] && $i[3] == $F[14]}' panTro5_crossmap_hg38_back_panTro5.${m}.lite.strand.bed > panTro5_crossmap_hg38_bothway.${m}.lite.strand.bed


# grep -v 'Fail' panTro5_crossmap_hg38.${m}.lite.strand.bed > panTro5_crossmap_pass_hg38.${m}.lite.strand.bed
# bedtools groupby -i panTro5_crossmap_pass_hg38.${m}.lite.strand.bed -g 1,2,3,4,5,6 -c 8,9,10,13 -o distinct,collapse,collapse,distinct > panTro5_crossmap_pass_hg38.${m}.region.bed
# # perl -MList::Util=min,max -lane 'BEGIN{$,="\t"}{@s=split(/,/, $F[7]);@e=split(/,/,$F[8]); print $F[0],$F[1],$F[2],$F[3],$F[4],$F[5],$F[6],min(@s),max(@e),$F[9] if $F[6]!~/,/ && $F[9]!~/,/ && ((min(@s)==$s[0] && max(@e)==$e[-1]) or (min(@e)==$e[-1] && max(@s)==$s[0])) && abs(max(@e)-min(@s)-$F[2]+$F[1])<50000}' panTro5_crossmap_pass_hg38.${m}.region.bed > panTro5_crossmap_pass_hg38.${m}.region.filter.bed
# perl -MList::Util=min,max -lane 'BEGIN{$,="\t"}{@s=split(/,/, $F[7]);@e=split(/,/,$F[8]); print $F[0],$F[1],$F[2],$F[3],$F[4],$F[5],$F[6],min(@s),max(@e),$F[9] if $F[6]!~/,/ && $F[9]!~/,/ && issort(@s) && issort(@e) && abs(max(@e)-min(@s)-$F[2]+$F[1])<50000; sub issort{$sort=1;$order=$_[0]<$_[1]?"ord":"rev";for $i(2..$#_){if (($_[$i-1] > $_[$i] && $order eq "ord") || ($_[$i-1] < $_[$i] && $order eq "rev")){$sort=0;last}}return $sort}}' panTro5_crossmap_pass_hg38.${m}.region.bed > panTro5_crossmap_pass_hg38.${m}.region.filter.bed
# perl -MList::Util=min,max -lane 'BEGIN{$,="\t"}{@s=split(/,/, $F[7]);@e=split(/,/,$F[8]); print $F[0],$F[1],$F[2],$F[3],$F[4],$F[5] unless $F[6]!~/,/ && $F[9]!~/,/ && issort(@s) && issort(@e) && abs(max(@e)-min(@s)-$F[2]+$F[1])<50000; sub issort{$sort=1;$order=$_[0]<$_[1]?"ord":"rev";for $i(2..$#_){if (($_[$i-1] > $_[$i] && $order eq "ord") || ($_[$i-1] < $_[$i] && $order eq "rev")){$sort=0;last}}return $sort}}' panTro5_crossmap_pass_hg38.${m}.region.bed > panTro5_crossmap_pass_hg38.${m}.region.filterout.bed

python3 ~/github/utility/crossmap_parser.py -b -i panTro5_crossmap_hg38_bothway.${m}.lite.strand.bed > panTro5_crossmap_hg38.${m}.filter.bed
perl -lane 'print $_ if length($F[6]) < 6 and $F[6] ne "chrM"' panTro5_crossmap_hg38.${m}.filter.bed > panTro5_crossmap_hg38.${m}.filter.lite.bed
# bigWigAverageOverBed /bar/genomes/hg38/hg38.mappability.75.bigwig <(perl -lane 'print "$F[6]\t$F[7]\t$F[8]\t$F[0]:$F[1]:$F[2]:$F[3]:$F[4]:$F[5]:$F[6]:$F[7]:$F[8]:$F[9]"' panTro5_crossmap_hg38.${m}.filter.bed) panTro5_${m}_peak_inhg38.mappablehg38.txt -bedOut=panTro5_crossmap_pass_hg38.${m}.region.filter.mappablehg38.bed
# bigWigAverageOverBed /bar/genomes/panTro5/panTro5.mappability.75.bigwig <(perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[0]:$F[1]:$F[2]:$F[3]:$F[4]:$F[5]:$F[6]:$F[7]:$F[8]:$F[9]"' panTro5_crossmap_hg38.${m}.filter.bed) panTro5_${m}_peak_inhg38.mappablepanTro5.txt -bedOut=panTro5_crossmap_pass_hg38.${m}.region.filter.mappablepanTro5.bed
# perl -lane 'BEGIN{$,="\t"}{if(scalar @ARGV == 1){$h{$F[3]}=$F[4]}else{print $F[3],$F[4],$h{$F[3]}}}' panTro5_crossmap_pass_hg38.${m}.region.filter.mappablehg38.bed panTro5_crossmap_pass_hg38.${m}.region.filter.mappablepanTro5.bed |sort -k1,1V -k2,2n -k3,3n > panTro5_crossmap_pass_hg38.${m}.region.filter.mappable.bed
# perl -lane 'BEGIN{$,="\t"}{@i=split(/:/, $F[0]);print @i if $F[1]>0.7 && $F[2]>0.7}' panTro5_crossmap_pass_hg38.${m}.region.filter.mappable.bed > panTro5_crossmap_pass_hg38.${m}.region.filter.mappable0.7.bed

perl -lane 'BEGIN{$,="\t"}{print $F[0],$F[1],$F[2],"chimp_$.",$F[5],$F[3]}' panTro5_crossmap_hg38.${m}.filter.lite.bed |sort -k1,1 -k2,2n -k3,3n > panTro5_peak_onpanTro5.${m}.lite.bed
perl -lane 'BEGIN{$,="\t"}{print $F[6],$F[7],$F[8],"chimp_$.",$F[5],$F[9]}' panTro5_crossmap_hg38.${m}.filter.lite.bed |sort -k1,1 -k2,2n -k3,3n > panTro5_peak_onhg38.${m}.lite.bed

# CrossMap.py bed /bar/genomes/hg38/hg38ToPanTro5.over.chain.gz hg38_crossmap_pass_panTro5.${m}.region.filterout.bed hg38_filterout.${m}.on_panTro5.bed
# CrossMap.py bed /bar/genomes/panTro5/panTro5ToHg38.over.chain.gz panTro5_crossmap_pass_hg38.${m}.region.filterout.bed panTro5_filterout.${m}.on_hg38.bed

# cat hg38_crossmap_pass_panTro5.${m}.region.filterout.bed panTro5_filterout.${m}.on_hg38.bed > hg38.${m}_filterout.bed
# cat panTro5_crossmap_pass_hg38.${m}.region.filterout.bed hg38_filterout.${m}.on_panTro5.bed > panTro5.${m}_filterout.bed

cat <(perl -lane 'BEGIN{$,="\t"}{print $F[6],$F[7],$F[8],$F[9],".",".",$F[0],$F[1],$F[2],$F[3],$F[4],$F[5]}' panTro5_crossmap_hg38.${m}.filter.lite.bed) <(perl -lpe '$_=$_."\t.\t."' hg38_crossmap_panTro5.${m}.filter.lite.bed) |sort -k1,1 -k2,2n -k3,3n> all_peaks_on_both.${m}.bed

# filter map0.7 and collapse
# bigWigAverageOverBed /bar/genomes/hg38/hg38.mappability.75.bigwig <(perl -lane 'BEGIN{$,="\t"}{print $F[0],$F[1],$F[2],"index_$."}' all_peaks_on_both.${m}.bed) hg38.${m}.mappability.tab
# bigWigAverageOverBed /bar/genomes/panTro5/panTro5.mappability.75.bigwig <(perl -lane 'BEGIN{$,="\t"}{print $F[6],$F[7],$F[8],"index_$."}' all_peaks_on_both.${m}.bed) panTro5.${m}.mappability.tab
# paste <(sort -k1,1V hg38.${m}.mappability.tab) <(sort -k1,1V panTro5.${m}.mappability.tab) | cut -f5,11 | paste all_peaks_on_both.${m}.bed - | perl -lane 'print $_ if $F[-1] > 0.7 && $F[-2] > 0.7' > all_peaks_on_both_mappability0.7.${m}.bed
# bedtools merge -i all_peaks_on_both_mappability0.7.${m}.bed -c 4,5,6,7,8,9,10,11,12,13,14 -o collapse > all_peaks_on_both_merged0.7.collapse.${m}.bed
# bedtools merge -i all_peaks_on_both_mappability0.7.${m}.bed -c 4,5,6,7,8,9,10,11,12,13,14 -o collapse,collapse,collapse,distinct,min,max,collapse,collapse,collapse,min,min > all_peaks_on_both_merged0.7.${m}.bed
# perl -lane 'print $_ unless $F[6]=~/,/' all_peaks_on_both_merged0.7.${m}.bed> all_peaks_on_both_merged0.7.region.${m}.bed
# perl -lan /bar/xzhuo/SV/CNCC/macs2_p005idr/K4me3_pilot/split_collapse.pl all_peaks_on_both_merged0.7.collapse.${m}.bed

# collapse and filter map0.7:
bedtools merge -i all_peaks_on_both.${m}.bed -c 4,5,6,7,8,9,10,11,12 -o collapse > all_peaks_on_both.collapse.${m}.bed
perl -lan ~/github/utility/split_collapse.pl all_peaks_on_both.collapse.${m}.bed
bigWigAverageOverBed /bar/genomes/hg38/hg38.mappability.50.bigwig <(perl -lane 'BEGIN{$,="\t"}{print $F[0],$F[1],$F[2],"index_$."}' all_peaks_on_both.collapse.${m}.out.txt) hg38.${m}.mappability.tab
bigWigAverageOverBed /bar/genomes/panTro5/panTro5.mappability.50.bigwig <(perl -lane 'BEGIN{$,="\t"}{print $F[6],$F[7],$F[8],"index_$."}' all_peaks_on_both.collapse.${m}.out.txt) panTro5.${m}.mappability.tab
paste <(sort -k1,1V hg38.${m}.mappability.tab) <(sort -k1,1V panTro5.${m}.mappability.tab) | cut -f5,11 | paste all_peaks_on_both.collapse.${m}.out.txt - | perl -lane 'BEGIN{$,="\t"}{@map = splice @F,-2,2; print @F if $map[0] > 0.7 && $map[1] > 0.7}' > all_peaks_on_both.collapse.${m}.mappability0.7.txt
paste <(sort -k1,1V hg38.${m}.mappability.tab) <(sort -k1,1V panTro5.${m}.mappability.tab) | cut -f3,4,5,9,10,11 | paste all_peaks_on_both.collapse.${m}.out.txt - > all_peaks_on_both.collapse.${m}.mappability.txt
perl -lane 'BEGIN{$,="\t"}{@map = splice @F,-6,6; print @F if $map[2] > 0.7 && $map[5] > 0.7}' all_peaks_on_both.collapse.${m}.mappability.txt > all_peaks_on_both.collapse.${m}.map_pass.txt
perl -lane 'BEGIN{$,="\t"}{print @F if $F[5]>830 or $F[11]>830}' all_peaks_on_both.collapse.${m}.map_pass.txt > all_peaks_on_both.collapse.${m}.map_pass.idr_peak.txt
