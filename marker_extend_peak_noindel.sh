m=$1

# merged MACS2 peak files:
hg_peak=/saloon/WangLibrary/xiaoyu/SV/CNCC/macs2_p005idr/human_${m}_all.Pvalueall.merge.lite.bed
pt_peak=/saloon/WangLibrary/xiaoyu/SV/CNCC/macs2_p005idr/chimp_${m}_all.Pvalueall.merge.lite.bed
# liftOver files:
hg38ToPanTro5=/bar/genomes/hg38/hg38ToPanTro5.over.chain.gz
panTro5ToHg38=/bar/genomes/panTro5/panTro5ToHg38.over.chain.gz
# mappability files:
hg38mappability=/bar/genomes/hg38/hg38.mappability.50.bigwig
panTro5mappability=/bar/genomes/panTro5/panTro5.mappability.50.bigwig

# run CrossMap+ from h to c:
perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[5]\t$F[4]-$F[3]\t+"' $hg_peak > human_${m}_all.merge.bestpeak.lite.strand.bed
python2 ~/github/indel_epi_landscape/CrossMap+.py bed $hg38ToPanTro5 human_${m}_all.merge.bestpeak.lite.strand.bed > hg38_crossmap_panTro5.${m}.lite.strand.bed
python2 ~/github/indel_epi_landscape/CrossMap+.py bed $panTro5ToHg38 <(awk 'BEGIN{OFS="\t"}{if($7!="Fail"){print $8,$9,$10,$13,$7,$1,$2,$3,$4,$5,$6}}' hg38_crossmap_panTro5.${m}.lite.strand.bed) > hg38_crossmap_panTro5_back_hg38.${m}.lite.strand.bed
perl -lane 'BEGIN{$,="\t"}{@i=split /:/, $F[4]; print $F[5],$F[6],$F[7],$F[8],$F[9],$F[10],$F[4],$F[0],$F[1],$F[2],$F[8],$F[9],$F[3] if $i[1] eq $F[12] && $i[2] == $F[13] && $i[3] == $F[14]}' hg38_crossmap_panTro5_back_hg38.${m}.lite.strand.bed > hg38_crossmap_panTro5_bothway.${m}.lite.strand.bed

# run crossmap parser and filter:
python3 ~/github/indel_epi_landscape/crossmap_parser.py -d 20 -m 50 -i hg38_crossmap_panTro5_bothway.${m}.lite.strand.bed > hg38_crossmap_panTro5.${m}.filter.bed
perl -lane 'print $_ if length($F[6]) < 6 and $F[6] ne "chrM"' hg38_crossmap_panTro5.${m}.filter.bed > hg38_crossmap_panTro5.${m}.filter.lite.bed
perl -lane 'BEGIN{$,="\t"}{print $F[0],$F[1],$F[2],"human_$.",$F[5],$F[3]}' hg38_crossmap_panTro5.${m}.filter.lite.bed |sort -k1,1 -k2,2n -k3,3n > hg38_peak_onhg38.${m}.lite.bed
perl -lane 'BEGIN{$,="\t"}{print $F[6],$F[7],$F[8],"human_$.",$F[5],$F[9]}' hg38_crossmap_panTro5.${m}.filter.lite.bed |sort -k1,1 -k2,2n -k3,3n > hg38_peak_onpanTro5.${m}.lite.bed

# run CrossMap+ from c to h:
perl -lane 'print "$F[0]\t$F[1]\t$F[2]\t$F[5]\t$F[4]-$F[3]\t+"' $pt_peak > chimp_${m}_all.merge.bestpeak.lite.strand.bed
python2 ~/github/indel_epi_landscape/CrossMap+.py bed $panTro5ToHg38 chimp_${m}_all.merge.bestpeak.lite.strand.bed > panTro5_crossmap_hg38.${m}.lite.strand.bed
python2 ~/github/indel_epi_landscape/CrossMap+.py bed $hg38ToPanTro5 <(awk 'BEGIN{OFS="\t"}{if($7!="Fail"){print $8,$9,$10,$13,$7,$1,$2,$3,$4,$5,$6}}' panTro5_crossmap_hg38.${m}.lite.strand.bed) > panTro5_crossmap_hg38_back_panTro5.${m}.lite.strand.bed
perl -lane 'BEGIN{$,="\t"}{@i=split /:/, $F[4]; print $F[5],$F[6],$F[7],$F[8],$F[9],$F[10],$F[4],$F[0],$F[1],$F[2],$F[8],$F[9],$F[3] if $i[1] eq $F[12] && $i[2] == $F[13] && $i[3] == $F[14]}' panTro5_crossmap_hg38_back_panTro5.${m}.lite.strand.bed > panTro5_crossmap_hg38_bothway.${m}.lite.strand.bed

# run crossmap parser and filter:
python3 ~/github/indel_epi_landscape/crossmap_parser.py -d 20 -m 50 -i panTro5_crossmap_hg38_bothway.${m}.lite.strand.bed > panTro5_crossmap_hg38.${m}.filter.bed
perl -lane 'print $_ if length($F[6]) < 6 and $F[6] ne "chrM"' panTro5_crossmap_hg38.${m}.filter.bed > panTro5_crossmap_hg38.${m}.filter.lite.bed
perl -lane 'BEGIN{$,="\t"}{print $F[0],$F[1],$F[2],"chimp_$.",$F[5],$F[3]}' panTro5_crossmap_hg38.${m}.filter.lite.bed |sort -k1,1 -k2,2n -k3,3n > panTro5_peak_onpanTro5.${m}.lite.bed
perl -lane 'BEGIN{$,="\t"}{print $F[6],$F[7],$F[8],"chimp_$.",$F[5],$F[9]}' panTro5_crossmap_hg38.${m}.filter.lite.bed |sort -k1,1 -k2,2n -k3,3n > panTro5_peak_onhg38.${m}.lite.bed

# cat all peaks on both species:
cat <(perl -lane 'BEGIN{$,="\t"}{print $F[6],$F[7],$F[8],$F[9],".",".",$F[0],$F[1],$F[2],$F[3],$F[4],$F[5]}' panTro5_crossmap_hg38.${m}.filter.lite.bed) <(perl -lpe '$_=$_."\t.\t."' hg38_crossmap_panTro5.${m}.filter.lite.bed) |sort -k1,1 -k2,2n -k3,3n> all_peaks_on_both.${m}.bed

# collapse and filter map0.7:
bedtools merge -i all_peaks_on_both.${m}.bed -c 4,5,6,7,8,9,10,11,12 -o collapse > all_peaks_on_both.collapse.${m}.bed
perl -lan ~/github/indel_epi_landscape/split_collapse.pl all_peaks_on_both.collapse.${m}.bed
bigWigAverageOverBed $hg38mappability <(perl -lane 'BEGIN{$,="\t"}{print $F[0],$F[1],$F[2],"index_$."}' all_peaks_on_both.collapse.${m}.out.txt) hg38.${m}.mappability.tab
bigWigAverageOverBed $panTro5mappability <(perl -lane 'BEGIN{$,="\t"}{print $F[6],$F[7],$F[8],"index_$."}' all_peaks_on_both.collapse.${m}.out.txt) panTro5.${m}.mappability.tab
paste <(sort -k1,1V hg38.${m}.mappability.tab) <(sort -k1,1V panTro5.${m}.mappability.tab) | cut -f5,11 | paste all_peaks_on_both.collapse.${m}.out.txt - | perl -lane 'BEGIN{$,="\t"}{@map = splice @F,-2,2; print @F if $map[0] > 0.7 && $map[1] > 0.7}' > all_peaks_on_both.collapse.${m}.mappability0.7.txt
paste <(sort -k1,1V hg38.${m}.mappability.tab) <(sort -k1,1V panTro5.${m}.mappability.tab) | cut -f3,4,5,9,10,11 | paste all_peaks_on_both.collapse.${m}.out.txt - > all_peaks_on_both.collapse.${m}.mappability.txt
perl -lane 'BEGIN{$,="\t"}{@map = splice @F,-6,6; print @F if $map[2] > 0.7 && $map[5] > 0.7}' all_peaks_on_both.collapse.${m}.mappability.txt > all_peaks_on_both.collapse.${m}.map_pass.txt
perl -lane 'BEGIN{$,="\t"}{print @F if $F[5]>830 or $F[11]>830}' all_peaks_on_both.collapse.${m}.map_pass.txt > all_peaks_on_both.collapse.${m}.map_pass.idr_peak.noindel.txt
