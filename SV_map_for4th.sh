bigWigAverageOverBed /bar/genomes/hg38/hg38.mappability.75.bigwig <(perl -lane 'BEGIN{$,="\t"}{($chr,$start,$end) = ($F[0],$F[1],$F[2]);print $chr,$start-200,$start,"index_$."}' $1) hg38.5end.mappability.tab
bigWigAverageOverBed /bar/genomes/hg38/hg38.mappability.75.bigwig <(perl -lane 'BEGIN{$,="\t"}{($chr,$start,$end) = ($F[0],$F[1],$F[2]);print $chr,$end,$end+200,"index_$."}' $1) hg38.3end.mappability.tab
bigWigAverageOverBed /bar/genomes/panTro5/panTro5.mappability.75.bigwig <(perl -lane 'BEGIN{$,="\t"}{($chr,$start,$end) = ($F[7],$F[8],$F[9]);print $chr,$start-200,$start,"index_$."}' $1) panTro5.5end.mappability.tab
bigWigAverageOverBed /bar/genomes/panTro5/panTro5.mappability.75.bigwig <(perl -lane 'BEGIN{$,="\t"}{($chr,$start,$end) = ($F[7],$F[8],$F[9]);print $chr,$end,$end+200,"index_$."}' $1) panTro5.3end.mappability.tab

paste <(sort -k1,1V hg38.5end.mappability.tab) <(sort -k1,1V hg38.3end.mappability.tab) <(sort -k1,1V panTro5.5end.mappability.tab) <(sort -k1,1V panTro5.3end.mappability.tab) | cut -f5,11,17,23 | paste $1 - | perl -lane 'print $_ if $F[-1] > 0.7 && $F[-2] > 0.7 && $F[-3] > 0.7 && $F[-4] > 0.7'
