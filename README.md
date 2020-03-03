# indel_epi_landscape
pipeline and script for human-chimp sv affecting cis regulatory elements manuscript.

I will document my script in the coming weeks, maybe I will rewrite my bash script with snakemake.

## refine_calledSV

### Usage

```
python3 refine_calledSV.py --sv SV --bam BAM > calledSV.bed

Options:

--sv SV			Annotated structural variations. DASVC pipeline output file.
--bam BAM		BAM alignment file between the two species. Can be found in the DASVC pipeline.

```



## OrthoINDEL

### Usage

```
python3 OrthoINDEL.py [-i INPUT] [-d MIN] [-m MAX] [-b] [-s] [-n] [-p PERC]

Options:
-i INPUT, --input INPUT  The input crossmap output file.
-d MIN, --distance MIN   The distance used to define continous regions. fragments within a region are considered continous if the gap is smaller than the distance. Default is 50bp.
-m MAX, --max MAX        The distance used as upper limit of tolerated INDELs. Fragments separated by INDELs larger than the max distance will not be merged. Default is 50,000bp.
-b --broad               If the input peaks are broadpeaks (broad peaks do not have summit defined).
-s --stringent           If false, return splitted regions after final combination process. Default is True.
-n --noindel             If True, only return regions without any indel > distance.
-p PERC --perc PERC      The threshold to remove non-primary fragments (fraction of total region length).
```

