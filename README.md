# BurdenRiskAB
This project describes an automated pipeline to calculate the relatevie burden risk between two populations.

## Prerequisites
```
Python (recommended version >= 3.7)
```

## Example
An example is prepared in the project **Example** folder, executable as
```
cd /path/to/brrAB/Example
sh run_example.sh 
```

## Input files
### Group information (group_infos.tsv)
```
Group  Sample
populationA sample1
populationA sample2
populationB sample3
populationB sample4
outgroupC sample5
```

## Run Brr

Brr can be easily run with:
```
python ../brrAB.py -f in.vcf.gz -A populationA -B populationB -C outgroupC -w work_dir -G group_infos.tsv
```
(Optional) if you got genotype frequence file (gt_freq_info.tsv), use:
```
python ../brrAB.py -f gt_freq_info.tsv -A populationA -B populationB -C outgroupC -w work_dir -G group_infos.tsv
```
(Optional) To specify a particular frequence (e.g. to filter the >= 0.6 in population), use:
```
--freq 0.6
```
(Optional) To specify a fix sites for jackknifes (default is 1/5 of total sites), use:
```
--fix_sites 100000
```
