# chraccd

This pipeline is still in production! Handle with care!

Multi-tool software pipeline for chromatin accessibility detection of Oxford Nanopore sequencing data. 

The software pipeline is based on the paper: Simultaneous profiling of chromatin accessibility and methylation on human cell lines with nanopore sequencing (https://www.nature.com/articles/s41592-020-01000-7#Sec9). Since the pipeline is designed for Oxford Nanopore Technologies r9.4.1 flowcell an updated pipeline needs to be created. 
The introduced sequencing technique nanoNOMe works similar to nucleosome occupancy and methylome sequencing (NOMe-seq) in Bisulfite sequencing which labels open chromatin by exposing the DNA to an exogenous M. CviPI GpC methyltransferase which methylates the GpC motif within open chromatin. 
To detect the methylation I run Oxford Nanopore Sequencing with one change, modified basecalling. I use within modified basecalling the all purpose 5mC-5hmC (dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mC_5hmC_v3). But I encourge you to use the newest available modified basecalling model to increase the accuracy. (https://github.com/nanoporetech/dorado)

Dorado basecaller can be run with the following code:
```
dorado basecaller dna_r10.4.1_e8.2_400bps_hac@v4.3.0 -x cuda:0 --modified-bases-models dna_r10.4.1_e8.2_400bps_hac@v4.3.0_5mC_5hmC@v1 <input_path_with_pod5> --reference <path_to_ref > <path_to_output (out.bam)>
```

Run the following samtools command for further processing: 
```
samtools sort -o sort.bam out.bam 

samtools index sort.bam
```
Dorado will write a samtag into the bamfile storing the methylated cytosines. To access the cytosines use the modbam2bed (https://github.com/epi2me-labs/modbam2bed) workflow by Oxford Nanopore.
```
modbam2bed --chh --chg -e -m 5mC --combine <ref> sort.bam
```
The output are two files in your working directory. 
1. mod-counts.chh.bed
2. mod-counts.chg.bed

For the next steps. We need to scrape down the modbam2bed output files.
First step is two merge the two files into one. The easiest way is:
```
cat mod-counts.chg.bed mod-counts.chh.bed > mod-counts.bed
```

The next step is to remove unused columns within mod-counts.bed. Therefore, run the script line.sh in your current working directory where you stored your mod-counts.bed

```
bash line.sh
```
The output file (out.bed) will be used as input for coverage2cytosine. 

Modbam2bed can currently not filter for GpC motif and can correctly assign the motif to its strand position. Therefore, I have to run coverage2cytosine a tool from Bismark (https://github.com/FelixKrueger/Bismark/tree/master). 
```
coverage2cytosine --genome_folder ref/ --dir . --output cov2cyt --gc --nome-seq <input.file>
```

To translate the non-CpG motif into a GpC motif I run the following lookup to remove all other motifs but GpC. (the input is cov2cyt.NOMe.GpC_report.txt the output file of the previous command.

```
Rscript lookup.R
```

The output file loci.bed needs to split into chromosomes like (loci.chr1.bed) and folders for each chromosome needs to be created (like chr1).

Then the final peak calling can be run which will create a bed file which can be viewed in e.g. IGV.

```
Rscript peak_calling.R
```







