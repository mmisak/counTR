# STRcount
![](https://img.shields.io/github/license/mmisak/STRcount)
[![trcount Release](https://img.shields.io/github/v/release/mmisak/STRcount)](https://github.com//mmisak/STRcount/releases/)
![](https://img.shields.io/github/repo-size/mmisak/STRcount)

## Overview
<p align="justify">
STRcount is a read mapping-free method to detect differential tandem repeat content between different groups of quantitative sequencing data samples. STRcount scans raw reads for repeats using Phobos, optionally filters them and then groups detected repeats and outputs their counts in a countstable, allowing for read count normalization and subsequent identification of differential tandem repeat content between groups of samples. Possible applications include the unbiased detection of tandem repeat enrichment in DNA/chromatin profiling sequencing (e.g. ChIP-seq, CUT&Run, CUT&Tag, DIP-seq) or comparison of tandem repeat content in whole genome sequencing samples.

Note: STRcount is not affiliated with Phobos or its developer Christoph Mayer.
</p>

## Features
1. STRcount offers unbiased de novo detection of differential tandem repeats in quantitative NGS data
2. Works independent of reference genomes, raw reads are the only input data
3. Repeats can be grouped by repeat length and/or perfection prior to differential comparison
4. Output can be used to generate PCA plots, heatmaps, volcano plots and others

## Requirements
- Linux - STRcount was tested on Ubuntu 18.04. Windows or macOS systems might work but were not tested thus far
- Phobos binary - STRcount was tested with the binary `phobos_64_libstdc++6` of the package `phobos-v3.3.12-linux`. Phobos is freely available for academic research and can be downloaded here: https://www.ruhr-uni-bochum.de/spezzoo/cm/cm_phobos_download.htm
- Python 3.2 or newer - STRcount was tested on Python 3.9.7

## Performing an STRcount analysis
### Running STRcount processrepeats
Run STRcount's processrepeats function on all samples of your sequencing data in FASTQ(.gz) or FASTA(.gz) format:
```
python STRcount.py processrepeats SAMPLE OUTPUT_FOLDER PHOBOS_BINARY [addtional parameters]
```

Example:
```
python STRcount.py processrepeats sample1.fastq.gz ./strcount_out/ /home/phobos/phobos-v3.3.12-linux/bin/phobos_64_libstdc++6 --grouping 'perfection:[0,100)[100,100] length:[0,40),[0,80),[80,inf]' --processes 40
```

### Running STRcount summarizecounts
Summarize the results into a count matrix that can be used for downstream analysis:
```
python STRcount.py summarizecounts SAMPLES [...] OUTPUT_FILE
```

Example:
```
python STRcount.py summarizecounts sample1.countstable.txt sample2.countstable.txt control1.countstable.txt control2.countstable.txt experiment.countmatrix.txt'
```



### Differential analysis using EdgeR
Last, we use edgeR to calculate differential tandem repeats. Here, we then use EnhancedVolcano to plot the results in a volcano plot

## Full options
### processrepeats function
```
strcount2.6_rewrite.py processrepeats [-h] [--outputprefix OUTPUT_PREFIX] [--outputtype OUTPUT_TYPE]
                                             [--processes PROCESSES_NUMBER] [--grouping GROUPING_SETTING]
                                             [--minperfection MIN_PERFECTION] [--maxperfection MAX_PERFECTION]
                                             [--minrepeatlength MIN_REP_REGION_LENGTH]
                                             [--maxrepeatlength MAX_REP_REGION_LENGTH] [--minunitsize MIN_UNIT_SIZE]
                                             [--maxunitsize MAX_UNIT_SIZE] [--minrepeatnumber MIN_REP_NUMBER]
                                             [--maxrepeatnumber MAX_REP_NUMBER]
                                             [--multirepreads MULTI_REP_READS_SETTING]
                                             [--readwhitelist READ_WHITELIST_FILE] [--readblacklist READ_BLACKLIST_FILE]
                                             [--readblocksize READ_BLOCK_SIZE]
                                             inputpath outputdirectory phobospath

positional arguments:
  inputpath             path to sequencing data file in fasta(.gz) or fastq(.gz) format
  outputdirectory       directory where the putput will be written to
  phobospath            path to Phobos executable

optional arguments:
  -h, --help            show this help message and exit
  --outputprefix OUTPUT_PREFIX
                        prefix of output files, prefix will be taken from input file, if empty string (default: )
  --outputtype OUTPUT_TYPE
                        output to generate, countstable.txt ('c'), repeatinfo.txt ('i'), repeatinfo.txt.gz ('g'),
                        concatenate the letters for multiple outputs, e.g. 'cg' for countstable.txt and
                        repeatinfo.txt.gz (default: c)
  --processes PROCESSES_NUMBER
                        number of parallel processes to be used, to automatically set to maximum number of available
                        logical cores, use 'auto' (default: auto)
  --grouping GROUPING_SETTING
                        repeat grouping settings, example: 'perfection:[0,100)[100,100] length:[0,30)[30,inf]' (default:
                        None)
  --minperfection MIN_PERFECTION
                        minimum perfection of a repeat to be considered (default: 0)
  --maxperfection MAX_PERFECTION
                        maximum perfection of a repeat to be considered (default: 100)
  --minrepeatlength MIN_REP_REGION_LENGTH
                        minimum repeat region length for a repeat to be considered (default: 0)
  --maxrepeatlength MAX_REP_REGION_LENGTH
                        maximum repeat region length for a repeat to be considered (for infinite, set value to 'inf')
                        (default: inf)
  --minunitsize MIN_UNIT_SIZE
                        minimum repeat unit size for a repeat to be considered (default: 0)
  --maxunitsize MAX_UNIT_SIZE
                        maximum repeat unit size for a repeat to be considered (for infinite, set value to 'inf')
                        (default: inf)
  --minrepeatnumber MIN_REP_NUMBER
                        minimum number of repeat units for a repeat to be considered (default: 0)
  --maxrepeatnumber MAX_REP_NUMBER
                        maximum number of repeat units for a repeat to be considered (default: inf)
  --multirepreads MULTI_REP_READS_SETTING
                        which repeat to consider in case of reads with multiple repeats (after other filters have been
                        applied), either 'all' (consider all repeats for each read), 'none' (ignore multi repeat reads),
                        'longest' (only consider the longest repeat) or 'unique_longest' (for each unique repeat unit,
                        only consider the longest) (default: all)
  --readwhitelist READ_WHITELIST_FILE
                        path to list of readnames that will not be filtered out, the rest is filtered (default: None)
  --readblacklist READ_BLACKLIST_FILE
                        path to list of readnames that will be filtered out, the rest is kept (default: None)
  --readblocksize READ_BLOCK_SIZE
                        approximate number of lines that are analyzed at once in a (parallel) process (default: 50000)
```
### summarizecounts function
```
strcount2.6_rewrite.py summarizecounts [-h] [--samplenames REPORTED_MOTIF [REPORTED_MOTIF ...]]
                                              [--reportmotif {detected,rc,combine}]
                                              inputpaths [inputpaths ...] outputfile

positional arguments:
  inputpaths            countstable.txt files to be summarized into a count matrix, can be either a directory (all
                        contained files ending in '.countstable.txt' will be summarized) or the paths to the individual
                        files
  outputfile            path to output count matrix.

optional arguments:
  -h, --help            show this help message and exit
  --samplenames REPORTED_MOTIF [REPORTED_MOTIF ...]
                        list of sample names to be used in the resulting header in the same order as input files. If
                        not set, input file names will be used (default: None)')
  --reportmotif {detected,rc,combine}
                        motif to report in the output. 'detected' to report the detected motif, 'rc' for it's reverse
                        complement (for reversely stranded data), 'combine' (combines forward and reverse complement to
                        alphebetically minimal motif) (default: detected)')
```
