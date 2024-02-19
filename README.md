# counTR
## Overview
<p align="justify">
counTR is a read mapping-free method to detect differential tandem repeat content between different groups of quantitative sequencing data samples. counTR scans raw reads for tandem repeats using Phobos, optionally filters them and then groups detected repeats and outputs their counts in a counts table, allowing for read count normalization and subsequent identification of differential tandem repeat content between groups of samples. Possible applications include the unbiased detection of tandem repeat enrichment in DNA/chromatin profiling sequencing (e.g. ChIP-seq, CUT&Run, CUT&Tag, DIP-seq) or comparison of tandem repeat content in whole genome sequencing samples.
</p>

## Features
- Unbiased de novo detection of differential short tandem repeat content in quantitative sequencing data
- Independent of reference genomes
- Repeats can be grouped by repeat length and/or perfection prior to differential comparison

## Requirements
- Linux - counTR was tested on Ubuntu 18.04, macOS systems might work but were not tested thus far
- Phobos binary - counTR was tested with the binary `phobos_64_libstdc++6` of the package `phobos-v3.3.12-linux`. Phobos is freely available for academic research and can be downloaded here: https://www.ruhr-uni-bochum.de/spezzoo/cm/cm_phobos_download.htm
- Python 3.2 or newer - counTR was tested on Python 3.9.7

## Installation
As of now, it is sufficient to download the counTR.py and use the program as described below. The file can be obtained from this repository or using this link: https://raw.githubusercontent.com/mmisak/counTR/main/counTR.py

## Quick start
### Running counTR processrepeats
counTR's processrepeats function detects repeats in sequencing data in in FASTQ(.gz) or FASTA(.gz) format. It can further group and count repeats and output a counts table ("Counts table files", see below, generated if `--outputtype` includes `c` as argument; default) as well as output a detailled list of detected repeats ("Repeat info files", see below, generated if `--outputtype` includes `i` as argument; not generated by default):
```
python counTR.py processrepeats SAMPLE OUTPUT_FOLDER PHOBOS_BINARY [addtional parameters]
```

**Example**:
```
python counTR.py processrepeats sample1.fastq.gz ./counTR_out/ /home/phobos/phobos-v3.3.12-linux/bin/phobos_64_libstdc++6 --grouping 'perfection:[0,100)[100,100] length:[0,40),[40,80),[80,inf]' --processes 40
```

In this example, we are running the `processrepeats` function on a sample called sample1.fastq.gz using 40 processes and grouping the detected repeats by perfection and length. The output is written into a folder called `counTR_out` in the current working directory, the folder will be created if not already present. Perfection grouping is done here by splitting the repeats in 2 groups, imperfect (perfection between 0 and lower than 100, as indicated by the brackets: `[0,100)`, square brackets indicate an inclusive range, while round brackets indicate an exclusive range) and perfect (`[100,100]`) as well as 3 length groups with the last group (`[80,inf]`) not having an upper length range limit.

**Important**:
- counTR runs much faster on machines with many cores. By default, counTR will determine the available number of cores by itself and use all of them. On a workstation PC, it can make sense to limit the number of cores by setting `--processes` to a lower number than the actually available cores.
- When using grouping, in most cases it will make sense to only use grouping ranges that a are not overlapping as illustrated in the above example (`length:[0,40),[40,80),[80,inf]`). Note the square and round brackets indicating inclusive and exclusive ranges, respecitvely. In case of overlapping ranges (e.g. `length:[0,40],[0,80],[80,inf]`), a repeat of length 40bp or 80bp would be grouped into two diffent groups and increase the repeat count of both groups by 1.
- counTR classifies tandem repeat units by their lexicographically minimal string rotation, e.g. the repeat unit of a CAGCAGCAGCAGCAGC repeat is AGC, the repeat unit of its reverse compliment is CTG
- Depending on the analysis, you might want to set the `--groupingmotif` parameter accordingly. In a non strand-aware sequencing (common WGS, ChIP-seq, CUT&Run, CUT&Tag, ..), it can make sense to set the parameter to `combine` and combine reverse complements of repeats (e.g. AGC/CTG) into a single group as these techniques commonly do not discriminate between strands. If you do not wish to combine reverse complements of repeats and group repeats by the repeat as it is detected in the repeat (e.g. for a forward-stranded sequencing), set the parameter to `detected`, to group by the reverse complement of the detected repeat (e.g. for a reversely-stranded sequencing), set the parameter to `rc`.
- If you wish to generate the repeatinfo.txt file to get detailed information (such as perfection, mismatches in comparison to a perfect repeat and more) on every single detected repeat, include `i` (for `.repeatinfo.txt`) or `g` (for `.repeatinfo.txt.gz`) in the `--outputtype` parameter. Warning: This file is usually relatively large. In our tests, its file size was usually comparable to the size of the raw reads in FASTQ format. 
  
### Running counTR summarizecounts
After obtaining `.countstable.txt` files for all of our samples, we summarize the results into a count matrix that can be used for downstream analysis:
```
python counTR.py summarizecounts OUTPUT_FILE SAMPLES [...]
```

**Example**:
```
python counTR.py summarizecounts experiment.countmatrix.txt sample1.countstable.txt sample2.countstable.txt control1.countstable.txt control2.countstable.txt
```

### Read count normalization and differential repeat content analysis
Last, we use our read count matrix for read count normalization and differential analysis. For this, analyses using edgeR, DESeq2, limma voom or Wilcoxon rank sum test are suitable.

## Full options
### counTR processrepeats parameters

| Parameter          | Description   |
| ------------------ | ------------- |
| inputpath          | Path to sequencing data file in fasta(.gz) or fastq(.gz) format. |
| outputdirectory    | Directory where the output will be written to. |
| phobospath         | Path to Phobos executable. |
| outputprefix       | Prefix of output files, prefix will be taken from input file, if empty string (default: ) |
| outputtype         | Output to generate, countstable.txt (c), repeatinfo.txt (i), repeatinfo.txt.gz (g), concatenate the letters for multiple outputs, e.g. ci (countstable.txt and repeatinfo.txt) (default: c) |
| processes          | Number of parallel processes to be used, to automatically set to maximum number of available logical cores, use 'auto' (default: auto) |
| grouping           | Repeat grouping settings, example: 'perfection:[0,100)[100,100] length:[0,30)[30,inf]' (note the single quotation marks), if 'None', repeats will be only grouped by their motif (default: None) |
| groupingmotif      | Motif to report for grouping, report the detected motif as is (detected), its reverse complement (rc), or combine forward and reverse complement (combine), all motifs are reported as their lexicographically minimal string rotation (default: detected)') |
| minperfection      | Minimum perfection of a repeat to be considered (default: 0) |
| maxperfection      | Maximum perfection of a repeat to be considered (default: 100) |
| minrepeatlength    | Minimum repeat region length for a repeat to be considered (default: 0) |
| maxrepeatlength    | Maximum repeat region length for a repeat to be considered (for infinite, set value to: inf) (default: inf) |
| minunitsize        | Minimum repeat unit size for a repeat to be considered (default: 0) |
| maxunitsize        | Maximum repeat unit size for a repeat to be considered (for infinite, set value to: inf) (default: inf) |
| mincopynumber      | Minimum number of repeat unit copies in a repeat for a repeat to be considered (default: 0) |
| maxcopynumber      | Maximum number of repeat unit copies in a repeat for a repeat to be considered (for infinite, set value to: inf) (default: inf) |
| multirepreads      | Which repeat to consider in case of reads with multiple repeats (after other filters have been applied), either 'all' (consider all repeats for each read), 'none' (ignore multi repeat reads), 'longest' (only consider the longest repeat) or 'unique_longest' (for each unique repeat unit, only consider the longest) (default: all) |
| readwhitelist      | Path to list of read names that will not be filtered out, the rest is filtered (default: None)
| readblacklist      | Path to list of read names that will be filtered out, the rest is kept (default: None)
| readchunksize      | Approximate number of lines that are analyzed at once in a (parallel) process (default: 50000)
| addphobosarguments | Add arguments to the default Phobos call (which is run with: --outputFormat 1 --reportUnit 1 --printRepeatSeqMode 2) <br><br>**Example:**  '--indelScore -4;--mismatchScore -5' (note the single quotation marks) to change the parameters Phobos uses to align detected repeats to ideal repeats <br><br>**Warning:** This command changes the way Phobos generates its output before it is passed to counTR and can result in unexpected behavior, use with caution (default: None)

### counTR summarizecounts parameters

| Parameter     | Description                                                                                                                                         |
| ------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| outputfile    | path to output count matrix |
| inputpaths    | countstable.txt files to be summarized into a count matrix |
| samplenames   | list of sample names to be used in the resulting header in the same order as input files. If not set, input file names will be used (default: None) |

## counTR output files
### Counts table files

Counts table files are generated by default by the the processrepeats function and cotain the computed repeat counts per group for one sample. They always consist of two tab-separated columns with the first column being the group and the second column being the counts of the respectiv group.

**Example:**
```
ACCG <perfection:[100.0,100.0] length:[0.0,40.0)>  1
AG <perfection:[100.0,100.0] length:[0.0,40.0)>  2
AT <perfection:[0.0,100.0) length:[40.0,80.0)>  1
AT <perfection:[100.0,100.0] length:[0.0,40.0)>  1
CTT <perfection:[100.0,100.0] length:[0.0,40.0)>  1
```
### Repeat info files

Repeat info files are relatively large tab-separated files (size is usually comparable to the size of FASTQ files used as input) containing details for each single detected repeat. These files can be generated by the processrepeats function and are only generated if specified via the `--outputtype` argument. 

**Example:**
```
unit  perfection  length  normalized_length  unit_offset  start_in_read  end_in_read  copy_number  alignment_score  mismatches  insertions  deletions  Ns  read_name  grouping  imperfections
ACCG  100.0  10  10  1  20  29  2.5  6  0  0  0  0  SRR3277954.2  ACCG <perfection:[100.0,100.0] length:[0.0,40.0)>  none
CTT  100.0  9  9  0  41  49  3.0  6  0  0  0  0  SRR3277954.8  CTT <perfection:[100.0,100.0] length:[0.0,40.0)>  none
AT  98.039  50  51  0  1  50  25.5  43  0  0  1  0  SRX3277954.15  AT <perfection:[0.0,100.0) length:[40.0,80.0)>  [6_1insG,18_19insG,27_27delG]
AG  100.0  8  8  1  1  8  4.0  6  0  0  0  0  SRR3277954.15	AG <perfection:[100.0,100.0] length:[0.0,40.0)>  none
AT  100.0  15  15  0  8  22  7.5  13  0  0  0  0  SRR3277954.15  AT <perfection:[100.0,100.0] length:[0.0,40.0)>  none
AG  100.0  29  29  0  22  50  14.5  27  0  0  0  0  SRR3277954.15  AG <perfection:[100.0,100.0] length:[0.0,40.0)>  none
```
Most of the values contained are directly taken from the Phobos output and their detailed definition can be looked up in the manual of the Phobos package. The following table briefly describes the meaning of each column:

| Column             | Description                                                                                                                          |
| ------------------ | ------------------------------------------------------------------------------------------------------------------------------------ |
| unit               | Repeat unit of the detected repeat. This will always correspond to the lexicographically minimal string rotation of the repeat as it was detected in the read. |
| perfection         | A value between 0 and 100 describing the perfection of the repeat. Perfection is calculated by Phobos using the following formula: ${ normalized\\_length − mismatches − deletions − insertions − Ns \over normalized\\_length}$ |
| length             | Length of the whole detected repeat in a read. |
| normalized_length  | Detected repeat length minus insertions, plus deletions. According to the Phobos manual, this value can be interpreted as _the number of nucleotides in the tandem repeat sequence before insertions and deletions lead to its degeneration_.
| unit_offset        | Offset of the first completed repeat unit in respect to its lexicographically minimal string rotation. <br><br>**Example:** a detected repeat might look like this: TTAGGGTTAGGGTTAGGG, the first fully completed repeat unit if read from the start is TTAGGG. counTR will label this as an AGGGTT repeat (in lexicographically minimal string rotation form). To go from TTAGGG to AGGGTT, the string has to be rotated forward 4 times, which is the value of the offset. |
| start_in_read      | Start position of the repeat in the read. |
| end_in_read        | End position of the repeat in the read. |
| copy_number        | Number of times a repeat unit is found in a repeat. Calulated by dividing the `normalized_length` by the length of the repeat unit, number can have decimals. <br><br>**Example:** The ATGC copy number in ATGCATGCAT is 2.5. |
| alignment_score    | Alignment score of the repeat when compared to an ideal repeat consisting of the detected repeat unit. |
| mismatches         | Number of mismatches in the repeat. |
| insertions         | Number of insertions in the repeat. |
| deletions          | Number of deletions in the repeat. |
| Ns                 | Number of bases that could not be called and are therefore denoted as 'N' in the read. |
| read_name          | Name of the read containing the detected repeat. |
| grouping           | Group to which a repeat got assigned. If run with `--groupingmotif` set to `detected` and without providing the `--grouping` parameter (i.e. both at default settings), this parameter will simply correspond to the detected repeat unit, since repeats will only be grouped by them. Changing the aforementioned parameters can change the group a repeat gets assigned to. This group corresponds to the group by which will have its count increased in the counts table due to the detection of the current repeat. |
| imperfections      | List of all imperfections (insertions/deletions/mutations) detected in a repeat in comparison to a perfect repeat of the same repeat unit. The resulting differences between perfect and detected repeat are reported in a style that is similar to HGVS nomenclature. <br><br> **Examples:** <br><ul><li>[6_1insG,18_19insG,27_27delG] - The first entry in the bracket (6_1insG) denotes the location change in the repeat unit. E.g. in a perfect AGGGTT repeat, a G is inserted between the T and A of the AGGGTT repeat unit (..TAGGGTT**G**AGGGTTA..). The second entry (18_19insG) is the location of the change in the ideal repeat, i.e. between base 18 and 19 of the perfect repeat. The third entry shows the location of the imperfection from the perspective of the repeat as it is detected in the read. I.e. in this case, the G at position 27 in the read would need to be deleted to match the perfect repeat. </li><li>[3C>G,38C>G,50G>C] - C to G mutation occured in the third base of the detected repeat unit, at base 38 of the perfect repeat and from the perspective of the detected repeat in the read, base 50 would need to be changed to a C to match the perfect repeat. </ul>|

### Count matrix files

Count matrix files are the result of merging multiple counts table files into a single matrix. It also adds zero counts for repeats that were only detected in a subset of all samples for the samples that do not contain them.

**Example:**
```
Repeat_group  sample1.countstable.txt  sample2.countstable.txt  control1.countstable.txt  control2.countstable.txt
ACCG <perfection:[100.0,100.0] length:[0.0,40.0)>  1  2  0  1
AG <perfection:[100.0,100.0] length:[0.0,40.0)>  2  4  1  1
AT <perfection:[0.0,100.0) length:[40.0,80.0)>  1  1  1  0
AT <perfection:[100.0,100.0] length:[0.0,40.0)>  1  0  1  2
CTT <perfection:[100.0,100.0] length:[0.0,40.0)>  1  3  0  1
```
