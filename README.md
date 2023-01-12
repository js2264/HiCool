[![](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![rworkflows](https://github.com/js2264/HiCool/actions/workflows/rworkflows.yml/badge.svg)](https://github.com/js2264/HiCool/actions/workflows/rworkflows.yml)

# HiCool

The `HiCool` R/Bioconductor package provides an **end-to-end interface** to 
process and normalize Hi-C paired-end fastq reads into `.(m)cool` files.

1. The heavy lifting (fastq mapping, pairs parsing and pairs filtering) is 
performed by the underlying lightweight `hicstuff` python library 
([https://github.com/koszullab/hicstuff](https://github.com/koszullab/hicstuff)).
2. Pairs filering is done using the approach described in 
[Cournac et al., 2012](https://doi.org/10.1186/1471-2164-13-436) and implemented
in `hicstuff`.
3. `Cooler` ([https://github.com/open2c/cooler](https://github.com/open2c/cooler)) 
library is used to parse pairs into a multi-resolution, balanced `.mcool` file. 
`.(m)cool` is a compact, indexed HDF5 file format specifically tailored 
for efficiently storing HiC-based data. The `.(m)cool`  file format was 
developed by Abdennur and Mirny and 
[published in 2019](https://doi.org/10.1093/bioinformatics/btz540).
4. Internally, all these external dependencies are automatically installed and 
managed in R by a `basilisk` environment.

![](https://raw.githubusercontent.com/js2264/HiCool/master/man/figures/pipeline.png)

## Processing `.fastq` paired-end files into a `.mcool` Hi-C contact matrix

The main processing function offered in this package is `HiCool()`.
One simply needs to specify: 

- The path to each fastq file;
- The genome reference, as a `.fasta` sequence, a pre-computed `bowtie2` index 
or a supported ID (`hg38`, `mm10`, `dm6`, `R64-1-1`, `WBcel235`, `GRCz10`, 
`Galgal4`);
- The restriction enzyme(s) used for Hi-C.

```r
library(HiCool)
x <- HiCool(
    r1 = '<PATH-TO-R1.fq.gz>', 
    r2 = '<PATH-TO-R2.fq.gz>', 
    restriction = 'DpnII,HinfI', 
    genome = 'R64-1-1'
)
```

```sh
## HiCool :: Recovering bowtie2 genome index from AWS iGenomes...
## HiCool :: Initiating processing of fastq files [tmp folder: /tmp/RtmpARIRQo/DZ28I8]...
## HiCool :: Mapping fastq files...
## HiCool :: Best-suited minimum resolution automatically inferred: 1000
## HiCool :: Remove unwanted chromosomes...
## HiCool :: Generating multi-resolution .mcool file...
## HiCool :: Balancing .mcool file...
## HiCool :: Tidying up everything for you...
## HiCool :: .fastq to .mcool processing done!
## HiCool :: Check /home/rsg/repos/HiCool/HiCool folder to find the generated files
## HiCool :: Generating HiCool report. This might take a while.
## HiCool :: Report generated and available @ sample^mapped-R64-1-1^DZ28I8.html
## HiCool :: All processing successfully achieved. Congrats!
```

```r
x
```

```sh
## CoolFile object
##   .mcool file: sample^mapped-R64-1-1^55IONQ.mcool
##   resolution: 1000
##   pairs file: sample^55IONQ.pairs
##   metadata(3): log args stats
```

## Output files

```sh
## HiCool/
## |-- sample^mapped-R64-1-1^55IONQ.html
## |-- logs
## |   |-- sample^mapped-R64-1-1^55IONQ.log
## |-- matrices
## |   |-- sample^mapped-R64-1-1^55IONQ.mcool
## |-- pairs
## |   |-- sample^mapped-R64-1-1^55IONQ.pairs
## `-- plots
##     |-- sample^mapped-R64-1-1^55IONQ_event_distance.pdf
##     |-- sample^mapped-R64-1-1^55IONQ_event_distribution.pdf
```

## Reporting 

On top of processing fastq reads, HiCool provides convenient reports for 
single/multiple sample(s).

```r
x <- importHiCoolFolder(output = 'HiCool/', hash = '55IONQ')
HiCReport(x)
```

## Installation

As an R/Bioconductor package, `HiCool` should be very easy to install. The only
dependency is R (>= 4.2). In R, one can run: 

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("HiCool")
```

The first time a `HiCool()` function is executed, a `basilisk` environment 
will be automatically set up. In this environment, few dependencies will be 
installed: 

- python (pinned 3.9.1)
- numpy (pinned 1.23.4)
- bowtie2 (pinned 2.4.5)
- samtools (pinned 1.7)
- **hicstuff** (pinned 3.1.5)
- **cooler** (pinned 0.8.11)

## HiCExperiment ecosystem

`HiCool` is integrated within the `HiCExperiment` ecosystem in Bioconductor. 
Read more about the `HiCExperiment` class and handling Hi-C data in R 
[here](https://github.com/js2264/HiCExperiment).

![](https://raw.githubusercontent.com/js2264/HiCExperiment/master/man/figures/HiCExperiment_ecosystem.png)

- [HiCExperiment](https://github.com/js2264/HiCExperiment): Parsing Hi-C files in R
- [HiCool](https://github.com/js2264/HiCool): End-to-end integrated workflow to process fastq files into .cool and .pairs files
- [HiContacts](https://github.com/js2264/HiContacts): Investigating Hi-C results in R
- [HiContactsData](https://github.com/js2264/HiContactsData): Data companion package
- [fourDNData](https://github.com/js2264/fourDNData): Gateway package to 4DN-hosted Hi-C experiments
