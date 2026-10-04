## pandas is pinned below 3: chromosight 1.6.3 cannot detect patterns with
## pandas >= 3 (copy-on-write).
HiCool_args <- list(
    pkg="HiCool",
    name="env1",
    version="0.3.0",
    packages=c(
        "python==3.12.11",
        "numpy==1.26.4",
        "pandas==2.3.3",
        "bowtie2==2.5.4",
        "chromosight==1.6.3", 
        "cooler==0.10.4", 
        "hicstuff==3.2.5", 
        "pairtools==1.1.3", 
        "samtools==1.22.1"
    ), 
    channels = c("conda-forge", "bioconda")
)
