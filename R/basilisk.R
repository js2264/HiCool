HiCool_args <- list(
    pkg="HiCool",
    name="env1",
    version="0.2.0",
    packages=c(
        "python==3.12.11",
        "numpy==1.26.4",
        "bowtie2==2.5.4",
        "chromosight==1.6.3", 
        "cooler==0.10.3", 
        "hicstuff==3.2.4", 
        "pairtools==1.1.3", 
        "samtools==1.22.1"
    ), 
    channels = c("conda-forge", "bioconda")
)
