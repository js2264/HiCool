.HiCool_dependencies <- c(
    "python==3.9.1", 
    "numpy==1.23.4", 
    "bowtie2==2.4.5", 
    "samtools==1.7"
)

#' @importFrom basilisk BasiliskEnvironment

env_HiCool <- basilisk::BasiliskEnvironment(
    "env", 
    pkgname = "HiCool",
    packages = .HiCool_dependencies, 
    pip = c("hicstuff==3.1.5", "cooler==0.8.11"),
    channels = c("bioconda", "conda-forge")
)
