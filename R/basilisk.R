.HiCool_dependencies <- c(
    "python==3.7.12", 
    "bowtie2==2.5.0", 
    "samtools==1.16.1", 
    "hicstuff==3.1.5", 
    "chromosight==1.6.3", 
    "cooler==0.9.1"
)

#' @importFrom basilisk BasiliskEnvironment

env_HiCool <- basilisk::BasiliskEnvironment(
    "env", 
    pkgname = "HiCool",
    packages = .HiCool_dependencies, 
    channels = c("conda-forge", "bioconda")
)
