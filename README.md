## Description
Pipeline in R for analysing FASTQ files using stream processing (`FastqStreamer`).Includes performance profiling (`profvis`).

## Dependencies
- R (>= 4.0)
- Packages: `BiocManager`, `ShortRead`, `ggplot2`, `dplyr`, `profvis`


## Installation
``R
install.packages(c("BiocManager", "ggplot2", "dplyr", "profvis"))
BiocManager::install(‘ShortRead’)
