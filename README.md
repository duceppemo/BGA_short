## BGA_short
Bacterial genome assembly pipeline for Short reads.

## Description
1. Trim reads with Fastp.
2. Assemble trimmed reads with SKESA or SPAdes.
3. Perform assemble QC with Qualimap and QUAST.

## Installation
You'll need conda (I like [minconda](https://docs.conda.io/en/latest/miniconda.html)) to create a virtual environment to manage all the dependencies. Refer to the [bioconda](https://bioconda.github.io/) for instructions on how to set your channels if not already done. I also like to use `mamba` ([here](https://mamba.readthedocs.io/en/latest/installation.html)) instead of `conda` as it's much faster and better are solving multiple dependencies.
```commandline
# Create conda enironment
conda create -n BGA_short falco=1.2.1 skesa=2.4 spades=3.15.5 bbmap=39.01 fastp=0.23.3 samtools=1.17 minimap2=2.26 \
    psutil=5.9.5 qualimap=2.2.2d quast=5.2.0 bandage=0.8.1

# Activate environment
conda activate BGA_short

# Clone repository
cd ~/prog
git clone https://github.com/duceppemo/BGA_short

# Test BGA_short
cd BGA_short
python bga_sort.py -h
```

## USAGE
```commandline
usage: python bga_short.py [-h] -i /path/to/short_reads [--ion-torrent] -o /path/to/output_folder/ [-a {skesa,spades}]
                           [-t 16] [-p 2] [-m 57] [-v]

Bacterial Genome assembly for short reads.

optional arguments:
  -h, --help            show this help message and exit
  -i /path/to/short_reads, --input /path/to/short_reads
                        Folder that contains the sort read files in fastq format. Gzipped or not. Single-end or paired-end.
                        Samples will be named according to everything before the first underscore ("_"). Mandatory
  --ion-torrent         Use this flag if fastq files are from Ion Torrent. Needed for SPAdes. Optional
  -o /path/to/output_folder/, --output /path/to/output_folder/
                        Folder to hold the result files. Mandatory.
  -a {skesa,spades}, --assembler {skesa,spades}
                        Assembly method. Default "skesa". Optional.
  -t 16, --threads 16   Number of threads. Default is maximum available(16). Optional.
  -p 2, --parallel 2    Number of samples to process in parallel. Keep low if your computer has low memory. Default is 2.
                        Optional.
  -m 57, --memory 57    Memory in GB. Default is 85% of total memory (57)
  -v, --version         show program's version number and exit
```

# Example
Assemble paired-end short reads, 3 samples in parallel, using spades:
```commandline
python bga_short.py \
    -i /home/marco/analyses/short_ass_test/fastq \
    -o /home/marco/analyses/short_ass_test/out_sapdes \
    -p 3 \
    -a spades
```
Assemble Ion Torrent single-end short reads, 2 samples in parallel, using 16 cores and 50GB of memory, with skesa:
```commandline
python bga_short.py \
    -i /home/marco/analyses/short_ass_test/fastq \
    -o /home/marco/analyses/short_ass_test/out_sapdes \
    -t 16 \
    -m 50 \
    -p 2 \
    -a skesa
```
