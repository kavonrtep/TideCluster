# TideCluster

TideCluster is wrapper for Tidehunter and is intended for identifying tandem repeats in genome assemblies. 
It utilizes Tidehunter for tandem repeat detection in smaller overlapping windows, and then clusters these repeats based on similarity using mmseqs2 and NCBI BLAST.
TideCluster run in four steps:
1. Tidehunter step - run Tidehunter on fasta file. This step generates GFF3 file with all tandem repeats as detected by Tidehunter. 
2. Clustering step - Clustering step is performed in two steps, first step is done using mmseqs2 and second step is done using graph based clustering based on all-to-all NCBI-BLAST comparison. GFF3 file from Tidehunter step is updated based on clustering results. 
3. Annotation step - Annotation step is performed using reference library of tandem repeats. Representative consensus sequences as reported by TideHunter are annotated using RepeatMasker and resulting annotaton for each tandem repeat is added to GFF3 file.
4. TAREAN step - Tandem Repeat Analyzer is used on origin sequences from reference and consensus sequences are estimated using k-mer based approach. This steps also generate HTML summary.

## Installation

TideCluster is available on Anaconda repository. To install TideCluster run we recomend to install it using [Mamba](https://github.com/mamba-org/mamba) an extremely fast replacement for the Conda package manager

In case you do not have Mamba installed, you can install it using conda to you base environment:

```bash
conda install -n base -c conda-forge mamba
```

Then install TideCluster using Mamba:
    
```bash
mamba create -n tidecluster -c conda-forge -c bioconda -c petrnovak tidecluster
```

```bash

## Usage

```help

usage: TideCluster.py [-h] [-c CPU] {tidehunter,clustering,annotation,tarean} ...

Wrapper of TideHunter
    This script enable to run TideHunter on large fasta files in parallel. It splits
    fasta file into chunks and run TideHunter on each chunk. Identified tandem repeat 
    are then clusters, annotated and representative consensus sequences are extracted.
    
     
    

positional arguments:
  {tidehunter,clustering,annotation,tarean}
                        TideHunter wrapper
    tidehunter          Run wrapper of TideHunter
    clustering          Run clustering on TideHunter output
    annotation          Run annotation on TideHunter output using reference library of tandem repeats
    tarean              Run TAREAN on custers to extract representative sequences

options:
  -h, --help            show this help message and exit
  -c CPU, --cpu CPU     Number of CPUs to use

    Example of usage:
    
    # first run tidehunter on fasta file to generate raw gff3 output
    TideCluster.py -c 10 tidehunter -f test.fasta -pr prefix 
    
    # then run clustering on the output from previous step to cluster similar tandem repeats
    TideCluster.py -c 10 clustering -f test.fasta -pr prefix -d -m 5000
    
    # then run annotation on the clustered output to annotate clusters with reference
    # library of tandem repeats in RepeatMasker format
    TideCluster.py -c 10 annotation -pr prefix -l library.fasta
    
    # then run TAREAN on the annotated output to extract representative consensus
    # and generate html report
    TideCluster.py -c 10 tarean -f test.fasta -pr prefix
    
    Recommended parameters for TideHunter:
    short monomers: -T "-p 10 -P 39 -c 5 -e 0.25"
    long monomers: -T "-p 40 -P 3000 -c 5 -e 0.25" (default)
    
    For parallel processing include -c option before command name. 
    
    For more information about TideHunter parameters see TideHunter manual.

```

## TideHunter step

```help
usage: TideCluster.py clustering [-h] -f FASTA [-m MIN_LENGTH] -pr PREFIX [-g GFF]

options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Reference fasta
  -m MIN_LENGTH, --min_length MIN_LENGTH
                        Minimum length of tandem repeat
  -pr PREFIX, --prefix PREFIX
                        Base name used for input and output files from
  -g GFF, --gff GFF     gff output file from tidehunter ster. If not provided the file named 'prefix_tidehunter.gff3' will be used
```

## Clustering step

```help

usage: TideCluster.py clustering [-h] -f FASTA [-m MIN_LENGTH] -pr PREFIX [-g GFF] [-nd]

options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Reference fasta
  -m MIN_LENGTH, --min_length MIN_LENGTH
                        Minimum length of tandem repeat
  -pr PREFIX, --prefix PREFIX
                        Base name used for input and output files from
  -g GFF, --gff GFF     gff output file from tidehunter ster. If not provided the file named 'prefix_tidehunter.gff3' will be used
  -nd, --no_dust        Do not use dust filter in blastn when clustering
     gff output file from tidehunter ster. If not provided the file named 'prefix_tidehunter.gff3' will be used

```

## Annotation step

```help
usage: TideCluster.py annotation [-h] -pr PREFIX [-g GFF] [-cd CONSENSUS_DIRECTORY] -l LIBRARY

options:
  -h, --help            show this help message and exit
  -pr PREFIX, --prefix PREFIX
                        Base name used for input and output files from
  -g GFF, --gff GFF     gff output file from clustering step. If not provided the file named 'prefix_clustering.gff3' will be used
  -cd CONSENSUS_DIRECTORY, --consensus_directory CONSENSUS_DIRECTORY
                        Directory with consensus sequences which are to be annotated. If not provided the directory named 'prefix_consensus' will be used
  -l LIBRARY, --library LIBRARY
                        Path to library of tandem repeats

```

LIBRARY is fasta file with tandem repeats reference library.
Required format for sequence names  in fasta is `>seqid#class/subclass`


## TAREAN step

```help
usage: TideCluster.py tarean [-h] [-g GFF] -f FASTA -pr PREFIX

options:
  -h, --help            show this help message and exit
  -g GFF, --gff GFF     gff output file from annotation or clustering stepIf not provided the file named 'prefix_annotation.gff3' will be used instead. If
                        'prefix_annotation.gff3' is not found, 'prefix_clustering.gff3' will be used
  -f FASTA, --fasta FASTA
                        Reference fasta
  -pr PREFIX, --prefix PREFIX
                        Base name used for input and output files

```

## Example for full pipeline

```bash
TideCluster.py -c 40 tidehunter -pr cen6_sat -f CEN6_ver_220406.fasta 
TideCluster.py -c 40 clustering -pr cen6_sat -f CEN6_ver_220406.fasta
TideCluster.py -c 40 annotation -pr cen6_sat -l library.fasta
TideCluster.py -c 40 tarean -pr cen6_sat -f CEN6_ver_220406.fasta 



```

## Credits

TideCluster utilizes Tidehunter [https://github.com/Xinglab/TideHunter] for tandem repeat detection.
