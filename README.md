# TideCluster

TideCluster is wrapper for Tidehunter and is intended for identifying tandem repeats in genome assemblies. It utilizes Tidehunter for tandem repeat detection in smaller overlapping windows, and then clusters these repeats based on similarity using mmseqs2 and NCBI BLAST

# Usage

```help

usage: TideCluster.py [-h] [-c CPU] {tidehunter,clustering,annotation} ...

Wrapper of TideHunter
This script enable to run TideHunter on large fasta files in parallel. It splits
fasta file into chunks and run TideHunter on each chunk. Then it merges results
and clusters similar tandem repeats.

positional arguments:
{tidehunter,clustering,annotation}
TideHunter wrapper
tidehunter          Run wrapper of TideHunter
clustering          Run clustering on TideHunter output
annotation          Run annotation on TideHunter output using reference library of tandem repeats

options:
-h, --help            show this help message and exit
-c CPU, --cpu CPU     Number of CPUs to use

Example of usage:

# first run tidehunter on fasta file to generate raw gff3 output
TideCluster.py tidehunter -f test.fasta -o test.gff3 -T "-p 40 -P 3000 -c 5 -e 0.25 
-t 46"

# then run clustering on the output to cluster similar tandem repeats
TideCluster.py clustering -f test.fasta -g test.gff3 -o test_clustered.gff3


Recommended parameters for TideHunter:
short monomers: -T "-p 10 -P 39 -c 5 -e 0.25"
long monomers: -T "-p 40 -P 3000 -c 5 -e 0.25"

For parallel processing include TideHunter -t option. 

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

TideCluster.py clustering --help
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

```

```

## Example for full pipeline

```bash
TideCluster.py -c 40 tidehunter -pr cen6_sat -f CEN6_ver_220406.fasta -T "-p 40 -P 3000 -c 5 -e 0.25"
TideCluster.py -c 40 clustering -pr cen6_sat -f CEN6_ver_220406.fasta
TideCluster.py -c 40 annotation -pr cen6_sat -l library.fasta


```

## Credits

TideCluster utilizes Tidehunter [https://github.com/Xinglab/TideHunter] for tandem repeat detection.
