# TideCluster

TideCluster is a software tool designed to identify tandem repeats in genome assemblies by utilizing Tidehunter for tandem repeat detection in smaller overlapping windows and clustering these repeats based on similarity using mmseqs2 and NCBI BLAST. The software runs in four steps as outlined below:

Tidehunter Step - In this first step, Tidehunter is run on a fasta file. This generates a GFF3 file that contains all the detected tandem repeats.

Clustering Step - This step is divided into two parts. Firstly, mmseqs2 is used to 
cluster the repeats, and then graph-based clustering is performed based on all-to-all NCBI-BLAST comparison. The GFF3 file from the Tidehunter step is updated based on the clustering results.

Annotation Step - The annotation step uses a reference library of tandem repeats. The representative consensus sequences, as reported by TideHunter, are annotated using RepeatMasker. The resulting annotation for each tandem repeat is added to the GFF3 file.

TAREAN Step - In this final step, the Tandem Repeat Analyzer is used to estimate 
consensus sequences using a k-mer-based approach on original sequences from the 
reference. This step also generates an HTML summary.

## Installation

TideCluster is available on Anaconda repository. To install TideCluster run we 
recommend to install it using [Mamba](https://github.com/mamba-org/mamba) an extremely fast replacement for the Conda package manager

In case you do not have Mamba installed, you can install it using conda to you base environment:

```bash
conda install -n base -c conda-forge mamba
```

Then install TideCluster using Mamba:
    
```bash
mamba create -n tidecluster -c conda-forge -c bioconda -c petrnovak tidecluster
```



## Usage

```help
usage: TideCluster.py [-h] [-v] {tidehunter,clustering,annotation,tarean,run_all} ...

Wrapper of TideHunter
    This script enable to run TideHunter on large fasta files in parallel. It splits
    fasta file into chunks and run TideHunter on each chunk. Identified tandem repeat 
    are then clustered, annotated and representative consensus sequences are extracted.
    
     
    

positional arguments:
  {tidehunter,clustering,annotation,tarean,run_all}
                        TideHunter wrapper
    tidehunter          Run wrapper of TideHunter
    clustering          Run clustering on TideHunter output
    annotation          Run annotation on output from clustering stepusing reference library of tandem repeats
    tarean              Run TAREAN on clusters to extract representative sequences
    run_all             Run all steps of TideCluster

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

    Example of usage:
    
    # first run tidehunter on fasta file to generate raw GFF3 output
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
    
    Library of tandem repeats for annotation step are sequences in RepeatMasker format
    where header is in format:
    
    >id#clasification
```

## TideHunter step

```help
usage: TideCluster.py tidehunter [-h] -f FASTA -pr PREFIX [-T [TIDEHUNTER_ARGUMENTS]] [-c CPU]

options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Path to reference sequence in fasta format
  -pr PREFIX, --prefix PREFIX
                        Base name for output files
  -T [TIDEHUNTER_ARGUMENTS], --tidehunter_arguments [TIDEHUNTER_ARGUMENTS]
                        additional arguments for TideHunter in quotes, default value: -p 40 -P 3000 -c 5 -e 0.25)
  -c CPU, --cpu CPU     Number of CPUs to use
```

## Clustering step

```help
usage: TideCluster.py clustering [-h] -f FASTA [-m MIN_LENGTH] -pr PREFIX [-g GFF] [-nd] [-c CPU]

options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Reference fasta
  -m MIN_LENGTH, --min_length MIN_LENGTH
                        Minimum length of tandem repeat, default (5000)
  -pr PREFIX, --prefix PREFIX
                        Prefix is used as a base name for output files.If --gff is not provided, prefix will be also usedto identify GFF file from previous tidehunter step
  -g GFF, --gff GFF     GFF3 output file from tidehunter step. If not provided the file named 'prefix_tidehunter.gff3' will be used
  -nd, --no_dust        Do not use dust filter in blastn when clustering
  -c CPU, --cpu CPU     Number of CPUs to use

```

## Annotation step

```help
usage: TideCluster.py annotation [-h] -pr PREFIX [-g GFF] [-cd CONSENSUS_DIRECTORY] -l LIBRARY [-c CPU]

options:
  -h, --help            show this help message and exit
  -pr PREFIX, --prefix PREFIX
                        Prefix is used as a base name for output files.If --gff is not provided, prefix will be also usedto identify GFF3 file from previous clustering step
  -g GFF, --gff GFF     GFF3 output file from clustering step. If not provided the file named 'prefix_clustering.gff3' will be used
  -cd CONSENSUS_DIRECTORY, --consensus_directory CONSENSUS_DIRECTORY
                        Directory with consensus sequences which are to be annotated. If not provided the directory named 'prefix_consensus' will be used
  -l LIBRARY, --library LIBRARY
                        Path to library of tandem repeats
  -c CPU, --cpu CPU     Number of CPUs to use
```

LIBRARY is fasta file with tandem repeats reference library.
Required format for sequence names  in fasta is `>seqid#class/subclass`


## TAREAN step

```help
usage: TideCluster.py tarean [-h] [-g GFF] -f FASTA -pr PREFIX [-c CPU]

options:
  -h, --help            show this help message and exit
  -g GFF, --gff GFF     GFF3 output file from annotation or clustering stepIf not provided the file named 'prefix_annotation.gff3' will be used instead. If 'prefix_annotation.gff3' is not
                        found, 'prefix_clustering.gff3' will be used
  -f FASTA, --fasta FASTA
                        Reference fasta
  -pr PREFIX, --prefix PREFIX
                        Prefix is used as a base name for output files.If --gff is not provided, prefix will be also usedto identify GFF3 files from previous clustering/annotation step
  -c CPU, --cpu CPU     Number of CPUs to use
```

## Run all steps

```help
usage: TideCluster.py run_all [-h] -f FASTA -pr PREFIX [-l LIBRARY] [-m MIN_LENGTH] [-T [TIDEHUNTER_ARGUMENTS]] [-nd] [-c CPU]

options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Reference fasta
  -pr PREFIX, --prefix PREFIX
                        Base name used for input and output files
  -l LIBRARY, --library LIBRARY
                        Path to library of tandem repeats
  -m MIN_LENGTH, --min_length MIN_LENGTH
                        Minimum length of tandem repeat (5000)
  -T [TIDEHUNTER_ARGUMENTS], --tidehunter_arguments [TIDEHUNTER_ARGUMENTS]
                        additional arguments for TideHunter in quotes, default value: -p 40 -P 3000 -c 5 -e 0.25)
  -nd, --no_dust        Do not use dust filter in blastn when clustering
  -c CPU, --cpu CPU     Number of CPUs to use
```

## Example for full pipeline

```bash
TideCluster.py tidehunter -c 40 -pr cen6_sat -f CEN6_ver_220406.fasta 
TideCluster.py clustering -c 40 -pr cen6_sat -f CEN6_ver_220406.fasta
TideCluster.py annotation -c 40 -pr cen6_sat -l library.fasta
TideCluster.py tarean -c 40 -pr cen6_sat -f CEN6_ver_220406.fasta

# or simpli run all steps at once
TideCluster.py run_all -c 40 -pr cen6_sat -f CEN6_ver_220406.fasta -l library.fasta 
```

## Output

### Tidehunter Step

- `prefix_tidehunter.gff3` - GFF3 file with tandem repeats detected by TideHunter.
- `prefix_chunks.bed` - BED file showing how the reference sequence was split into chunks for parallel processing.

### Clustering Step

- `prefix_clustering.gff3` - GFF3 file with tandem repeats identified by `mmseqs2` and `BLASTN`. 
  Tandem repeat regions in the GFF3 file are labeled by **T**andem **R**epeat **C**luster ID (TRC1, TRC2, etc.). Each TRC is described by the `repeat_type` attribute. `repeat_type` can be either TR (Tandem Repeat) or SSR (Simple Sequence Repeat).
- `prefix_clustering.gff3_1.gff3` - Intermediate file with tandem repeats clustered by `mmseqs2`.
- `prefix_consensus` - Directory with consensus sequences for each cluster as identified by TideHunter. There is one FASTA file per cluster. Each FASTA file contains all consensus sequences identified by TideHunter for a given cluster. 
- `prefix_consensus_1` - Intermediate directory with consensus sequences for each cluster as identified by `mmseqs2`. 

### Annotation Step

- `prefix_annotation.gff3` - GFF3 file with tandem repeats annotated by RepeatMasker. 
  Annotations are shown as additional attributes in the GFF3 file.
- `prefix_annotation.tsv` - Summarized annotation for each TRC cluster in a tab-delimited format. 

### TAREAN Step

- `prefix_tarean_report.html` - HTML report with tandem repeat annotations. 
- `prefix_tarean_report.tsv` - File with tandem repeat annotations in a tab-delimited format.
- `prefix_tarean` - Directory containing subdirectories with detailed TAREAN output for each TRC cluster.

## Credits

TideCluster utilizes Tidehunter [https://github.com/Xinglab/TideHunter] for tandem repeat detection.

