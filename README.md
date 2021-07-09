# SWAM - Smartly weighted averaging across multiple tissues

## Overview

SWAM is an gene expression imputation method which combines information from multiple sources to boost accuracy of imputed gene expression levels 

## Getting Started

### Prerequisites

To use SWAM, you need to have the following tools installed in your system.
* `perl` (version 5 or later recommended)
* `python` (version 2.7 is recommended. `numpy` package is also required)
* `R` (version 3.4 or later recommended)
* `htslib` : The binary `tabix` should be include in your `$PATH`. Type `tabix` in your command line to check.
* `sqlite3` : (3.22 or later recommended)

(The software `prediXcan`, which is currently deprecated, was originally downloaded from https://github.com/hakyimlab/PrediXcan and copied in this repository)


### Cloning the repository

Next, clone the repository using the following command:

```
git clone https://github.com/aeyliu/SWAM.git
```

### Preparing Input Files

* genotypes : the genotype file should be dosages with the following columns - chromosome, RSID, position, reference allele, alternate allele, minor allele frequency, dosages
```
#example of genotypes file (tab separated)
#chromosome   rsid          position  ref alt   maf   dosage_individual1  dosage_individual2    ...
chr22         rs141944226   16286442  C    G    0.08       0.5                1  
chr22         rs150703810   16286465  G    C    0.07       0.5                2
...
```
* samples : a text file with sample IDs should be included corresponding to the columns of the genotype file
```
#example of samples file (tab separated)
ID1     ID1
ID2     ID2
...
```
* single-tissue database files : these .db files are compatible with the format generated by prediXcan (for example, GTEx derived tissues can be obtained from http://predictdb.org/). If using your own expression/genotype data, refer to prediXcan pipeline on how to generate predictDB-style models (https://github.com/hakyimlab/PrediXcan)
* expression file for target tissue : measured expression file with following columns - gene ensemble ID, expression values (for each sample)
```
#example of expression file
Id                    ID1                    ID2                    ID3                 ...
ENSG00000177663.9     -0.120170828099355     -0.905468215149992     1.62413112667095 
ENSG00000069998.8     -1.35973738393861      0.640666889919105      -1.00885646145564
...
```

### Running SWAM (example)
* An example to run SWAM is included in the /SWAM/examples folder
* Example of input files can be examined in the /SWAM/examples/sample folder
* To run the example, simply replace {} with the directory where you cloned this repository
* Also remember to specify the correct prediXcan installation path
```
## Modify these environment variables to conform your settings
export SWAMDIR=/path/to/SWAM
## Run this command to run example code
${SWAMDIR}/scripts/swam \
--directory ${SWAMDIR}/examples/sample/GTEx-V6p-1KG-2016-11-16 \
--name TW_Cells_EBV-transformed_lymphocytes_0.5_1KG \
--expr ${SWAMDIR}/examples/sample/Cells_EBV-transformed_lymphocytes_Analysis.chr22.expr.txt \
--geno ${SWAMDIR}/examples/sample/genotypes \
--PrediXcan-path ${SWAMDIR}/scripts/PrediXcan.py \
--num-cpu 4 \
--out ${SWAMDIR}/examples/lcl
```
* To run SWAM, either --directory or --index file must be specified (either are fine)
* The file format for the index file can be examined in the output {}/SWAM/examples/lcl/index.txt, and contains two columns: first column is the name of each tissue, and second column is the file path of its corresponding prediction model

### Commands
Current list of commands for SWAM:

* _--index_

  Index file containing the list of tissue-specific training models. Each line should have two columns: 1) Tissue name and 2) file path for prediction model
* _--directory_

  Directory containing all prediction models to be used  by SWAM. If index file is not specified, will be generated automatically
* _--name_

  Name of the target tissue. Must be included in the index file
* _--expr_

  Measured expression data for the target tissue, see input files for further details
* _--geno_

  Genotype files in gzipped dosage format, see input files for further details
* _--out_

  Prefix of output files

Additional options:
* _--num-cpu_
  
  Assign number of CPUs for parallelization (this will be helpful when calculating covariance file)
* _--PrediXcan-path_

  Path to PrediXcan software tool
* _--sqlite3-path_

  Path to sqlite3 software tool
* _--Rscript-path_

  Path to Rscript software tool
* _--tabix-path_

  Path to tabix software tool
* _--keep-files_

  Option to keep intermediate files
* _--cal-cov_

  Calculate covariate matrix, which is needed for metaXcan

## Command-line documentation of SWAM

A brief usage on how to use SWAM can be obtained by running it without any arguments.

```sh
${SWAMDIR}/scripts/swam

ERROR: Missing required arguments. Please see the usage below
Usage:
    /net/fantasia/home/hmkang/code/git/SWAM/scripts/swam [options]

     General Options:
      -help             Print out brief help message [OFF]
      -man              Print the full documentation in man page style [OFF]

     Required Options:
      -index STR        Index file containing the list of tissue-specific training models. Each line should have [TISSUE_NAME] [Path to PredictDB-formatted file []
      -directory STR    Directory containing all prediction models to be used  by SWAM. If index file is not specified, will be generated from directory []
      -name STR         Name of the target tissue. Must be included in the index file []
      -expr STR         Measured expression data for the target tissue (in PrediXcan format). First line has sample IDs, and from the second line [GENE_NAME] [EXPR_FOR_SAMPLE_1] [EXPR_FOR_SAMPLE_2] ... []
      -geno STR         Genotype files in gzipped dosage format in PrediXcan format []
      -out STR          Prefix of output files []

     Additional Options:
      -num-cpu STR      Assign number of CPUs for parallelization [1]
      -PrediXcan-path STRPath to PrediXcan software tool [PrediXcan.py]
      -Rscript-path STR Path to Rscript tool [Rscript]
      -sqlite3-path STR Path to sqlite3 tool [sqlite3]
      -tabix-path STR   Path to tabix tool [tabix]
      -keep-files       Option to keep intermediate files [OFF]
      -cal-cov          Calculate covariate matrix [OFF]
```


The full command line documentation of SWAM can be obtained using `--help` option as follows

```sh
${SWAMDIR}/scripts/swam --help

Usage:
    /net/fantasia/home/hmkang/code/git/SWAM/scripts/swam [options]

     General Options:
      -help             Print out brief help message [ON]
      -man              Print the full documentation in man page style [OFF]

     Required Options:
      -index STR        Index file containing the list of tissue-specific training models. Each line should have [TISSUE_NAME] [Path to PredictDB-formatted file []
      -directory STR    Directory containing all prediction models to be used  by SWAM. If index file is not specified, will be generated from directory []
      -name STR         Name of the target tissue. Must be included in the index file []
      -expr STR         Measured expression data for the target tissue (in PrediXcan format). First line has sample IDs, and from the second line [GENE_NAME] [EXPR_FOR_SAMPLE_1] [EXPR_FOR_SAMPLE_2] ... []
      -geno STR         Genotype files in gzipped dosage format in PrediXcan format []
      -out STR          Prefix of output files []

     Additional Options:
      -num-cpu STR      Assign number of CPUs for parallelization [1]
      -PrediXcan-path STRPath to PrediXcan software tool [PrediXcan.py]
      -Rscript-path STR Path to Rscript tool [Rscript]
      -sqlite3-path STR Path to sqlite3 tool [sqlite3]
      -tabix-path STR   Path to tabix tool [tabix]
      -keep-files       Option to keep intermediate files [OFF]
      -cal-cov          Calculate covariate matrix [OFF]

Options:
    -help   Print a brief help message and exits

    -man    Prints the manual page and exits

    --help [ON]
            Print a help message and exits

    --man [OFF]
            Prints a manual page and exits upon typing 'q'

    --index STR []
            Index file containing the list of tissue-specific training
            models. Each line should have [TISSUE_NAME] [Path to
            PredictDB-formatted file

    --directory STR []
            Directory containing all prediction models to be used by SWAM.
            If index file is not specified, will be generated from directory

    --name STR []
            Name of the target tissue. Must be included in the index file

    --expr STR []
            Measured expression data for the target tissue (in PrediXcan
            format). First line has sample IDs, and from the second line
            [GENE_NAME] [EXPR_FOR_SAMPLE_1] [EXPR_FOR_SAMPLE_2] ...

    --geno STR []
            Genotype files in gzipped dosage format in PrediXcan format

    --out STR []
            Prefix of output files

    --num-cpu STR [1]
            Assign number of CPUs for parallelization

    --PrediXcan-path STR [PrediXcan.py]
            Path to PrediXcan software tool

    --Rscript-path STR [Rscript]
            Path to Rscript tool

    --sqlite3-path STR [sqlite3]
            Path to sqlite3 tool

    --tabix-path STR [tabix]
            Path to tabix tool

    --keep-files [OFF]
            Option to keep intermediate files

    --cal-cov [OFF]
            Calculate covariate matrix
```

## Citing SWAM

* Our paper is currently in pre-print and can be found at:

 **Liu A., Kang H.M., Meta-imputation of transcriptome from genotypes across multiple datasets using summary-level data. bioRxiv, 442575.** [Link](https://www.biorxiv.org/content/10.1101/2021.05.04.442575v1)


