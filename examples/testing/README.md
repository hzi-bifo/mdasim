# Testing Pipeline for MDAsim 2

This pipeline was originally developed for testing and debugging the output of MDAsim 2 for the presence of single nucleotide errors and to check the correctness of the respective log file. This is done via alignment of the created amplicons against the original input sequence with [`minimap2`](https://github.com/lh3/minimap2) and identification of differences to that reference with [`samtools mpileup`](http://www.htslib.org/doc/samtools.html). We have found it very useful and might use it for future testing, so we provide a cleaned up and documented version.


## Installation

To use this pipeline, tools will need to be installed via the conda channel [bioconda](https://bioconda.github.io/). For the required one time setup of conda with the bioconda channel, please refer to [bioconda's installation instructions](https://bioconda.github.io/#using-bioconda). Once this is done, all you need to do is issue the following command in the `examples/testing` folder to create the environment necessary for this pipeline:

```
conda create --name mdasim_analysis -c bioconda --file requirements.txt
```


## Usage

### Run MDAsim with logging

Make sure you have installed MDAsim, for instructions see the [Installation section of the main `README.md`](../../README.md#installation).

To obtain amplicons and a respective single nucleotide error log file to test, make sure to run MDAsim with a log file specified. E.g. runs the command specified in the [Usage section of the main `README.md`](../../README.md#usage).

### Set up pipeline input

Enter the `testing` directory and set up the input file section of the `Snakefile`. E.g. if you edit with vim:

```
cd examples/testing
vim Snakefile
```

The Snakefile is set up to work with the output of the usage example from the main `README.md` recommended above. For other inputs and locations, please edit the strings of the variable `ref`, `amp`, `log` and `out_folder`.

### run pipeline

Activate the conda environment with all the dependencies installed:

```
source activate mdasim_analysis
```

To run the pipeline, just issue:
```
snakemake
```

## Output

Unless you edited the `Snakefile`, all intermediate files and results should be written to a new subfolder `output`. With the wildcars as specified in the `Snakefile`, the files contained should be:

### `{outfolder}/{reference}.{amplicons}.aln_mismatches.tsv`

This includes all positions of the input sequence where a single nucleotide mismatch has been found after the MDAsim amplicons were mapped onto the original MDAsim input sequence. The tab-separated format is:
```
#pos	ref	alt	r_cnt	a_cnt
35678	C	T	890	    167
103604	T	C	359	    166
                ...
```

### `{outfolder}/{reference}.{amplicons}.{log}.cnt_comparison.tsv`

This includes two tab-separated tables that compare the single nucleotide substitutions found in the alignment to those logged by MDAsim.

The first one compares the total counts per type of substitution in the following format:
```
#substitution	log_cnt		aln_cnt
A to T 		    42		    46
A to C 		    46		    43
                ...
```

The second one compares every position found in at least one of either the MDAsim single nucleotide error log or the alignment mismatches summary in the following format:
```
#ref_pos  mdasim_ref  mdasim_alt  aln_pos  aln_ref  aln_alt  ref=aln?
36128     C           A           36128    C        A        true 
44775     C           G           44775    C        G        true 
                ...
```

## Credits

The pipeline uses [Snakemake](https://snakemake.readthedocs.io/en/stable/) as the pipelining framework and [bioconda](https://bioconda.github.io/) for software management.
Apart from our own custom python scripts for creating readable summaries, the pipeline makes use of [minimap2](https://github.com/lh3/minimap2) ([Li, H. (2017). Minimap2: fast pairwise alignment for long nucleotide sequences](https://arxiv.org/abs/1708.01492)) and [samtools](http://www.htslib.org/doc/samtools.html).
Exact versions used by the pipeline are specified in [requirements.txt](requirements.txt).
