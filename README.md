MDAsim 2 extends the original published MDAsim 1.2, including single nucleotide copy errors and outputs a respective log file of the introduced errors.

The original MDAsim 1.2 can be [found on sourceforge](https://sourceforge.net/projects/mdasim/) and we have kept its [README_mdasim1-2.txt](README_mdasim1-2.txt) for reference.
For credits for different features, please refer to the [CREDITS_mdasim.txt](CREDITS_mdasim.txt) and the [change logs of the releases](https://github.com/hzi-bifo/mdasim/releases).
The license, as set by the original MDAsim 1.2, is provided in [LICENSE.txt](LICENSE.txt).

Information on how to use and cite MDAsim 2:
1. [Installation](#installation)
2. [Usage](#usage)
3. [Citation](#citation)

# Installation

## bioconda

The easiest way to install MDAsim is via [bioconda](https://bioconda.github.io/). For the required one time setup of conda and bioconda, please refer to [bioconda's installation instructions](https://bioconda.github.io/#using-bioconda). Once this is done, all you need to do is issue the following command:

```
conda install mdasim
```

## build from source

To build MDAsim from source, download the latest [tagged release](https://github.com/hzi-bifo/mdasim/releases) and unpack it. Alternatively, you can install the latest version of the repository by cloning it with:

```
git clone https://github.com/hzi-bifo/mdasim.git
```

Enter the created directory:

```
cd mdasim
```

To install in a `bin` folder in this source directory, just run:
```
make
```

To install to a `bin` folder at a custom location, run:
```
make prefix=path/to/desired/build/folder
```

To be able to run the software from anywhere on your system, make sure that the created `bin` directory is in your `$PATH`.

# Usage

Once installed, you can get the full command line usage message with:
```
mdasim --help
```

For an initial quick test run, you can use the provided example files (please note: the `=` between command line arguments and their respective values are required):
```
cd examples
mdasim --input=example_input.fa --primers=primerList.fasta --coverage=15 \
       --output=example_mdasim_out_prefix_ --log=example_mdasim_errors.log >example_mdasim_run.log
```
Note that the file provided under `--input` must contain exactly one sequence, since MDAsim does not process input files with more than one sequence.

## RAM requirements and runtime

The required memory is linearly proportional to the size of the input genome and to the final coverage requested. E.g., for the S. aureus sample for which the size of the genome is in the order of 3M, to get 50x average coverage 6G RAM is needed. So please note that for a genome of size 3G to get average coverage of 50x, approximately 6T of RAM is needed. One suggestion to reduce the required memory size for large genomes is to break the genome to smaller pieces (may be with some overlaps), then apply MDAsim on each piece separately.

The runtime of MDAsim increases exponentially depending on the length of the simulated sequence. Simulating the amplification of 60M bases at a target coverage of 20 takes around 8h.

## Output formats

### `<out_prefix>Amplicons.fasta`

The format is fasta, with the ID line of each amplicon as follows:

```
>R<IOA> | length = <LA> | ref = <POS> | strand = <S>
```
* `R<IOA>`: amplicon name consisting of `R` and an amplicon index counter that starts at 1 (i.e. the last one shows the total number of amplicons in the file). Within the output file, the name of a fragment is unique.
* `<LA>`: length of the amplicon
* `<POS>`: position on the original input sequence where this fragment can be aligned to. Positions start at 0.
* `<S>`: `+` or `-`, indicating the positive or negative strand. In order to align a negative strand with the original input sequence, it must be reverted and complemented.

### `--log errors.log`

The format of the log file for single nucleotide substitution errors is tab separated as follows:
```
#pos\tref\tsub
```

* `pos`: position on the original input sequence (0-based)
* `ref`: reference nucleotide in the original input sequence that was replaced
* `sub`: nucleotide that is generated in the strand of the original input sequence

For consistent reference back to the original input sequence, both `ref` and `sub` will be reported as if incorporated into the input sequence's strand. I.e., if the substitution happens in the complementary strand, both nucleotides will be complemented before logging.

## caveats

### `Floating point exception` on small inputs

Smaller sizes of input fasta sequences combined with a low target coverage can give a `Floating point exception`. In the `example` run provided above, setting `--coverage=10` chowcases this. While we assume numerical issues, we have resisted the urge of [being nerd-sniped](https://www.xkcd.com/356/) in this particular case. But if you want to investigate, please contribute to the [issue we use to track this](https://github.com/hzi-bifo/mdasim/issues/10).

# Citation

MDAsim 2 extends the original MDAsim 1.2. Whenever you use MDAsim 2, please cite both versions:
* MDAsim 1.2 citation: [Taghavi, Z. and Draghici, S. (2012). MDAsim: a multiple displacement amplification simulator. In IEEE Conference on Bioinformatics and Biomedicine, pages 575–578.](https://doi.org/10.1109/BIBM.2012.6392622).
* MDAsim 2 citation: Until there is a publication to cite, please cite the link to the tagged version that you use (e.g.: <https://github.com/hzi-bifo/mdasim/releases/v2.1.1>) or the exact bioconda version of the package that you use.
