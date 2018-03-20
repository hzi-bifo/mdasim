MDAsim 2 extends the original published MDAsim 1.2 to include single nucleotide copy errors and outputs a respective log file of the introduced errors.
The original MDAsim 1.2 can be [found on sourceforge](https://sourceforge.net/projects/mdasim/) and we have kept its [README_mdasim1-2.txt](README_mdasim1-2.txt) for reference.

For credits for different features, please refer to the [CREDITS_mdasim.txt](CREDITS_mdasim.txt) and the [change logs of the releases](https://github.com/hzi-bifo/mdasim/releases).
The license, as set by the original MDAsim 1.2, is provided in [LICENSE.txt](LICENSE.txt).

# SUMMARY

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
mdasim --input=Staphylococcus_aureus_USA300_FPR3757.fa --primers=primersList.fasta --coverage=4 --output=test_mdasim_out_prefix_ >test_mdasim.log
```

# Citation

MDAsim 2 extends the original MDAsim 1.2. Whenever you use MDAsim 2, please cite both versions:
* MDAsim 1.2 citation: [Taghavi, Z. and Draghici, S. (2012). MDAsim: a multiple displacement amplification simulator. In IEEE Conference on Bioinformatics and Biomedicine, pages 575–578.](https://doi.org/10.1109/BIBM.2012.6392622).
* MDAsim 2 citation: Until there is a publication to cite, please cite the link to the tagged version that you use (e.g.: <https://github.com/hzi-bifo/mdasim/releases/v2.0.0>) or the exact bioconda version of the package that you use.

