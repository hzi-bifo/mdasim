# MDAsim 1.3
----------------------------------------------------------------------------------
*In order to incorporate single nucleotide copy errors, the functionality of MDAsim 1.3 adds to that of 1.2. The files that have been altered were the Makefile, mdasim.h (see patch-file), commonmda.c and mdasim.c (see files itself). Currently, the error rate is hard coded 2.95^-6, but we are planning on making it a cmd-line parameter in the future. Other than that, the functionality of MDAsim remains the same.*

## SUMMARY

+ [Credits and License](https://github.com/hzi-bifo/mdasim/blob/master/README.md#credits-and-license)
+ [Usage](https://github.com/hzi-bifo/mdasim/blob/master/README.md#usage)
+ [Examples](https://github.com/hzi-bifo/mdasim#examples)
+ [Accompaying files](https://github.com/hzi-bifo/mdasim/blob/master/README.md#accompaying-files)
+ [Change Log](https://github.com/hzi-bifo/mdasim/blob/master/README.md#change-log)
+ [Reference](https://github.com/hzi-bifo/mdasim#references)

## Credits and License

MDAsim has been described in more detail in [Tagliavi, Zeinab, and Sorin Draghici. "MDAsim: A multiple displacement amplification simulator." Bioinformatics and Biomedicine (BIBM), 2012 IEEE International Conference on. IEEE, 2012](https://doi.org/10.1109/BIBM.2012.6392622).

The original software of 1.2 can be downloaded from [Sourceforge](https://sourceforge.net/projects/mdasim/).

For credits, please refer to the [CREDITS_mdasim.txt](https://github.com/hzi-bifo/mdasim/blob/master/CREDITS_mdasim.txt).
More information about the GNU General Public License can be found in [LICENSE.txt](https://github.com/hzi-bifo/mdasim/blob/master/LICENSE.txt).
This repository also contains the original [README_mdasim1-2.txt](https://github.com/hzi-bifo/mdasim/blob/master/README_mdasim1-2.txt) for your convenience.

## USAGE
To build MDAsim 1.3, please do
```
	$ cd mdasim
	$ make clean
	$ make prefix=./
```
or
```
$ make prefix=Path/to/desired/build/folder
```
Please note that the target folder has to contain (empty) `obj` and `bin`-subfolders as well as the `lib`-folder.
The cmd-line arguments work just as in MDAsim 1.2, please refer to the section [D/ Usage](https://github.com/hzi-bifo/mdasim#d-usage).

To run mdasim, type `mdasim [options]`. For options, see below.
```
	-V,--version      prints the version
	-h,--help      shows this help
	-v,--verbose      extended verbose for debug mode
	-I,--input      =file name of reference DNA sequence (default: reference.fasta)
	-O,--output      =output files prefix (default: out)
	-o,--outputfragments      writes the lists of fragments and primer positions at the end of simulation in two txt files suffixed by Fragments.txt and PrimerPositions.txt
	-P,--primers      =file name of input primers in fasta format (default: primerList.fasta)
	-p,--primerNo      =average number of initial available primers (default: input reference length * coverage / frgLngth * 1000)
	-L,--frgLngth      =average number of synthesized bases per phi29 (default: 70,000 nt; synthesized bases per phi29 has uniform distribution; variance = frgLngth^2 / 1200)
	-C,--coverage      =expected average coverage (default: 1000)
	-s,--stepSize      =number of synthesized bases per phi29 in each step (default: 10000)
	-A,--alpha      =normalized number of primers attached in each step (default: 0.5e-11)
	-a,--attachNum      =number of primers attached in the first step (overrides -A; alpha = attachNum / (input reference length * primerNo))
	-R,--readLength      =minimum length of output amplicons (default: 10)
	-S,--single      Input reference is amplified as a single strand sequence
```

For more detailed description of the parameters refer to [1].

## Example

The following command is an example of usage of mdasim
    `$ ./bin/mdasim -I=Staphylococcus_aureus_USA300_FPR3757.fa -C=50 > mdasim_S_aureus.log`

The accompanying example files with the package are
- primerList.fasta: list of input primers.
- Staphylococcus_aureus_USA300_FPR3757.fa: reference input dna sequence.
- mdasim_S_aureus.log: log file of the output of the command.


## Accompaying files

+ Source code of MDAsim 1.3
+ Example files
    + primerList.fasta: list of input primers.
    + Staphylococcus_aureus_USA300_FPR3757.fa: reference input dna sequence.
    + mdasim_S_aureus.log: log file of the output of an example command.

## Change Log

| DATE         | CHANGES                                                                    | NOTES                 |
| ------------:|:---------------------------------------------------------------------------|:----------------------|
| **01/16/18** | README.md with change log                                                  |                       |
| **01/10/18** | Initialisation of MDAsim 1.3 Repository, applied patch                     |                       |
|              | License and credit correction                                              |                       |
|              | Implemented single nucleotide copy errors with fixed, hardcoded error rate |                       |
| **08/09/17** | Makefile: Differentiation between BASE and TARGET dirs                     | Patch for MDAsim 1.2  |
| **08/08/17** | Makefile: Use prefix in installation                                       | Patch for MDAsim 1.2  |
| **07/19/17** | Bug fix: Include `getopt.h` in `mdasim.h` (for argument parsing)           | Patch for MDAsim 1.2  |
|              | Removed absolute openmpi paths, use any present openmpi                    | Patch for MDAsim 1.2  |
|              | Compile using openmp (mpic++) and pthread                                  | Patch for MDAsim 1.2  |

### References

[1] Taghavi, Z. and Draghici, S. (2012). MDAsim: a multiple displacement amplification simulator. In IEEE Conference on Bioinformatics and Biomedicine, pages 575–578.
