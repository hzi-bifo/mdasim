# MDAsim 1.3

January 16, 2018
David Lähnemann (david.laehnemann@helmholtz-hzi.de)
Helmholtz Centre for Infection Research
Braunschweig, Germany
----------------------------------------------------------------------------------
In order to incorporate single nucleotide copy errors, the functionality of MDAsim 1.3 adds to that of 1.2. The files that have been altered were the Makefile, mdasim.h (see patch-file), commonmda.c and mdasim.c (see files itself). Currently, the error rate is hard coded 2.95^-6, but we are planning on making it a cmd-line parameter in the future. Other than that, the functionality of MDAsim remains the same. For the original contents of the README of MDAsim1.2, see below.

The original software of 1.2 can be downloaded from [Sourceforge](https://sourceforge.net/projects/mdasim/)

MDAsim has been described in more detail in [Tagliavi, Zeinab, and Sorin Draghici. "MDAsim: A multiple displacement amplification simulator." Bioinformatics and Biomedicine (BIBM), 2012 IEEE International Conference on. IEEE, 2012](https://doi.org/10.1109/BIBM.2012.6392622).


## Usage of MDAsim 1.3
For building MDAsim 1.3, please do
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
The cmd-line arguments work just as in MDAsim 1.2, please refer to the section **D/ Usage**.

## MDAsim 1.3 Change Log
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

----------------------------------------------------------------------------------
----------------------------------------------------------------------------------

## MDAsim 1.2
May 27, 2014
Zeinab Taghavi
ztaghavi@wayne.edu
Wayne State University
Detroit, MI
----------------------------------------------------------------------------------
### SUMMARY
+ A/ Credits
+ B/ Accompaying files
+ C/ Compile
+ D/ Usage
+ E/ Examples
+ F/ Q&A
+ G/ Reference

----------------------------------------------------------------------------------
### A/ Credits

Most part of the code in MDAsim 1.2 was written by Zeinab Taghavi (ztaghavi@wayne.edu). Some part of the code is borrowed from package HyDA (http://sourceforge.net/projects/hyda/).

Many thanks to Hamidreza Chitsaz and Sorin Draghici for helpful discussions.

If you find this tool helpful, please cite the following paper:

Zeinab Taghavi, and Sorin Draghici, (2012). MDAsim: a multiple displacement amplification simulator. In IEEE Conference on Bioinformatics and Biomedicine, pages 575–578.

----------------------------------------------------------------------------------
### B/ Accompaying files

+ Source code of MDAsim 1.2
+ Example files
    + primerList.fasta: list of input primers.
	+ Staphylococcus_aureus_USA300_FPR3757.fa: reference input dna sequence.
	+ mdasim_S_aureus.log: log file of the output of an example command.

----------------------------------------------------------------------------------
### C/ Compile
```
	$ cd mdasim
	$ make clean
	$ make
```
----------------------------------------------------------------------------------

### D/ Usage

Usage: mdasim [options]

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

----------------------------------------------------------------------------------
### E/ Examples

The following command is an example of usage of mdasim
    `$ ./bin/mdasim -I=Staphylococcus_aureus_USA300_FPR3757.fa -C=50 > mdasim_S_aureus.log`

The accompanying example files with the package are
- primerList.fasta: list of input primers.
- Staphylococcus_aureus_USA300_FPR3757.fa: reference input dna sequence.
- mdasim_S_aureus.log: log file of the output of the command.

------------------------------------------------------------------------------------
### F/ Q&A

**Q1** - What is the format of outAmplicons.fasta?
**A1** - The format is fasta. The ID line of each amplicon is in the following format:
```>XX | name = RXX | fragment = XX | position = XX```

The first two values of the ID line for each amplicon is the length of the amplicon and the name of the amplicon. The name of the amplicon is in the format of R(Index of the amplicon). Indexes of amplicons start from 1 and increase for each amplicon by one. The last one shows the total number of amplicons in the file.
The last two values in each line are mostly used for debugging the program. The number that comes after the term `fragment = ` is the ID of the fragment that the respective amplicon comes from. In fact, as it is described in [1], the amplicons are generated by cleaving fragments from the positions the double- and single-starnded regions meet (look at 3rd step of the algorithm in [1]). The number that comes after the term `position = ` is the position of the last nucleotide of the amplicon in the respective fragment.

**Q2** - How can we find the position of the amplicons in the genome?
**A2** - The program does not store the position of the amplicons in the genome. However, to find these values aligner softwares, for example bwa (http://bio-bwa.sourceforge.net/), can be used to align the amplicons to the genome.

**Q3** - How much RAM is needed to run the simulator?
**A3** - The required memory is linearly proportinal to the size of input genome and the final coverage. For example for the S. aureus sample for which the size of the genome is in the order of 3M, to get 50x average coverage 6G RAM is needed. Also, for a genome of size 3G to get average coverage of 50x, approximately 6T of RAM is needed.
One suggestion to reduce the required memory size for large genomes is to break the genome to smaller pieces (may be with some overlaps), then apply MDAsim on each piece, separately.

Please contact me at (ztaghavi@wayne.edu) if you have more questions or suggestions about the software.

------------------------------------------------------------------------------------
### G/ References

[1] Taghavi, Z. and Draghici, S. (2012). MDAsim: a multiple displacement amplification simulator. In IEEE Conference on Bioinformatics and Biomedicine, pages 575–578.
