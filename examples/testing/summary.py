"""
Reads a pileup-file from Samtools and fetches all SNPs


Format of a line in the pileup:

chromosome:GRCh38:17:0:270000:1	2	N	7	AAaaaaa	~~~~~~~	2,2,2,2,2,2,2
<sequencename> <position> <?> <number_of_matched_sequences> <matched_chars> <?> <?>
"""

import argparse
import os
import sys

# parse cmd-line arguments
parser = argparse.ArgumentParser()

#get cmd-line arguments
parser.add_argument('--input', action='store', help='Samtools pileup output-file name');
parser.add_argument('--output', action='store', help='Filename of summary-file');
args = parser.parse_args();

# handle files
file_in = args.input;
file_out = args.output;

if not os.path.exists(file_in):
    os.makedirs(file_in);

inputfile = open(file_in, 'rt');
allSNPs = [];

sys.stdout.write('Analysing pileup');

for line in inputfile:
    #print(line);
    contents = line.split();
    position = contents[1];
    bases = contents[4];
    ref_base = '';
    alt_base = '';
    count = 0;
    for base in bases:
        base = base.upper();
        if(base.isalpha()):
            if(ref_base is ''):
                ref_base = base;
            elif (ref_base != base and ref_base is not '' and base):            # TODO: what happens to those GggTttGggGgGGgggg lines? How are those counted?
                if(alt_base is ''):
                    alt_base = base;
                count = count+1;
    if(count > 0):
        allSNPs.append({'position':position, 'ref_base':ref_base, 'alt_base':alt_base, 'count':count});
        sys.stdout.write('.');
        sys.stdout.flush();

sys.stdout.write(' Printing summary to file\n');
sys.stdout.flush();
allSNPs.sort(key=lambda x: x['count'], reverse=True);

outfile = open(file_out, 'wt');
outfile.write("#pos\tref\talt\tcnt\n");
for snp in allSNPs:
    outfile.write("" + str(snp['position']) + "\t" + snp['ref_base'] + "\t" + snp['alt_base'] + "\t" + str(snp['count']) + "\n");

outfile.close();
