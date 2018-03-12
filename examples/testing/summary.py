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
currentBases = {};

sys.stdout.write('Analysing pileup');

for line in inputfile:
    #print(line);
    contents = line.split();
    position = int(contents[1]);
    bases = contents[4];
    isQualityChar = False;
    for base in bases:
        if(base == '^'):                                                        #beginning of a read with quality char following the '^'
            isQualityChar = True;                                               #flag next char as no base
        else:                                                                   #if this char does not mark the beginning of a read
            if(isQualityChar):                                                  #check if its predecessor did
                isQualityChar = False;                                          #in this case, the next char will be a proper base
            else:                                                               #if the char is a normal base
                base = base.upper();                                            #cast letter to upper case
                if(base.isalpha()):
                    if(base in currentBases):
                        currentBases[base] = currentBases[base]+1;                                           #increase counter
                    else:
                        currentBases[base] = 1;
    if(len(currentBases) > 1):
        ref_base = ' ';
        ref_count = 0;
        for base in currentBases:                                               #find reference base (with highest count)
            if(ref_count < currentBases[base]):
                ref_base = base;
                ref_count = currentBases[base];

        for base in currentBases:
            if(base != ref_base):
                allSNPs.append({'position':position, 'ref_base':ref_base, 'alt_base':base, 'ref_count':ref_count, 'alt_count':currentBases[base]});

        sys.stdout.write('.');
        sys.stdout.flush();

    currentBases.clear();

sys.stdout.write(' Printing summary to file\n');
sys.stdout.flush();
allSNPs.sort(key=lambda x: x['position'], reverse=False);

outfile = open(file_out, 'wt');
outfile.write("#pos\tref\talt\tr_cnt\ta_cnt\n");
for snp in allSNPs:
    outfile.write("" + str(snp['position']) + "\t" + snp['ref_base'] + "\t" + snp['alt_base'] + "\t" + str(snp['ref_count']) + "\t" + str(snp['alt_count']) + "\n");

outfile.close();
