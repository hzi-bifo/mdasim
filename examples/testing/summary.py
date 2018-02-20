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
    alt_count = 0;
    ref_count = 0;
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
                    if(ref_base is ''):                                         #set reference base if there is none so far
                        ref_base = base;
                        ref_count += 1;                                         #increase counter
                    elif (ref_base == base):
                        ref_count += 1;
                    else:
                        if(alt_base is ''):                                     #set alternative base if one is found
                            alt_base = base;
                        alt_count +=1                                           #increase counter
    if(alt_count > 0):
        allSNPs.append({'position':position, 'ref_base':ref_base, 'alt_base':alt_base, 'ref_count':ref_count, 'alt_count':alt_count});
        sys.stdout.write('.');
        sys.stdout.flush();

sys.stdout.write(' Printing summary to file\n');
sys.stdout.flush();
allSNPs.sort(key=lambda x: x['alt_count'], reverse=True);

outfile = open(file_out, 'wt');
outfile.write("#pos\tref\talt\tr_cnt\ta_cnt\n");
for snp in allSNPs:
    outfile.write("" + str(snp['position']) + "\t" + snp['ref_base'] + "\t" + snp['alt_base'] + "\t" + str(snp['ref_count']) + "\t" + str(snp['alt_count']) + "\n");

outfile.close();
