/*
Copyright 2012- Zeinab Taghavi (ztaghavi@wayne.edu)
Copyright 2018- Victoria Sack (victoria.sack@helmholtz-hzi.de) and David Lähnemann (david.laehnemann@helmholtz-hzi.de)

    This file is part of MDAsim.

    MDAsim is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    MDAsim is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MDAsim; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/***************************************************************************
 * Title:          mdasim.C
 * Author:         Zeinab Taghavi
 * Created:        2012
 * Last modified:  5/27/2014
 * Author:         Victoria Sack
 * Last modified:  see https://github.com/hzi-bifo/mdasim
 *
 * Copyright (c) 2012- Zeinab Taghavi
 * Copyright (c) 2018- Victoria Sack
 * All Rights Reserved
 * See file LICENSE for details.
 * Major code changes from 1.2 to 2.0 are marked with #2.0
 ***************************************************************************/
// -a can be non-integer and any positive number
//add save and load fragment instead of fragmentasInput

#include <malloc.h>
#include <unistd.h>

#include "getopt.C"
#include "mdasim.h"
#include "commonMDA.h"
#include "commonMDA.C"
#include "sequence.C"

#define MAIN_COMMAND (char *)"mdasim"
#define VERSION (char *)PACKAGE_STRING
#define MAXFILECHAR	2000

Option OPTIONS[] = {
    Option('l', (char *)"log", NEEDS_ARG, (char *)"          = file name for a log file of all single nucleotide errors that happen during amplification"), //#2.0
    Option('m', (char *)"mutationrate", NEEDS_ARG, (char *)" = chance of a nucleotide substitution"),       //#2.0
	Option('V', (char *)"version", NO_ARG, (char *)"        prints the version"),
	Option('h', (char *)"help", NO_ARG, (char *)"           shows this help"),
	Option('v', (char *)"verbose", NO_ARG, (char *)"        extended verbose for debug mode"),
	Option('I', (char *)"input", NEEDS_ARG, (char *)"        = file name of reference DNA sequence (default: reference.fasta)"),
	Option('O', (char *)"output", NEEDS_ARG, (char *)"       = output files prefix , `Amplicons.fasta` will be appended to the prefix (default: out)"),
	Option('o', (char *)"outputfragments", NO_ARG, (char *)"writes the lists of fragments and primer positions at the end of simulation in two txt files suffixed by Fragments.txt and PrimerPositions.txt"),
	Option('P', (char *)"primers", NEEDS_ARG, (char *)"      = file name of input primers in fasta format (default: primerList.fasta)"),
	Option('p', (char *)"primerNo", NEEDS_ARG, (char *)"     = average number of initial available primers (default: input reference length * coverage / frgLngth * 1000)"),
	Option('L', (char *)"frgLngth", NEEDS_ARG, (char *)"     = average number of synthesized bases per phi29 (default: 70,000 nt; synthesized bases per phi29 has uniform distribution; variance = frgLngth^2 / 1200)"),
	Option('C', (char *)"coverage", NEEDS_ARG, (char *)"     = expected average coverage (default: 1000)"),
	Option('s', (char *)"stepSize", NEEDS_ARG, (char *)"     = number of synthesized bases per phi29 in each step (default: 10000)"),
	Option('A', (char *)"alpha", NEEDS_ARG, (char *)"        = normalized number of primers attached in each step (default: 0.5e-11)"),
	Option('a', (char *)"attachNum", NEEDS_ARG, (char *)"    = number of primers attached per single strand of reference sequence in the first step. It can be any positive number. (overrides -A; alpha = attachNum / (input reference length * primerNo))"),
	Option('R', (char *)"readLength", NEEDS_ARG, (char *)"   = minimum length of output amplicons (default: 10)"),
	Option('S', (char *)"single", NO_ARG, (char *)"         Input reference is amplified as a single strand sequence"),
	//Option('f', (char *)"fragmentasInput", NO_ARG, (char *)"if the input is a fragment list"),
	Option(0, NULL, 0, NULL)
};

const char* errorLogFileName = "";        //#2.0
double mutationRate = 0.00000295;         //#2.0
Coord dnaLength = 0;
Coord frgAveLength = 70000;
Coverage aveCoverage = 1000;
Coverage avePrimerNum = 0;
Coord minimumReadLength = 10;
Coverage phi29StepSize = 10000;  // check if the Coverage is the proper type
double phi29AttachmentNum = 0;
double alpha = 0.5e-11;
bool FragmentasInput = false;
bool Writefragmentasoutput = false;
bool verbose = false;
PrimerList primerList;
DNAType originalDNA;
FragmentList fragmentList;
FreePrimerPositionSets freePrimerPositionSets;
ReadsList readsList;
Coverage singleStrandCoverage;
Coverage wholeCoverage;
Phi29List phi29List;    //list of copying processes
Coverage freePrimerPositionTotal = 0;
Coverage primerCurrentNoTotal = 0;
bool doubleStranded = true;

//****
//#2.0
FILE *errorLog;
bool printLog = false;

/**
 * @brief fprintfSeq prints a DNAType to a file base by base
 * @param filename file to write to
 * @param format for example "%c" for chars (=bases)
 * @param seq DNAType to print to file
 */
void fprintfSeq(FILE * filename, const char * format, DNAType seq)
{
	for (Coord i = 0; i < seq.size(); i ++)
		fprintf (filename, format, seq.at(i));
}

/**
 * @brief writeSeq print DNAType to std-output base by base
 * @param seq DNAType to print
 */
void writeSeq(DNAType seq)
{
	for (Coord i = 0; i < seq.size(); i ++)
		cout << seq.at(i);
}

/**
 * @brief reverseComplement
 * @param base
 * @return the complementary base (A<->T, C<->G, U->A)
 */
inline Base reverseComplement(Base base)
{
	switch(toupper(base))
	{
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T':
		case 'U': return 'A';
		default: return 'T';
	}
}

/** Inline method that has been added to the software of MDAsim 1.2
 * to provide occasional single nucleotide copyerrors for version 2.0
 * #2.0
 * @brief mutateBase
 * @param base      original base
 * @return          another base than the input base, with a chance of 1/3rd per possible alternative
 */
inline Base mutateBase(Base base)
{
        double r = (double)rand()/(double)(RAND_MAX);
        switch(toupper(base))
        {
                case 'C':   if(r <= (1.0/3.0))
                            {
                                return 'G';
                            }
                            else if(r > (1.0/3.0) && r <= (2.0/3.0))
                            {
                                return 'T';
                            }
                            else
                            {
                                return 'A';
                            }
                case 'G':   if(r <= (1.0/3.0))
                            {
                                return 'C';
                            }
                            else if(r > (1.0/3.0) && r <= (2.0/3.0))
                            {
                                return 'T';
                            }
                            else
                            {
                                return 'A';
                            }
                case 'T':   if(r <= (1.0/3.0))
                            {
                                return 'G';
                            }
                            else if(r > (1.0/3.0) && r <= (2.0/3.0))
                            {
                                return 'C';
                            }
                            else
                            {
                                return 'A';
                            }
                case 'U':   if(r <= (1.0/3.0))
                            {
                                return 'G';
                            }
                            else if(r > (1.0/3.0) && r <= (2.0/3.0))
                            {
                                return 'C';
                            }
                            else
                            {
                                return 'A';
                            }
                case 'A':   if(r <= (1.0/3.0))
                            {
                                return 'G';
                            }
                            else if(r > (1.0/3.0) && r <= (2.0/3.0))
                            {
                                return 'C';
                            }
                            else
                            {
                                return 'T';
                            }
                default: return 'T';
        }
}

/**
 * @brief reverseComplementSeq
 * @param seq sequence to create complement for
 * @return complementary sequence
 */
inline DNAType reverseComplementSeq(DNAType seq)
{
	DNAType seqRC;
	for (Coord i = seq.size() - 1; i >= 0; i--)
		seqRC.push_back(reverseComplement(seq.at(i)));
	return seqRC;
}

/**
 * TODO check if this can be just replaced at the neccessary spots instead of using a method for it
 * @brief close_file
 * @param stream
 * @return
 */
int close_file( FILE * stream)
{
	return (fclose(stream));
}

void loadOriginalSequence(string inputFileName)

{
	FILE * dnaFile = NULL;
	dnaFile = open_file(inputFileName.c_str(), "rt");
        //**********************************************************************
        //#2.0
        //bug fix: 1.2 does not load more than one sequence. Since it's better to
        //let the user know that only the first sequence will be amplified, feeding
        //mdasim 2.0+ an input file with more than one sequence will lead to an error
        //
        bool moreseq = true;
        bool seqLoaded = false;
        while(moreseq) {
            Sequence *originalDNAseq = new Sequence(dnaFile);
            //if a sequence can be loaded from the file
            if (originalDNAseq->isLoaded())
            {
                //and if there has not been a sequence loaded yet
                if(!seqLoaded)
                {
                    //read and store that sequence
                    seqLoaded = true;
                    char *DNAseq =  originalDNAseq->getString();
                    for (dnaLength = 0 ; DNAseq[dnaLength] !=0; dnaLength++)
                            originalDNA.push_back(Base(DNAseq[dnaLength]));
                } else
                //else, if there was already a sequence, exit with an error!
                {
                    exitMsg((char *)"Error: MDAsim cannot process more than one input sequence.", INPUT_ARG_ERROR);
                }
            } else
            //meanwhile, if no sequence could be loaded
            {
                    //and if there hasn't been a sequence yet
                    if(!seqLoaded)
                        //exit with error
                        exitMsg((char *)"Error: Input sequence cannot be loaded.", INPUT_ARG_ERROR);
                    else
                        //otherwise, exit loop
                        moreseq = false;

            }
            delete originalDNAseq;
        }
        //**************************************************************************

	if(dnaFile)
		close_file(dnaFile);
}

/**
 * @brief loadPrimers reads a file of primers, parses them and pushes them into the primerList
 * @param inputPrimerFileName   fasta-file containing primer sequences
 */
void loadPrimers(string inputPrimerFileName)
{
    // open the primer file
	FILE * primerFile = NULL;
	primerFile = open_file(inputPrimerFileName.c_str(), "rt");
    // step through the primers
	bool moreseq = true;
	while (moreseq)
	{
        // for each primer, create a Sequence object
		Sequence *primerseq = new Sequence(primerFile);
		if (primerseq->isLoaded())
		{
            // if the primer could be loaded, parse the primer to the primerList
			char *primerseqString =  primerseq->getString();
			DNAType primerTemp;
			for (int i = 0 ; primerseqString[i] !=0; i++)
				primerTemp.push_back(Base(primerseqString[i]));
			primerList.push_back(primerTemp);
		}
		else
            // if no more primer could be loaded, exit the while-loop
			moreseq = false;
		delete primerseq;
	}

	if (primerList.size() == 0)
		exitMsg((char *)"Error: Input primer list cannot be loaded.", INPUT_ARG_ERROR);
	if(primerFile)
		close_file(primerFile);
}

/**
 * @brief addToPrimerPositionSets checks whether a position is valid as a primer position
 * (not occupied, index great enough to put a primer in front of it on the respective fragment)
 * finds the primer that fits to it and stores the pair in freePrimerPositionSets
 * @param pos   Position to check if suitable as a Primer position
 */
void addToPrimerPositionSets(Position pos)
{
    // get the length of the primers, assuming that all primers have equal lengths
    Coord primerLength = primerList.at(0).size();
    // if the position is not behind the end of the primer, return, since
    // the primer should produce a fragment of size at least primerLength + 1
    if (pos.pos < primerLength)
        return;

	bool occupied = false;
	Primer primerTemp;

    // starting from the position going backwards as long as the index lies within the potential primer ending at pos
    // and while the index is not negative
    for (Coord posIndex = pos.pos; (posIndex >= pos.pos - primerLength + 1) && (posIndex >= 0) ; posIndex--)
	{
        // append the reverse complement of the base at the passed position to the primer
        primerTemp.push_back(reverseComplement(fragmentList.at(pos.fragmentNo1 - 1).dna.at(posIndex).base));
        // check if the position is occupied
        occupied = occupied || (fragmentList.at(pos.fragmentNo1 - 1).dna.at(posIndex).occupancy.fragmentNo1 != 0);
	};

    //if all positions were free
	if (!occupied)
	{
        // find the primer fitting at this position
        FreePrimerPositionSets::iterator it = freePrimerPositionSets.find(primerTemp);
        if (it != freePrimerPositionSets.end())
		{
            // add the position to the primer found in the list
            pair<PositionSet::iterator, bool> ret;
			ret = it->second.set.insert(pos);
			if (ret.second == false)
				exitMsg((char *)"Error: primer position duplication in the function addToPrimerPositionSets.", INTERNAL_WOW_ERROR);
			else
				freePrimerPositionTotal ++;
		}
	}
}

/**
 * @brief deleteFromPrimerPositionSets
 * @param pos   Position for which to find the primer to be removed from set of free primers
 */
void deleteFromPrimerPositionSets(Position pos)
{
    //assuming that all primers have equal lengths
    Coord primerLength = primerList.at(0).size();
	bool occupied = false;
	Primer primerTemp;
    // for the length of the primer at the passed position ...
    for (Coord posIndex = pos.pos; (posIndex >= pos.pos - primerLength + 1) && (posIndex >= 0) ; posIndex--)
	{
        // construct a temporary primer according to the given position and check if the fragment is occupied at that position
        primerTemp.push_back(reverseComplement(fragmentList.at(pos.fragmentNo1 - 1).dna.at(posIndex).base));
        occupied = occupied || (fragmentList.at(pos.fragmentNo1 - 1).dna.at(posIndex).occupancy.fragmentNo1 != 0);
	}
	if (occupied)
	{
        // if the position is occupied, find the newly constructed primer in the set of free primers
		FreePrimerPositionSets::iterator it = freePrimerPositionSets.find(primerTemp);
        if (it != freePrimerPositionSets.end())
		{
            // if the primer could be found in the set of free primers, erase the primer from that set
            size_t erasecount;
			erasecount = it->second.set.erase(pos);
			freePrimerPositionTotal = freePrimerPositionTotal - erasecount;
			if (erasecount > 1)
                // Debug output
				cout << "erasecount = " << erasecount << endl;
		}
	}
}

/**
 * @brief initializePrimerPosition creates a FreePrimerPosition for each primer, that means,
 * pair a Primer with a Position on a Fragment and add the pair to the freePrimerPositionSets
 */
void initializePrimerPosition()
{
    Coord averageCountPerPrimer = avePrimerNum / primerList.size();
	Coord primerLength = primerList.at(0).size(); //it is assumed that the length of all primers are equal
	if (avePrimerNum % primerList.size() != 0)
        averageCountPerPrimer ++;
	cout << "\nPrimers length = " << primerLength << endl << "Number of primers = " << (int) primerList.size() << endl << endl;

    // step through all primers
	for (int i = 0; i < primerList.size(); i++)
	{
        // for each primer, create a FreePrimerPosition
		FreePrimerPosition primerPositionTemp;
        // calculate the coverage for the primer by using the average count +/- a random factor
        primerPositionTemp.currentNo = averageCountPerPrimer + (rand() % ( averageCountPerPrimer /10)) - averageCountPerPrimer /20;
        // add the coverage of the current primer to the total value
		primerCurrentNoTotal = primerCurrentNoTotal + primerPositionTemp.currentNo;
        // add the FreePrimerPosition and the current primer as a PrimerPositionPair to the free primer position sets
		pair<FreePrimerPositionSets::iterator, bool> ret;
		ret = freePrimerPositionSets.insert(PrimerPositionPair(primerList.at(i), primerPositionTemp));
        // if the pair could not be inserted into the list, an error occured (duplicate primers)
		if (ret.second == false)
			exitMsg((char *)" error in initializing primer position list. There are dublicate primers in the input primer list.", INTERNAL_WOW_ERROR);
	}

    // for all fragments in the fragmentList
    for (FragmentID fragIndex = 0; fragIndex < fragmentList.size(); fragIndex++)
	{
        //step through the positions on the DNA from end of primer to the end of the dna
        for (Coord fragCoord = primerLength - 1 ; fragCoord < fragmentList.at(fragIndex).dna.size(); fragCoord++)
		{
            // create a position
			Position pos;
            // set the fragment number according to the fragment index in the list
            pos.fragmentNo1 = fragIndex + 1;
            // set the position on the fragment
            pos.pos = fragCoord;
            // check if the position works as a primer position and if so, add it
			addToPrimerPositionSets(pos);
		}
	}
}

/**
 * @brief savePrimerPositionSets TODO
 * @param Filename
 */
void savePrimerPositionSets(string Filename) // add error message
{
	FILE * primerPosFile = NULL;
	primerPosFile = open_file(Filename.c_str(), "wt");
	Coverage counter = 0;
	fprintf (primerPosFile,">> Total available primers =  %ld | Total available positions = %ld \n", primerCurrentNoTotal, freePrimerPositionTotal);
	for (FreePrimerPositionSets::iterator it = freePrimerPositionSets.begin(); it != freePrimerPositionSets.end(); it++)
	{
		fprintf (primerPosFile, "> %ld | primer = ",counter);
		fprintfSeq(primerPosFile, "%c", it->first);
		fprintf (primerPosFile, " | RC = ");
		fprintfSeq(primerPosFile, "%c", reverseComplementSeq(it->first));
		fprintf (primerPosFile, " | ");
		fprintf (primerPosFile, "available primers = %ld | ",it->second.currentNo);
		fprintf (primerPosFile, "available free positions = %ld | ",it->second.set.size());
		fprintf (primerPosFile, "\n");
		for (PositionSet::iterator setit = it->second.set.begin(); setit != it->second.set.end(); setit++)
			fprintf (primerPosFile, "(%ld, %ld), ",setit->fragmentNo1, setit->pos);
		fprintf (primerPosFile, "\n");
		counter ++;
	};
	if(primerPosFile)
		close_file(primerPosFile);
}

/**
 * @brief writePrimerPositionSets
 * @deprecated unused
 */
void writePrimerPositionSets() // add error message
{
	Coverage counter = 0;
	cout << ">> Total available primers =  " << primerCurrentNoTotal <<"| Total available positions = " << freePrimerPositionTotal << endl;
	for (FreePrimerPositionSets::iterator it = freePrimerPositionSets.begin(); it != freePrimerPositionSets.end(); it++)
	{
		cout << "> " << counter << " | org = ";
		writeSeq(it->first);
		cout << " | RC = ";
		writeSeq(reverseComplementSeq(it->first));
		cout << " | ";
		cout <<"available primers = " << it->second.currentNo << " | ";
		cout <<"available free positions = " << it->second.set.size() << " | "<<endl;
		for (PositionSet::iterator setit = it->second.set.begin(); setit != it->second.set.end(); setit++)
			cout << "(" << setit->fragmentNo1 << ", " << setit->pos << ")";
		cout << endl;
		counter ++;
	};
}

/**
 * @deprecated Unused
 */
void writephi29List()
{
	FragmentID phicounter = 0;
	for (list<Phi29>::iterator it = phi29List.begin(); it != phi29List.end(); it++)
	{
		cout << phicounter << "> FragmentID  = " << it->fragmentNo1 << ", expectedLength = " << it->expectedLength << ",  currentPosition = (" << it->currentPosition.fragmentNo1 << ", " << it->currentPosition.pos << ")" << endl;
		phicounter ++;
	}
}

void saveFragmentList(string Filename) // add error message
{
	FILE * fragmentFile = NULL;
	fragmentFile = open_file(Filename.c_str(), "wt");
	for (FragmentID k = 0; k < fragmentList.size(); k++)
	{
		fprintf (fragmentFile, ">%ld\n",k);
		for (int i = 0; i < fragmentList.at(k).dna.size(); i++)
		{
			fprintf (fragmentFile, "(%c, %ld, %ld), ",fragmentList.at(k).dna.at(i).base, fragmentList.at(k).dna.at(i).occupancy.fragmentNo1, fragmentList.at(k).dna.at(i).occupancy.pos);
		};
		fprintf (fragmentFile, "\n");
	};
	if(fragmentFile)
		close_file(fragmentFile);
}

/**
 * @deprecated unused
 */
void writeFragmentList() // add error message
{
	for (FragmentID k = 0; k < fragmentList.size(); k++)
	{
		cout << ">" << k <<endl;
		for (int i = 0; i < fragmentList.at(k).dna.size(); i++)
		{
            if(i <= 8 || i >= fragmentList.at(k).dna.size() - 8)
                cout << "(" << i << ": " <<fragmentList.at(k).dna.at(i).base << ", " << fragmentList.at(k).dna.at(i).occupancy.fragmentNo1 <<", " << fragmentList.at(k).dna.at(i).occupancy.pos << ", " << fragmentList.at(k).dna.at(i).occupancy.original << "), ";
		};
		cout <<endl;
	};
}

/**
 * @deprecated unused
 */
void writeDNA()
{
	cout << "Length of loaded DNA sequence: " << dnaLength << endl;
	for (int i = 0; i < originalDNA.size(); i++)
		cout << originalDNA.at(i);
	cout << endl;
}

/**
 * @deprecated unused
 */
void writePrimers()
{
	cout << "Number of primers: " << primerList.size() << endl;
	for (int k = 0; k < primerList.size(); k++)
	{
		cout << "> " << k << ": ";
		for (int i = 0; i < primerList.at(k).size(); i++)
			cout << primerList.at(k).at(i);
		cout << endl;
	}
}

void initializeFragmentList()
{
	Fragment dnaTemp;
	Position occ;
	occ.fragmentNo1 = 0;
	occ.pos = 0;
	FragmentBase fbaseTemp;
	fbaseTemp.occupancy = occ;
 	for (int i = 0; i < originalDNA.size(); i++)
	{
		fbaseTemp.base = originalDNA.at(i);
                fbaseTemp.occupancy.original = i+1; //#2.0
		dnaTemp.dna.push_back(fbaseTemp);
	};
	fragmentList.push_back(dnaTemp);
	singleStrandCoverage = originalDNA.size();
	wholeCoverage = singleStrandCoverage;
	if (doubleStranded)
	{
		DNAType originalDNARC = reverseComplementSeq(originalDNA);
 		for (int i = 0; i < originalDNARC.size(); i++)
		{
			fbaseTemp.base = originalDNARC.at(i);
                        fbaseTemp.occupancy.original = (-1) * (originalDNARC.size()-i); //#2.0
			dnaTemp.dna.push_back(fbaseTemp);
		};
		fragmentList.push_back(dnaTemp);
		singleStrandCoverage = singleStrandCoverage + originalDNARC.size();
		wholeCoverage = singleStrandCoverage;
	}
}

void writeAttachmentMeritList(vector<AttachmentMerit> attachmentMeritList)
{
	for (int i = 0; i < attachmentMeritList.size(); i++)
	{
		cout <<"( ";
		writeSeq(attachmentMeritList.at(i).primer);
		cout <<   ", "<< attachmentMeritList.at(i).coverage << "), ";
	}
	cout << endl;

}

void addToPrimerAvailability(Position pos)
{
	Coord primerLength = primerList.at(0).size(); // assuming length of all of the primers are equal
	Coord posBegin, posEnd;
	posBegin = pos.pos + primerLength - 1;
	posEnd = pos.pos - primerLength + 1;
	if (posBegin > fragmentList.at(pos.fragmentNo1 - 1).dna.size() - 1)
		posBegin = fragmentList.at(pos.fragmentNo1 - 1).dna.size() - 1;
	if (posEnd < 1)
		posEnd = 1; // primer should generate a fragment of size at least primerLength +1
	for (Coord i = posBegin; (i >= pos.pos) && ( i >= posEnd + primerLength - 1); i --)
	{
		Position newpos = pos;
		newpos.pos = i;
		addToPrimerPositionSets(newpos);
	}
}

/**
 * @brief deleteFromPrimerAvailability removes a primer from the set of free primers
 * @param pos Position of the primer on a fragment
 */
void deleteFromPrimerAvailability(Position pos)
{
	Coord primerLength = primerList.at(0).size(); /* assuming length of all of the primers are equal */
	Coord posBegin, posEnd;
	posBegin = pos.pos + primerLength - 1;
	posEnd = pos.pos - primerLength + 1;
	if (posBegin > fragmentList.at(pos.fragmentNo1 - 1).dna.size() - 1)
		posBegin = fragmentList.at(pos.fragmentNo1 - 1).dna.size() - 1;
	if (posEnd < 0)
		posEnd = 0;
    // from the end of the primer to its beginning (the termination clause here is double I think, see definition of posEnd) TODO check
	for (Coord i = posBegin; (i >= pos.pos) && ( i >= posEnd + primerLength - 1); i --)
	{
        //create a new position relative to the passed one
		Position newpos = pos;
		newpos.pos = i;
        // delete position from set
		deleteFromPrimerPositionSets(newpos);
	}
}

pair<Primer, Position> findPrimerPositionInSet(Coverage randomPrimerPosition)
{
	FreePrimerPositionSets::iterator primerit = freePrimerPositionSets.begin();
	Coverage counter = 0;
	do
	{
		counter = counter + primerit->second.set.size();
		primerit ++;
	}
	while ((counter <= randomPrimerPosition) && (primerit != freePrimerPositionSets.end()));

	if (counter < randomPrimerPosition)
		exitMsg((char *)" Error in findPrimerPositionInSet. randomPrimerPosition is larger than the total number of available free positions.", INTERNAL_WOW_ERROR);
	else
	{
		primerit --;
		counter = counter - primerit->second.set.size(); // check if primerit->second.currentNo > 0
		if  (primerit == freePrimerPositionSets.end())
			exitMsg((char *)" Error in findPrimerPositionInSet. An error occured related to primerit.", INTERNAL_WOW_ERROR);
		Coord randomIndex = randomPrimerPosition - counter;

		if  ((randomIndex < 0) || (randomIndex >= primerit->second.set.size()))
			exitMsg((char *)" Error in findPrimerPositionInSet. An error occured related to counter.", INTERNAL_WOW_ERROR);
		PositionSet::iterator setit = primerit->second.set.begin();
		for (int i = 0 ; ((i < randomIndex) && (setit != primerit->second.set.end())) ; i++)
			setit ++;
		if (setit == primerit->second.set.end())
			exitMsg((char *)" Error in findPrimerPositionInSet. An error occured related to setit.", INTERNAL_WOW_ERROR);
		Position fragpos = *setit;
		Primer primerTemp = primerit->first;
		primerit->second.currentNo --;
		if (primerit->second.currentNo <= 0)
		{
			freePrimerPositionTotal = freePrimerPositionTotal - primerit->second.set.size();
			if (verbose)
			{
				cout << "randomPrimerPosition = " << randomPrimerPosition << " :: ";
				writeSeq(primerTemp);
				cout << ", freePrimerPositionTotal = " << freePrimerPositionTotal << ", primerit->second.set.size() = " << primerit->second.set.size() << endl;
			};
			freePrimerPositionSets.erase(primerTemp);

		}
		primerCurrentNoTotal --;
		return pair<Primer, Position> (primerTemp, fragpos);
	}
}

/**
 * @brief attachPhi29 attaches primers to the existing fragments and creates a new fragment for each primer
 * @param phi29AttachmentNumCurrent     number of primers to be attached
 */
void attachPhi29(int phi29AttachmentNumCurrent)
{
    // primers that have been attached so far
	Coverage primerIndex = 0;
    //while there have to be more primers added and there are still free primers ...
    while ((primerIndex < phi29AttachmentNumCurrent) && (freePrimerPositionTotal > 0))
	{
        //select a random primer from the set
        Coverage randomPrimerPosition = rand() % freePrimerPositionTotal;
		pair<Primer, Position> primerpositionpair;
		primerpositionpair = findPrimerPositionInSet(randomPrimerPosition);
		Position fragPos = primerpositionpair.second;
		Primer primerAttached = primerpositionpair.first;
		Coord primerLength = primerAttached.size();

        //check if the primer can be attached here (this should not happen since it's checked at initialization of the primerPositionPairs
		if (fragPos.pos < primerLength)
			exitMsg((char *) "Error in attachPhi29. There is a primer position in the list which will generate a fragment of size 0 or less.", INTERNAL_WOW_ERROR);
		else
		{
            //check if the fragment is occupied:
            bool fragmentOccupancy = false;
            //for the length of the primer, check if the fragment is occupied by doing the following:
            for ( Coord k =  fragPos.pos;  k >= fragPos.pos - primerLength + 1; k--)
			{
                // in the fragmentlist, check the fragment that this one came from (?)
                // in the dna of that fragment, get the base at position k  <= fragPos.pos
                // and check if the occupancy has anything other for the attached fragment than 0 (== nothing)
                fragmentOccupancy = fragmentOccupancy  || (fragmentList.at(fragPos.fragmentNo1 - 1).dna.at(k).occupancy.fragmentNo1 != 0);
			};
            // if the fragment is occupied, the primer position is in the list of free positions but it is not free in real, so exit with an error
            if (fragmentOccupancy == true)
				exitMsg((char *) "Error in attachPhi29. The primer position is in the list of free positions but it is not free on the fragment.", INTERNAL_WOW_ERROR);
            else
			{
                // the fragment is free at this position for this primer

                // get the position relative to the reference sequence, #2.0
                Coord original = (-1) * fragmentList.at(fragPos.fragmentNo1 - 1).dna.at(fragPos.pos).occupancy.original;

				primerIndex ++;

                // get the fragment index of end of list
				FragmentID fragmentIndex = FragmentID(fragmentList.size());
                // create a new copying process struct
				Phi29 newPhi29;
                // new fragmentIndex
				newPhi29.fragmentNo1 = fragmentIndex + 1;
                // set a (random) expected length for this fragment
				if ((frgAveLength / 10) == 0)
					newPhi29.expectedLength = frgAveLength;
				else
					newPhi29.expectedLength = frgAveLength + (frgAveLength / 20) - (rand() % (frgAveLength / 10));
                // get the current coordinates of the copying process
                Coord currentCoord = fragPos.pos - primerLength + 1;
                // length of the new fragment can not be more than remaining part of the copied fragment
				if (newPhi29.expectedLength > currentCoord)
					newPhi29.expectedLength = currentCoord;
                // set same position parameters for fragment and copying process
				newPhi29.currentPosition.fragmentNo1 = fragPos.fragmentNo1;
				newPhi29.currentPosition.pos = currentCoord;

                // add copying process to list
				phi29List.push_back(newPhi29);

                // create a new fragment
				Fragment newFragment;
                // for the length of the primer
				for (Coord k = 0 ; k < primerLength ; k ++)
				{
                    // create a base
					FragmentBase baseTemp;
                    // set the base according to the attached primer
					baseTemp.base = primerAttached.at(k);
                    // set the occupancy to the position on the fragment as stored by the primer
					baseTemp.occupancy = fragPos;

                    baseTemp.occupancy.original = original + k; //#2.0

                    // add the new base to the new fragment
					newFragment.dna.push_back(baseTemp);
                    // set the fragment occupancy of the preceding (??) fragment to the new fragment index
					fragmentList.at(fragPos.fragmentNo1 - 1).dna.at(fragPos.pos).occupancy.fragmentNo1 = -(fragmentIndex + 1);
                    // set the position of the base on the preceeding fragment relative to the new fragment (???)
                    fragmentList.at(fragPos.fragmentNo1 - 1).dna.at(fragPos.pos).occupancy.pos = k;
                    // delete the primer from the list of available primers
                    deleteFromPrimerAvailability(fragPos);
                    fragPos.pos --;
				}
                // add the new fragment to list of fragments
				fragmentList.push_back(newFragment);
                // update single strand coverage and whole coverage
				singleStrandCoverage = singleStrandCoverage - primerLength;
				wholeCoverage = wholeCoverage + primerLength;
			}
		}
	}
}

/* a predicate for removing finished phi29 implemented as function */
bool zeroExpectedLength (Phi29 &phi29Value) { return (phi29Value.expectedLength == 0); }

/**
 * @brief OneStepAheadPhi29 TODO
 */
void OneStepAheadPhi29()
{
	for (int k = 0; k < phi29StepSize; k++)
	{
        //this counter just counts the iterations but is never used elsewhere
        int phi29Counter = 0;
        //for all active copying processes ...
		for (list<Phi29>::iterator it = phi29List.begin(); it != phi29List.end(); it++)
		{
            // create a position both for the new and the reference fragment
			Position positionOnNewFrag, positionOnOriginalFrag;
            // set occupancy of reference
			positionOnOriginalFrag.fragmentNo1 = it->currentPosition.fragmentNo1;
            // set position of reference but decrease by one
            positionOnOriginalFrag.pos = it->currentPosition.pos;
			positionOnOriginalFrag.pos--;
			if (positionOnOriginalFrag.pos < 0)
				exitMsg((char *)"Error: error in OneStepAheadPhi29. positionOnOriginalFrag.pos is negative.", INTERNAL_WOW_ERROR);
            // set occupancy of new fragment
			positionOnNewFrag.fragmentNo1 = it->fragmentNo1;
            // set position to end of the fragment
			positionOnNewFrag.pos = fragmentList.at(it->fragmentNo1 - 1).dna.size();
            // create a reference base and a new base
			FragmentBase newBase, originalBase;
			originalBase = fragmentList.at(positionOnOriginalFrag.fragmentNo1 - 1).dna.at(positionOnOriginalFrag.pos);
			newBase.base = reverseComplement(originalBase.base);

            //*******************************************************************
            // Single nucleodite errors are generated here for version 2.0
            // #2.0
            newBase.occupancy.original = (-1) * originalBase.occupancy.original; //#2.0

            double r = (double)rand()/(double)(RAND_MAX);
            if(r <= mutationRate)
            {
                if(printLog) {
                    errorLog = fopen(errorLogFileName, "a+");

                    if(newBase.occupancy.original > 0)
                    {
                        fprintf(errorLog, "%ld\t", ((long int) newBase.occupancy.original)-1);
                        fprintf(errorLog, "%c\t", newBase.base);
                    }
                    else
                    {
                        fprintf(errorLog, "%ld\t", ((long int) ((-1)*newBase.occupancy.original))-1);
                        fprintf(errorLog, "%c\t", reverseComplement(newBase.base));
                    }

                    newBase.base = mutateBase(newBase.base);

                    if(newBase.occupancy.original > 0)
                    {
                        fprintf(errorLog, "%c\n", newBase.base);
                    }
                    else
                    {
                        fprintf(errorLog, "%c\n", reverseComplement(newBase.base));
                    }

                    fclose(errorLog);
                } else {
                    newBase.base = mutateBase(newBase.base);
                }

            }

            //*******************************************************************

			newBase.occupancy.fragmentNo1 = positionOnOriginalFrag.fragmentNo1;
			newBase.occupancy.pos = positionOnOriginalFrag.pos;
			fragmentList.at(positionOnNewFrag.fragmentNo1 - 1).dna.push_back(newBase);
			/* If the original fragment is attached to some other fragments. In other words, if the original fragment is still connected This is the position*/
			if (originalBase.occupancy.fragmentNo1 !=0)
			{
				/* This is the position */
				Position occupancyTemp = originalBase.occupancy;
				if (occupancyTemp.fragmentNo1 < 0 )
					occupancyTemp.fragmentNo1 = - occupancyTemp.fragmentNo1;
				cout.flush();
				FragmentID fragmentNameTemp = fragmentList.at(occupancyTemp.fragmentNo1 - 1).dna.at(occupancyTemp.pos).occupancy.fragmentNo1;
				if (!((positionOnOriginalFrag.fragmentNo1 == fragmentNameTemp) || (positionOnOriginalFrag.fragmentNo1 == - fragmentNameTemp)))
					exitMsg((char *)"Error: fragmentList has inconsistency in occupancy variable.", INTERNAL_WOW_ERROR);
				/*detach fragment occupancyTemp.fragmentNo1 from fragment positionOnOriginalFrag.fragmentNo1 */
				fragmentList.at(occupancyTemp.fragmentNo1 - 1).dna.at(occupancyTemp.pos).occupancy.fragmentNo1 = 0;
				/*one of the double-stranded copy of DNA is released as single-strand*/
				singleStrandCoverage++;
				//add to the primer free availability
				addToPrimerAvailability(occupancyTemp);
			} else
			{
				/*one of the free single-stranded copy of DNA is now occupied*/
				singleStrandCoverage--;
				fragmentList.at(positionOnOriginalFrag.fragmentNo1 - 1).dna.at(positionOnOriginalFrag.pos).occupancy.fragmentNo1 = MaxFragmentID; /* temporary value to make this base occupied, the original value is given 3 lines ahead */
				//delete the primer free availability
				deleteFromPrimerAvailability(positionOnOriginalFrag);
			}
			/* a new base has been generated */
			wholeCoverage++;
			fragmentList.at(positionOnOriginalFrag.fragmentNo1 - 1).dna.at(positionOnOriginalFrag.pos).occupancy.fragmentNo1 = -positionOnNewFrag.fragmentNo1;
			fragmentList.at(positionOnOriginalFrag.fragmentNo1 - 1).dna.at(positionOnOriginalFrag.pos).occupancy.pos = positionOnNewFrag.pos;
			/* update the phi29 position */
			it->currentPosition.pos--;
			it->expectedLength--;
			phi29Counter++;
		}

		phi29List.remove_if(zeroExpectedLength);
	}
}

/**
 * @brief sequenceFragments TODO
 * @param outputFragmentsFile
 */
void sequenceFragments(string outputFragmentsFile)
{
    // coverage from previous round
	Coverage oldCoverage = 0;
    // checks if coverage is changeing at all, if not, the loop will terminate eventually
	Coverage coverageNotChangingCounter = 0;
    // counts while-loops
	Coverage counter = 0;

    // while the desired coverage has not been reached yet and it is still increasing ...
	while ((wholeCoverage < dnaLength * aveCoverage) && (coverageNotChangingCounter < 100))
	{
		counter ++;
		oldCoverage = wholeCoverage;
		Coverage phi29AttachmentNumCurrent = Coverage(phi29AttachmentNum * singleStrandCoverage / dnaLength * primerCurrentNoTotal / avePrimerNum);
		if (phi29AttachmentNumCurrent < 1)
			phi29AttachmentNumCurrent = 1;

        // attach primers to fragments to create new fragments
		attachPhi29(phi29AttachmentNumCurrent);
        // TODO
		OneStepAheadPhi29();

		double wholeCoverageRate = (double) wholeCoverage / (double) dnaLength;
		double singleCoverageRate = (double) singleStrandCoverage / (double) dnaLength;
		wholeCoverageRate = (double) wholeCoverage / (double) dnaLength;
		singleCoverageRate = (double) singleStrandCoverage / (double) dnaLength;
		if (wholeCoverage == oldCoverage)
			coverageNotChangingCounter ++;
		else
			coverageNotChangingCounter = 0;
		if (counter % 20 == 0)
		{
			if (Writefragmentasoutput && (fragmentList.size()  % 20000 == 0))
				saveFragmentList(outputFragmentsFile);
			cout << "Step = " << counter << " > frag No = " << fragmentList.size() << ", cov = " << wholeCoverageRate;
			if (verbose)
				cout << ", single-stranded cov = " << singleCoverageRate << ", newly activated phi29 = " << phi29AttachmentNumCurrent << ", total active phi29  = " << phi29List.size() << ", available primer No = " << primerCurrentNoTotal;
			cout << endl;
		}
	}
	if (Writefragmentasoutput)
	{
		cout << "Saving fragments ... " << endl;
		saveFragmentList(outputFragmentsFile);
		cout << "Saving fragments finished." << endl;
	};
	double wholeCoverageRate = (double) wholeCoverage / (double) dnaLength;
	double singleCoverageRate = (double) singleStrandCoverage / (double) dnaLength;
	wholeCoverageRate = (double) wholeCoverage / (double) dnaLength;
	singleCoverageRate = (double) singleStrandCoverage / (double) dnaLength;
	cout << "Step = " << counter << " > frag No = " << fragmentList.size() << ", cov = " << wholeCoverageRate;
	if (verbose)
		cout << ", single-stranded cov = " << singleCoverageRate << ", total active phi29 = " << phi29List.size() << ", current available primer No = " << primerCurrentNoTotal;
	cout << endl;

}

void cleaveFragments(string filename, double averageReadLength, FragmentID readNumber)
{
	FILE * readFile = NULL;
	readFile = open_file(filename.c_str(), "wt");
	(averageReadLength) = 0;
	(readNumber) = 0;
    Coord refPos = 0;
	Coverage singleStrandCounterAtTheEnd = 0;
	for (FragmentID fragIndex = 0; fragIndex < fragmentList.size(); fragIndex++)
	{
		bool occupied = (fragmentList[fragIndex].dna[0].occupancy.fragmentNo1 !=0);
		for (Coord posIndex = 0; posIndex < fragmentList[fragIndex].dna.size();)
		{
			DNAType readTemp;
			while (occupied && (posIndex < fragmentList[fragIndex].dna.size()))
			{
				readTemp.push_back(fragmentList[fragIndex].dna[posIndex].base);
				posIndex++;
				if (posIndex < fragmentList[fragIndex].dna.size())
				{
					occupied = (fragmentList[fragIndex].dna[posIndex].occupancy.fragmentNo1 !=0);
					if (fragmentList[fragIndex].dna[posIndex].occupancy.fragmentNo1 == 0)
						singleStrandCounterAtTheEnd ++;
				}
			}
			if (readTemp.size() >= minimumReadLength)
			{
				readsList.push_back(readTemp);
				readNumber ++;

                /****************** #2.0 */
                // get the position of this base on the original reference sequence
                Coord original = fragmentList[fragIndex].dna[posIndex].occupancy.original;
                if(original == 0) {
                    // if there is no value for the original position, infer it from the neighbour base
                    Coord prec = fragmentList[fragIndex].dna[posIndex-1].occupancy.original;
                    original = prec + 1;
                }

                char strand = '+';
                if(original > 0) {
                    refPos = original - readTemp.size() - 1;
                } else {
                    refPos = (-1)*original ;
                    strand = '-';
                }
                /*******************/

				averageReadLength = averageReadLength + readTemp.size();
                fprintf(readFile, ">R%ld | length = %ld | ref = %ld | strand = %c \n", readNumber, readTemp.size(), refPos, strand); //#2.0
				fprintfSeq(readFile,"%c",readTemp);
				fprintf(readFile, "\n");
			}
			DNAType readTemp2;
			while (!occupied && (posIndex < fragmentList[fragIndex].dna.size()))
			{
				readTemp2.push_back(fragmentList[fragIndex].dna[posIndex].base);
				posIndex++;
				if (posIndex < fragmentList[fragIndex].dna.size())
				{
					occupied = (fragmentList[fragIndex].dna[posIndex].occupancy.fragmentNo1 !=0);
					if (fragmentList[fragIndex].dna[posIndex].occupancy.fragmentNo1 == 0)
						singleStrandCounterAtTheEnd ++;

				}
			}
			if (readTemp2.size() >= minimumReadLength)
			{
				readsList.push_back(readTemp2);
				cout.flush();
                readNumber ++;
                /****************** #2.0 */
                // get the position of this base on the original reference sequence
                Coord original = fragmentList[fragIndex].dna[posIndex].occupancy.original;
                if(original == 0) {
                    // if there is no value for the original position, infer it from the neighbour base
                    Coord prec = fragmentList[fragIndex].dna[posIndex-1].occupancy.original;
                    original = prec + 1;
                }

                char strand = '+';
                if(original > 0) {
                    refPos = original - readTemp2.size() - 1;
                } else {
                    refPos = (-1)*original ;
                    strand = '-';
                }
                /*******************/

				averageReadLength = averageReadLength + readTemp2.size();
                fprintf(readFile, ">R%ld | length = %ld | ref = %ld | strand = %c \n", readNumber, readTemp2.size(), refPos, strand); //#2.0
                fprintfSeq(readFile,"%c",readTemp2);
				fprintf(readFile, "\n");
			}
		}
	}
	double totalCov = (double) averageReadLength / double (dnaLength);
	(averageReadLength) = (averageReadLength) / double (readNumber);
	cout << "\nAverage length of amplicons = " << averageReadLength << endl << "Total Coverage = " << totalCov << endl;
	if(readFile)
		close_file(readFile);
}

/**
 * @brief printErrorLogHeader opens the error log file and prints header info
 */
void printErrorLogHeader(string inputFileName)
{
    errorLog = fopen(errorLogFileName, "w");
    if(errorLog == NULL)
    {
        string msg ("Error: Failed to open logfile. Please make sure all folders in the path exist.");
        exitMsg((char *) msg.c_str(), INPUT_ARG_ERROR);
    }
    fprintf(errorLog, "#Generating software: %s\n",VERSION);
    fprintf(errorLog, "#Input sequence file: %s\n", inputFileName.c_str());
    fprintf(errorLog, "#Position count starts at: 0\n");
    fprintf(errorLog, "#pos\tref\tsub\n");
    fclose(errorLog);
}

int main(int argc, char *argv[])
{
	GetOpt opts(argc, argv, OPTIONS);

        // default file names
	string outputName ("out");
	string inputFileName ("reference.fasta");
	string inputPrimerFileName ("primerList.fasta");
	string outputFragmentsFile ("outFragments.txt");
	string outputReadsName ("outAmplicons.fasta");
	string outputPrimerPositionsFile ("outPrimerPositions.txt");

        // read and set options
	while (opts.hasNext())
	{
		Option *current = opts.next();
		char count = current->getShortForm();
        if (count == 'V')
        {
            // print software version and exit
			version(VERSION);
			exit(EXIT_SUCCESS);
		}
        else if (count == 'h')
		{
            // print help and exit
			printf("\nUsage: ");
			printf(MAIN_COMMAND);
			printf(" [optional args] --input=<input.fa> --output=<mda-amplified_fasta_prefix> --primers=<primers.fasta>\n");
			printf("\nNote: The above used arguments have defaults, but it is recommended to explicitly set them.\n");
			printf("Note: Arguments that require a value are marked with an '=' sign below. This needs to be used \n");
            printf("      between the argument and the value on the command line.\n");
			printf("%s\n", opts.help());
			exitMsg(NULL, NO_ERROR);
		}
        else if (count == 'l')
        {
            // set log file name
            //#2.0
            errorLogFileName = current->getArg();
            printLog = true;
        }
		else if (count =='v')
			verbose = true;
        /*else if (count == 'f')
            FragmentasInput = true;         //old code: this is not in Options anymore!!*/
		else if (count == 'I')
			inputFileName = current->getArg();
        else if (count == 'O')
            outputName = current->getArg();
        else if (count == 'm')
            mutationRate = strtod(current->getArg(), NULL);          //#2.0
		else if (count == 'o')
			Writefragmentasoutput = true;
		else if (count == 'P')
			inputPrimerFileName = current->getArg();
		else if (count == 'L')
			frgAveLength = strtoul(current->getArg(), NULL, 0);
		else if (count == 'C')
			aveCoverage = strtoul(current->getArg(), NULL, 0); //atoi all ints
		else if (count == 'p')
			avePrimerNum = strtoul(current->getArg(), NULL, 0);
		else if (count == 'R')
			minimumReadLength = strtoul(current->getArg(), NULL, 0);
		else if (count == 's')
			phi29StepSize = strtoul(current->getArg(), NULL, 0);
		else if (count == 'a')
			phi29AttachmentNum = strtod(current->getArg(), NULL);
		else if (count == 'A')
			alpha = strtod(current->getArg(), NULL);
		else if (count == 'S')
			doubleStranded = false;

	}

    // error log file for version 2.0, print header #2.0
    if(printLog)
    {
        printErrorLogHeader(inputFileName);
    }

    // set output filenames
	outputFragmentsFile = outputName;
	outputReadsName = outputName;
    outputFragmentsFile += "Fragments.txt";
	outputReadsName += "Amplicons.fasta";
	outputPrimerPositionsFile = outputName;
	outputPrimerPositionsFile += "PrimerPositions.txt";

    // check if output file can be opened BEFORE costly amplification process starts! #2.0
    FILE * outputCheck = fopen(outputReadsName.c_str(), "w");
    if(outputCheck == NULL)
    {
        string msg ("Error: Failed to open output file. Please make sure all folders in the path exist.");
        exitMsg((char *) msg.c_str(), INPUT_ARG_ERROR);
    }
    else
    {
        fclose(outputCheck);
    }

    // print input arguments (why though?) and software version
	string executedcommand;
	for(int i = 0; i < argc; i++)
	{
		executedcommand += argv[i];
		executedcommand += " ";
	}
	executedcommand += "\n \n";

    time_t now;
    time(&now);

	version(VERSION);
	cout << "**************************************"<< endl;
	cout << executedcommand;

    // print file name and amplification execution infos
	cout << "Output amplicons file name = " << outputReadsName << endl;
	if (Writefragmentasoutput)
	{
		cout << "Output fragments file name = " << outputFragmentsFile << endl << "Output PrimerPositions file name = " << outputPrimerPositionsFile << endl;
	}
	cout << "Average number of synthesized bases per phi29 = " << frgAveLength << " nt\n" << "Average expected coverage = " << aveCoverage << endl << "Minimum length of output amplicons = " << minimumReadLength << endl << "Number of synthesized bases per phi29 in each step = " << phi29StepSize << " nt" << endl;
 	if  (avePrimerNum != 0)
	{
		cout << "Average initial primer counts (primerNo) = " << avePrimerNum << endl;
	}
	cout << endl;
	cout.flush();

    cout << "**************************************" << endl;
    cout << "Loading reference DNA sequence ..." << endl;
    if(inputFileName != "")
        loadOriginalSequence(inputFileName);
    //writeDNA();
    initializeFragmentList();
    if (Writefragmentasoutput)
        saveFragmentList(outputFragmentsFile);
    cout << "DNA sequence loaded." << endl;

    // set average primer number based on dna length, coverage and fragment length
    // CoverageToPrimerConst = 1000 defined in header
    if (avePrimerNum == 0)
    {
		avePrimerNum = dnaLength * aveCoverage / frgAveLength * CoverageToPrimerConst;
		cout << "Average initial primer counts (primerNo) = " << avePrimerNum << endl;
	}

    // check if number of primers attached per single strand of reference sequence in the first step is valid and set it
	if (phi29AttachmentNum != 0)
	{
		if (phi29AttachmentNum < 0)
			exitMsg((char *)"Error: \"a\" should be a positive number.", INTERNAL_WOW_ERROR);
		alpha = (double) phi29AttachmentNum / (double) (dnaLength * avePrimerNum);
	}
	else
	{
        // if the number of primers attached in first step is zero, check normalized number of primers attached in each step
        // and calculate number of primers in first step accordingly
		if (alpha < 0)
			exitMsg((char *)"Error: \"alpha\" should be a positive number.", INTERNAL_WOW_ERROR);
		phi29AttachmentNum = (int) (alpha * dnaLength * avePrimerNum);
	}

	cout << "\nNormalized number of primers attached in each step (alpha) = " << alpha << endl << "Number of primers attached in the first step (attachNum) = " << phi29AttachmentNum << endl << "attachNum = alpha * input DNA length * primerNo\n\n";

    // load primers
	cout << "\nLoading primers ..." << endl;
	if (inputPrimerFileName !="")
		loadPrimers(inputPrimerFileName);
	cout << "Primer list loaded." << endl;

    // for all primers, check where they can go on the fragments
	initializePrimerPosition();

    // if needed write an output file about where the primers have been positioned
	if (Writefragmentasoutput)
		savePrimerPositionSets(outputPrimerPositionsFile);

	cout << "**************************************"<< endl;
	cout << "Amplification started ..." << endl;
	cout << "**************************************"<< endl;

    // for debugging
	if (verbose)
		cout << " RAND_MAX = " << RAND_MAX << endl;

	srand(time(NULL));
    sequenceFragments(outputFragmentsFile);
	if (Writefragmentasoutput)
		savePrimerPositionSets(outputPrimerPositionsFile);
	double averageReadLength;
	FragmentID readNumber;
	cout << "\nCleavage started ..." << endl;
	cout << "Writing amplicons in " << outputReadsName << endl;
    cleaveFragments(outputReadsName, averageReadLength, readNumber);
	cout << "\nCleavage finished." << endl;
	cout << "**************************************"<< endl;
	cout << "Amplification finished." << endl;
	cout << "**************************************"<< endl;
	return 0;
}
