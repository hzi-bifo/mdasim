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
#include "commonMDA.C"
#include "sequence.C"

#define FILE_STRING (char *)"mdasim"
#define FILE_VERSION (char *)"mdasim 2.0.1"
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
Phi29List phi29List;
Coverage freePrimerPositionTotal = 0;
Coverage primerCurrentNoTotal = 0;
bool doubleStranded = true;

//****
//#2.0
FILE *errorLog;
bool printLog = false;

int fprintfSeq(FILE * filename, const char * format, DNAType seq)
{
	for (Coord i = 0; i < seq.size(); i ++)
		fprintf (filename, format, seq.at(i));
}

void writeSeq(DNAType seq)
{
	for (Coord i = 0; i < seq.size(); i ++)
		cout << seq.at(i);
}

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


inline DNAType reverseComplementSeq(DNAType seq)
{
	DNAType seqRC;
	for (Coord i = seq.size() - 1; i >= 0; i--)
		seqRC.push_back(reverseComplement(seq.at(i)));
	return seqRC;
}


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
            if (originalDNAseq->isLoaded())                                     //if a sequence can be loaded from the file
            {
                if(!seqLoaded)                                                  //and if there has not been a sequence loaded yet
                {
                    seqLoaded = true;                                           //read and store that sequence
                    char *DNAseq =  originalDNAseq->getString();
                    for (dnaLength = 0 ; DNAseq[dnaLength] !=0; dnaLength++)
                            originalDNA.push_back(Base(DNAseq[dnaLength]));
                } else {                                                        //else, if there was already a sequence, exit with an error!
                    exitMsg((char *)"Error: MDAsim cannot process more than one input sequence.", INPUT_ARG_ERROR);
                }
            }
            else                                                                //meanwhile, if no sequence could be loaded
            {
                    if(!seqLoaded)                                              //and if there hasn't been a sequence yet
                        exitMsg((char *)"Error: Input sequence cannot be loaded.", INPUT_ARG_ERROR);    //exit with error
                    else
                        moreseq = false;                                        //otherwise, exit loop

            }
            delete originalDNAseq;
        }
        //**************************************************************************

	if(dnaFile)
		close_file(dnaFile);
}

void loadPrimers(string inputPrimerFileName)
{
	FILE * primerFile = NULL;
	primerFile = open_file(inputPrimerFileName.c_str(), "rt");
	bool moreseq = true;
	while (moreseq)
	{
		Sequence *primerseq = new Sequence(primerFile);
		if (primerseq->isLoaded())
		{
			char *primerseqString =  primerseq->getString();
			DNAType primerTemp;
			for (int i = 0 ; primerseqString[i] !=0; i++)
				primerTemp.push_back(Base(primerseqString[i]));
			primerList.push_back(primerTemp);
		}
		else
			moreseq = false;
		delete primerseq;
	}
	if (primerList.size() == 0)
		exitMsg((char *)"Error: Input primer list cannot be loaded.", INPUT_ARG_ERROR);
	if(primerFile)
		close_file(primerFile);
}

void addToPrimerPositionSets(Position pos)
{
	Coord primerLength = primerList.at(0).size(); //assuming that all primers have equal lengths
	bool occupied = false;
	Primer primerTemp;
	if (pos.pos < primerLength) // the primer should produce a fragment of size at least primerLength + 1
		return;
	for (Coord pindex = pos.pos; (pindex >= pos.pos - primerLength + 1) && (pindex >= 0) ; pindex--)
	{
		primerTemp.push_back(reverseComplement(fragmentList.at(pos.fragmentNo1 - 1).dna.at(pindex).base));
		occupied = occupied || (fragmentList.at(pos.fragmentNo1 - 1).dna.at(pindex).occupancy.fragmentNo1 != 0);
	};
	if (!occupied)
	{
		FreePrimerPositionSets::iterator it = freePrimerPositionSets.find(primerTemp);
		if (it == freePrimerPositionSets.end())
		{
		//	writeSeq(primerTemp);
		//	cout << " can not be found in freePrimerPositionList." << endl;
		//	exitMsg((char *)"Error: primer can not be found in freePrimerPositionSets in the function addToPrimerPositionSets.", INTERNAL_WOW_ERROR);
		}
		else
		{
		//	writeSeq(primerTemp);
		//	cout << " is found in primerPositionList." << endl;
			pair<PositionSet::iterator, bool> ret;
			ret = it->second.set.insert(pos);
			if (ret.second == false)
				exitMsg((char *)"Error: primer position duplication in the function addToPrimerPositionSets.", INTERNAL_WOW_ERROR);
			else
				freePrimerPositionTotal ++;
		}
	}
}


void deleteFromPrimerPositionSets(Position pos)
{
	Coord primerLength = primerList.at(0).size(); //assuming that all primers have equal lengths
	bool occupied = false;
	Primer primerTemp;
	for (Coord pindex = pos.pos; (pindex >= pos.pos - primerLength + 1) && (pindex >= 0) ; pindex--)
	{
		primerTemp.push_back(reverseComplement(fragmentList.at(pos.fragmentNo1 - 1).dna.at(pindex).base));
		occupied = occupied || (fragmentList.at(pos.fragmentNo1 - 1).dna.at(pindex).occupancy.fragmentNo1 != 0);
	}
	if (occupied)
	{
		FreePrimerPositionSets::iterator it = freePrimerPositionSets.find(primerTemp);
		if (it == freePrimerPositionSets.end())
		{
		//	writeSeq(primerTemp);
		//	cout << " can not be found in freePrimerPositionList." << endl;
		//	exitMsg((char *)"Error: primer can not be found in freePrimerPositionSets in the function deleteFromPrimerPositionSets.", INTERNAL_WOW_ERROR);
		}
		else
		{
		//	writeSeq(primerTemp);
		//	cout << " is found in primerPositionList." << endl;
			size_t erasecount;
			erasecount = it->second.set.erase(pos);
			freePrimerPositionTotal = freePrimerPositionTotal - erasecount;
			if (erasecount > 1)
				cout << "erasecount = " << erasecount << endl;
		}
	}
}



void initializePrimerPosition() // now it is assumed that the length of all primers are equal
{
	Coord AveCountPerPrimer = primerList.size();
	AveCountPerPrimer = avePrimerNum / AveCountPerPrimer;
	Coord primerLength = primerList.at(0).size(); //it is assumed that the length of all primers are equal
	if (avePrimerNum % primerList.size() != 0)
		AveCountPerPrimer ++;
	cout << "\nPrimers length = " << primerLength << endl << "Number of primers = " << (int) primerList.size() << endl << endl;
	for (int i = 0; i < primerList.size(); i++)
	{
		FreePrimerPosition primerPositionTemp;
		primerPositionTemp.currentNo = AveCountPerPrimer + (rand() % ( AveCountPerPrimer /10)) - AveCountPerPrimer /20;
		primerCurrentNoTotal = primerCurrentNoTotal + primerPositionTemp.currentNo;
		pair<FreePrimerPositionSets::iterator, bool> ret;
		ret = freePrimerPositionSets.insert(PrimerPositionPair(primerList.at(i), primerPositionTemp));
		if (ret.second == false)
			exitMsg((char *)" error in initializing primer position list. There are dublicate primers in the input primer list.", INTERNAL_WOW_ERROR);
	}
	for (FragmentID fn = 0; fn < fragmentList.size(); fn++)
	{
		for (Coord fc = primerLength - 1 ; fc < fragmentList.at(fn).dna.size(); fc++)
		{
			Position pos;
			pos.fragmentNo1 = fn + 1;
			pos.pos = fc;
			addToPrimerPositionSets(pos);
		}
	}
}

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

void writephi29List()
{
	FragmentID phicounter = 0;
	for (list<Phi29>::iterator it = phi29List.begin(); it != phi29List.end(); it++)
	{
		cout << phicounter << "> FragmentID  = " << it->fragmentNo1 << ", expectedLength = " << it->expectedLength << ",  currentPosition = (" << it->currentPosition.fragmentNo1 << ", " << it->currentPosition.pos << ")" << endl;
		phicounter ++;
	}
}


void loadFragmentList(string inputFilename)
{
	FILE * fragmentFile = NULL;
	fragmentFile = open_file(inputFilename.c_str(), "rt");
// fill fragmentList out
//initialize singleStrandCoverage, wholeCoverage; dnaLength
	if(fragmentFile)
		close_file(fragmentFile);
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

void writeFragmentList() // add error message
{
	for (FragmentID k = 0; k < fragmentList.size(); k++)
	{
		cout << ">" << k <<endl;
		for (int i = 0; i < fragmentList.at(k).dna.size(); i++)
		{
			cout << "(" << i << ": " <<fragmentList.at(k).dna.at(i).base << ", " << fragmentList.at(k).dna.at(i).occupancy.fragmentNo1 <<", " << fragmentList.at(k).dna.at(i).occupancy.pos <<"), ";
		};
		cout <<endl;
	};
}


void writeDNA()
{
	cout << "Length of loaded DNA sequence: " << dnaLength << endl;
	for (int i = 0; i < originalDNA.size(); i++)
		cout << originalDNA.at(i);
	cout << endl;
}

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
	for (Coord i = posBegin; (i >= pos.pos) && ( i >= posEnd + primerLength - 1); i --)
	{
		Position newpos = pos;
		newpos.pos = i;
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


void attachPhi29(int phi29AttachmentNumCurrent)
{
	Coverage primerIndex = 0;
        while ((primerIndex < phi29AttachmentNumCurrent) && (freePrimerPositionTotal > 0))          //while there have to be more primers added and there are still free primers
	{
                Coverage randomPrimerPosition = rand() % freePrimerPositionTotal;                   //select a random primer from the set
		pair<Primer, Position> primerpositionpair;
		primerpositionpair = findPrimerPositionInSet(randomPrimerPosition);

		Position fragPos = primerpositionpair.second;
		Primer primerAttached = primerpositionpair.first;
		Coord primerLength = primerAttached.size();
		/* length of Fragment should be larger than a primer */
		if (fragPos.pos < primerLength)
			exitMsg((char *) "Error in attachPhi29. There is a primer position in the list which will generate a fragment of size 0 or less.", INTERNAL_WOW_ERROR);
		else
		{
                        bool fragmentOccupancy = false;                                             //check if the fragment is occupied
                        for ( Coord k =  fragPos.pos;  k >= fragPos.pos - primerLength + 1; k--)    //for the length of the primer, check if the fragment is occupied
			{
                                fragmentOccupancy = fragmentOccupancy  ||
                                        (fragmentList.at(fragPos.fragmentNo1 - 1).dna.at(k).occupancy.fragmentNo1 != 0);
                                // in the fragmentlist
                                // check the fragment that this one came from (?)
                                // and in the dna, get the base at position k = fragPos.pos
                                // and check if the occupancy has anything other for the attached fragment than 0 (== nothing)
			};
			if (fragmentOccupancy == true) /* The primer position is in the list of free positions but it is not free in real */
				exitMsg((char *) "Error in attachPhi29. The primer position is in the list of free positions but it is not free on the fragment.", INTERNAL_WOW_ERROR);
			else /* the fragment is free at this position for this primer */
			{
                                Coord original = fragmentList.at(fragPos.fragmentNo1 - 1).dna.at(fragPos.pos).occupancy.original; //#2.0

				primerIndex ++;
				FragmentID fragmentIndex = FragmentID(fragmentList.size());
				Phi29 newPhi29;
				newPhi29.fragmentNo1 = fragmentIndex + 1;
				if ((frgAveLength / 10) == 0)
					newPhi29.expectedLength = frgAveLength;
				else
					newPhi29.expectedLength = frgAveLength + (frgAveLength / 20) - (rand() % (frgAveLength / 10));
 					Coord currentCoord = fragPos.pos - primerLength + 1;
				/*length of the new fragment can not be more than remaining part of the copied fragment */
				if (newPhi29.expectedLength > currentCoord)
					newPhi29.expectedLength = currentCoord;
				newPhi29.currentPosition.fragmentNo1 = fragPos.fragmentNo1;
				newPhi29.currentPosition.pos = currentCoord;
				phi29List.push_back(newPhi29);
				Fragment newFragment;
				for (Coord k = 0 ; k < primerLength ; k ++)
				{
					FragmentBase baseTemp;
					baseTemp.base = primerAttached.at(k);
					baseTemp.occupancy = fragPos;
                                        baseTemp.occupancy.original = original; //#2.0
					newFragment.dna.push_back(baseTemp);
					fragmentList.at(fragPos.fragmentNo1 - 1).dna.at(fragPos.pos).occupancy.fragmentNo1 = -(fragmentIndex + 1);
                                        fragmentList.at(fragPos.fragmentNo1 - 1).dna.at(fragPos.pos).occupancy.pos = k;
                                        //fragmentList.at(fragPos.fragmentNo1 - 1).dna.at(fragPos.pos).occupancy.original = original; //#2.0
                                        deleteFromPrimerAvailability(fragPos);
					fragPos.pos --;
				}
				fragmentList.push_back(newFragment);
				singleStrandCoverage = singleStrandCoverage - primerLength;
				wholeCoverage = wholeCoverage + primerLength;
			}
		}
	}
}


/* a predicate for removing finished phi29 implemented as function */
bool zeroExpectedLength (Phi29 &phi29Value) { return (phi29Value.expectedLength == 0); }

void OneStepAheadPhi29()
{
	for (int k = 0; k < phi29StepSize; k++)
	{
		int phi29Counter = 0;
		for (list<Phi29>::iterator it = phi29List.begin(); it != phi29List.end(); it++)
		{
			Position positionOnNewFrag, positionOnOriginalFrag;
			positionOnOriginalFrag.fragmentNo1 = it->currentPosition.fragmentNo1;
			positionOnOriginalFrag.pos = it->currentPosition.pos; //
			positionOnOriginalFrag.pos--;
			if (positionOnOriginalFrag.pos < 0)
				exitMsg((char *)"Error: error in OneStepAheadPhi29. positionOnOriginalFrag.pos is negative.", INTERNAL_WOW_ERROR);
			positionOnNewFrag.fragmentNo1 = it->fragmentNo1;
			positionOnNewFrag.pos = fragmentList.at(it->fragmentNo1 - 1).dna.size();
			FragmentBase newBase, originalBase;
			originalBase = fragmentList.at(positionOnOriginalFrag.fragmentNo1 - 1).dna.at(positionOnOriginalFrag.pos);
			newBase.base = reverseComplement(originalBase.base);
                        newBase.occupancy.original = (-1) * originalBase.occupancy.original; //#2.0

                        //*******************************************************************
                        // Single nucleodite errors are generated here for version 2.0
                        // #2.0

                        double r = (double)rand()/(double)(RAND_MAX);
                        if(r <= mutationRate)
                        {
                            if(printLog) {
                                errorLog = fopen(errorLogFileName, "a+");

                                if(newBase.occupancy.original > 0)
                                {
                                    fprintf(errorLog, "%ld\t", ((long int) newBase.occupancy.original));
                                    fprintf(errorLog, "%c\t", newBase.base);
                                }
                                else
                                {
                                    fprintf(errorLog, "%ld\t", ((long int) ((-1)*newBase.occupancy.original)));
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


void sequenceFragments(string outputFragmentsFile, string outputPrimerPositionsFile)
{
	Coverage oldCoverage = 0;
	Coverage coverageNotChangingCounter = 0;
	Coverage counter = 0;
	while ((wholeCoverage < dnaLength * aveCoverage) && (coverageNotChangingCounter < 100))
	{
		counter ++;
		oldCoverage = wholeCoverage;
		Coverage phi29AttachmentNumCurrent = Coverage(phi29AttachmentNum * singleStrandCoverage / dnaLength * primerCurrentNoTotal / avePrimerNum);
		if (phi29AttachmentNumCurrent < 1)
			phi29AttachmentNumCurrent = 1;
		attachPhi29(phi29AttachmentNumCurrent);


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
                                Coord original = fragmentList[fragIndex].dna[posIndex].occupancy.original;
                                Coord prec = fragmentList[fragIndex].dna[posIndex-1].occupancy.original;
                                if(original == 0) {
                                    if(prec > 0) {
                                        original = prec - 1;  //#2.0 hackey hack
                                        cout << "+ " << prec << " -> " << original << endl;
                                    } else if(prec < 0) {
                                        original = prec + 1;  //#2.0 hackey hack
                                        cout << "- " << prec << " -> " << original << endl;
                                    } else {
                                        cout << "There're still gaps" << endl;
                                    }
                                }

                                if(original == 0) {
                                    cout << "I failed" << endl;
                                }

                                char strand = '+'; //#2.0
                                if(original > 0) {
                                    refPos = original - readTemp.size();  //#2.0
                                    if(refPos == 0) cout << "Upsi Daisy +" << endl;
                                } else {
                                    refPos = (-1)*original ; //#2.0
                                    strand = '-';            //#2.0
                                    if(refPos == 0) cout << "Upsi Daisy -" << endl;
                                }
                                /*******************/

				averageReadLength = averageReadLength + readTemp.size();
                                fprintf(readFile, ">R%ld | length = %ld |fragment = %ld | position = %ld | ref = %ld | strand = %c | 1\n", readNumber, readTemp.size(), fragIndex, posIndex, refPos, strand); //#2.0
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
                                Coord original = fragmentList[fragIndex].dna[posIndex].occupancy.original;
                                Coord prec = fragmentList[fragIndex].dna[posIndex-1].occupancy.original;
                                if(original == 0) {
                                    if(prec > 0) {
                                        original = prec - 1;  //#2.0 hackey hack
                                        cout << "+ " << prec << " -> " << original << endl;
                                    } else if(prec < 0) {
                                        original = prec + 1;  //#2.0 hackey hack
                                        cout << "- " << prec << " -> " << original << endl;
                                    } else {
                                        cout << "There're still gaps" << endl;
                                    }
                                }

                                if(original == 0) {
                                    cout << "I failed" << endl;
                                }

                                char strand = '+'; //#2.0
                                if(original > 0) {
                                    refPos = original - readTemp2.size();  //#2.0
                                    if(refPos == 0) cout << "Upsi Daisy +" << endl;
                                } else {
                                    refPos = (-1)*original ; //#2.0
                                    strand = '-';            //#2.0
                                    if(refPos == 0) cout << "Upsi Daisy -" << endl;
                                }
                                /*******************/

				averageReadLength = averageReadLength + readTemp2.size();
                                fprintf(readFile, ">R%ld | length = %ld |fragment = %ld | position = %ld | ref = %ld | strand = %c | 2\n", readNumber, readTemp2.size(), fragIndex, posIndex, refPos, strand); //#2.0
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

int main(int argc, char *argv[])
{
	GetOpt opts(argc, argv, OPTIONS);

	int files = 0;
	string outputName ("out");
	string inputFileName ("reference.fasta");
	string inputPrimerFileName ("primerList.fasta");
	string outputFragmentsFile ("outFragments.txt");
	string outputReadsName ("outAmplicons.fasta");
	string outputPrimerPositionsFile ("outPrimerPositions.txt");

	while (opts.hasNext())
	{
		Option *current = opts.next();
		char count = current->getShortForm();
      		if (count == 'V')
		{
			version(FILE_VERSION);
			exit(EXIT_SUCCESS);
		}
      		else if (count == 'h')
		{
			printf("\nUsage: ");
			printf(FILE_STRING);
			printf(" [optional args] --input=<input.fa> --output=<mda-amplified_fasta_prefix> --primers=<primers.fasta>\n");
			printf("\nNote: The above used arguments have defaults, but it is recommended to explicitly set them.\n");
			printf("Note: Arguments that require a value are marked with an '=' sign below. This needs to be used \n");
                        printf("      between the argument and the value on the command line.\n");
			printf("%s\n", opts.help());
			exitMsg(NULL, NO_ERROR);
		}
                else if (count == 'l') {                                     //#2.0
                        errorLogFileName = current->getArg();
                        printLog = true;
                }
		else if (count =='v')
			verbose = true;
		else if (count == 'f')
			FragmentasInput = true;
		else if (count == 'I')
			inputFileName = current->getArg();
                else if (count == 'O')
                        outputName = current->getArg();
                else if (count == 'm')                                        //#2.0
                        mutationRate = strtod(current->getArg(), NULL);
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

        //*******************************************************************
        // error log file for version 2.0
        // #2.0
        if(printLog) {
            errorLog = fopen(errorLogFileName, "w");
            if(errorLog == NULL) {
                string msg ("Error: Failed to open logfile. Please make sure all folders in the path exist.");
                exitMsg((char *) msg.c_str(), INPUT_ARG_ERROR);
            }
            fprintf(errorLog, "#Generating software: %s\n",FILE_VERSION);
            fprintf(errorLog, "#Input sequence file: %s\n", inputFileName.c_str());
            fprintf(errorLog, "#Position count starts at: 1\n");
            fprintf(errorLog, "#pos\tref\tsub\n");
            fclose(errorLog);
        }
        //*******************************************************************

	outputFragmentsFile = outputName;
	outputReadsName = outputName;
        outputFragmentsFile += "Fragments.txt";
	outputReadsName += "Amplicons.fasta";
	outputPrimerPositionsFile = outputName;
	outputPrimerPositionsFile += "PrimerPositions.txt";

        //*******************************************************************
        // check if output file can be opened BEFORE costly amplification process starts!
        // #2.0
        FILE * outputCheck = fopen(outputReadsName.c_str(), "w");
        if(outputCheck == NULL) {
            string msg ("Error: Failed to open output file. Please make sure all folders in the path exist.");
            exitMsg((char *) msg.c_str(), INPUT_ARG_ERROR);
        } else {
            fclose(outputCheck);
        }
        //*******************************************************************

	string executedcommand;
	for(int i = 0; i < argc; i++)
	{
		executedcommand += argv[i];
		executedcommand += " ";
	}
	executedcommand += "\n \n";

	time_t now;
	time(&now);

	version(FILE_VERSION);
	cout << "**************************************"<< endl;
	cout << executedcommand;
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

	if (!FragmentasInput) // DNA sequence as input
	{
		cout << "**************************************" << endl;
		cout << "Loading reference DNA sequence ..." << endl;
		if(inputFileName != "")
			loadOriginalSequence(inputFileName);
		//writeDNA();
		initializeFragmentList();
		if (Writefragmentasoutput)
			saveFragmentList(outputFragmentsFile);
		cout << "DNA sequence loaded." << endl;

	}
	else
	{
		cout << "**************************************"<< endl;
		cout << "Loading fragments ..." << endl;
		if(inputFileName != "")
			loadFragmentList(inputFileName);
		cout << "Fragments loaded." << endl;
	};
	if (avePrimerNum == 0) // put it after reading the input dna
	{
		avePrimerNum = dnaLength * aveCoverage / frgAveLength * CoverageToPrimerConst;
		cout << "Average initial primer counts (primerNo) = " << avePrimerNum << endl;
	}
	if (phi29AttachmentNum != 0)
	{
		if (phi29AttachmentNum < 0)
			exitMsg((char *)"Error: \"a\" should be a positive number.", INTERNAL_WOW_ERROR);
		alpha = (double) phi29AttachmentNum / (double) (dnaLength * avePrimerNum);
	}
	else
	{
		if (alpha < 0)
			exitMsg((char *)"Error: \"alpha\" should be a positive number.", INTERNAL_WOW_ERROR);
		phi29AttachmentNum = (int) (alpha * dnaLength * avePrimerNum);
	}

	cout << "\nNormalized number of primers attached in each step (alpha) = " << alpha << endl << "Number of primers attached in the first step (attachNum) = " << phi29AttachmentNum << endl << "attachNum = alpha * input DNA length * primerNo\n\n";

	cout << "\nLoading primers ..." << endl;
	if (inputPrimerFileName !="")
		loadPrimers(inputPrimerFileName);
	cout << "Primer list loaded." << endl;
	initializePrimerPosition();
	if (Writefragmentasoutput)
		savePrimerPositionSets(outputPrimerPositionsFile);
	cout << "**************************************"<< endl;
	cout << "Amplification started ..." << endl;
	cout << "**************************************"<< endl;
	if (verbose)
		cout << " RAND_MAX = " << RAND_MAX << endl;
	srand(time(NULL));
	sequenceFragments(outputFragmentsFile, outputPrimerPositionsFile);
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
