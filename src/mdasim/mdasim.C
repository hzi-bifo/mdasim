/*
Copyright 2012- Zeinab Taghavi (ztaghavi@wayne.edu)

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
 * Last modified:  12/21/2012
 *
 * Copyright (c) 2012- Wayne State University
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include <malloc.h>

#include "mdasim.h"
#include "commonMDA.C" // ??
#include "sequence.C"
#define FILE_STRING (char *)"mdasim"
Option OPTIONS[] = {
	Option('V', (char *)"version", NO_ARG, (char *)"prints the version"),
	Option('h', (char *)"help", NO_ARG, (char *)"shows this help"),
	Option('f', (char *)"fragmentasInput", NO_ARG, (char *)"if the input is a fragment list"), 
	Option('I', (char *)"input", NEEDS_ARG, (char *)"=input file (mandatory)"), // default value "dna.fasta"
	Option('O', (char *)"output", NEEDS_ARG, (char *)"=output files prefix (mandatory)"),// default value "outFragments.txt"
	Option('o', (char *)"outputfile", NEEDS_ARG, (char *)"=write final fragment list and primer position list in a txt output file"), // default do not write 
	Option('P', (char *)"primers", NEEDS_ARG, (char *)"=primer input file (default primerList.fasta)"),
	Option('L', (char *)"frgLngth", NEEDS_ARG, (char *)"=fragment average length (default 70,000 nt)"), 
	Option('C', (char *)"coverage", NEEDS_ARG, (char *)"=maximum expected coverage (default 1000)"),
	Option('p', (char *)"primerNo", NEEDS_ARG, (char *)"=average input primer number (default: input dna length / 70,000 * covergae * 10)"),
	Option('R', (char *)"readLength", NEEDS_ARG, (char *)"=minimum output read length (default 10)"),
	Option('s', (char *)"stepSize", NEEDS_ARG, (char *)"=phi29 step size (default 20)"),
	Option('a', (char *)"attachNum", NEEDS_ARG, (char *)"=phi29 attach number (default 10)"),
//	Option('m', (char *)"machines", NEEDS_ARG, (char *)"number of machines (only needed when distributing among multiple machines)"),
//	Option('s', (char *)"slice", NEEDS_ARG, (char *)"0-base slice number (only needed when distributing among multiple machines)"),
//	Option('S', (char *)"single", NO_ARG, (char *)"single strand"),
 //	Option(1, (char *)"exportkmers", NO_ARG, (char *)"export kmers with their multiplicities"),
	Option(0, NULL, 0, NULL)
};

//add save and load fragment instead of fragmentasInput

Coord dnaLength = 0;
Coord frgAveLength = 7000; //change it to 70000;
Coverage aveCoverage = 600; //change it to 1000;
Coverage avePrimerNum = 0; //change it if it is not in the input
//input dna length / 70,000 * covergae * 10
//int primerLength = 0;
Coord minimumReadLength = 10;
Coverage phi29StepSize = 30; // ? Is Coverage proper type?
Coverage phi29AttachmentNum = 20*100;
bool FragmentasInput = false;
bool Writefragmentasoutput = false;
PrimerList primerList;
DNAType originalDNA;
FragmentList fragmentList;
FreePrimerPositionSets freePrimerPositionSets;
//d PrimerAvailabilityList primerAvailabilityList;
ReadsList readsList;
Coverage singleStrandCoverage;
Coverage wholeCoverage;
Phi29List phi29List;
//vector<AttachmentMerit> attachmentMeritList;
Coverage freePrimerPositionTotal = 0;
Coverage primerCurrentNoTotal = 0;


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
//	cout << "dnafile = " << dnaFile;
	Sequence *originalDNAseq = new Sequence(dnaFile);
	if (originalDNAseq->isLoaded())
	{
		char *DNAseq =  originalDNAseq->getString();
		for (dnaLength = 0 ; DNAseq[dnaLength] !=0; dnaLength++)
			originalDNA.push_back(Base(DNAseq[dnaLength]));
	}
	else
	{
		exitMsg((char *)"Error: Input sequence cannot be loaded.", INPUT_ARG_ERROR);
	}
	delete originalDNAseq;
	if(dnaFile)
		close_file(dnaFile);
}

void loadPrimers(string inputPrimerFileName)
{
	FILE * primerFile = NULL;
	primerFile = open_file(inputPrimerFileName.c_str(), "rt");
//	cout << "inputPrimerFileName = " << inputPrimerFileName;
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
	//	cout << " pindex = " << pindex << endl;
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
	//	cout << " pindex = " << pindex << endl;
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
	cout << " primerLength = " << primerLength << ", AveCountPerPrimer = " << AveCountPerPrimer << endl;
	//cout << "length of primer list " << primerList.size() << " AveCountPerPrimer =  " << AveCountPerPrimer << endl;
	for (int i = 0; i < primerList.size(); i++)
	{
	//	cout <<	"primer No = " << i << endl;
	//	writeSeq(primerList.at(i));
	//	cout << endl;
	//	cout.flush();
		FreePrimerPosition primerPositionTemp;
	//	cout << endl;
	//	cout.flush();
		primerPositionTemp.currentNo = AveCountPerPrimer + (rand() % ( AveCountPerPrimer /10)) - AveCountPerPrimer /20;
		primerCurrentNoTotal = primerCurrentNoTotal + primerPositionTemp.currentNo;
		pair<FreePrimerPositionSets::iterator, bool> ret;
		ret = freePrimerPositionSets.insert(PrimerPositionPair(primerList.at(i), primerPositionTemp));
		if (ret.second == false)
			exitMsg((char *)" error in initializing primer position list. There are dublicate primers in the input primer list.", INTERNAL_WOW_ERROR);
	}
	//cout << endl;
	for (FragmentID fn = 0; fn < fragmentList.size(); fn++)
	{	
		for (Coord fc = primerLength - 1 ; fc < fragmentList.at(fn).dna.size(); fc++)
		{
		//	cout << "fn = " << fn << ", fc = " << fc << endl;
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
	cout << ">> Total CurrentNo =  " << primerCurrentNoTotal <<"| Total available positions = " << freePrimerPositionTotal << endl;
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
	for (FragmentID k = 0; k < fragmentList.size(); k++) //change it
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
	cout << "DNA sequence loaded of length: " << dnaLength << endl;
	for (int i = 0; i < originalDNA.size(); i++)
		cout << originalDNA.at(i); 
	cout << endl;
}

void writePrimers()
{
	cout << "Primers number length: " << primerList.size() << endl;
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
		dnaTemp.dna.push_back(fbaseTemp);	
	};
	fragmentList.push_back(dnaTemp);
	singleStrandCoverage = originalDNA.size();
	wholeCoverage = singleStrandCoverage;
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

//	cout << "counter = " << counter << endl;
	
	if (counter < randomPrimerPosition) 
		exitMsg((char *)" Error in findPrimerPositionInSet. randomPrimerPosition is larger than the total number of available free positions.", INTERNAL_WOW_ERROR);
	else
	{
		primerit --;
		counter = counter - primerit->second.set.size(); // ?????????? must check if primerit->second.currentNo > 0 
		if  (primerit == freePrimerPositionSets.end())
			exitMsg((char *)" Error in findPrimerPositionInSet. An error occured related to primerit.", INTERNAL_WOW_ERROR);
		Coord randomIndex = randomPrimerPosition - counter;
	//	cout << "counter = " << counter << ", primerit->second.set.size() = " << primerit->second.set.size() << ", randomIndex = " << randomIndex << endl;
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
			cout << "randomPrimerPosition = " << randomPrimerPosition << " :: "; 
			writeSeq(primerTemp);		
			cout << ", freePrimerPositionTotal = " << freePrimerPositionTotal << ", primerit->second.set.size() = " << primerit->second.set.size() << endl;
			freePrimerPositionSets.erase(primerTemp);
			
		}
		primerCurrentNoTotal --;
		return pair<Primer, Position> (primerTemp, fragpos);
	}
}


void attachPhi29(int phi29AttachmentNumCurrent)
{
	Coverage primerIndex = 0;
	while ((primerIndex < phi29AttachmentNumCurrent) && (freePrimerPositionTotal > 0))  
	{
		Coverage randomPrimerPosition = rand() % freePrimerPositionTotal;
		pair<Primer, Position> primerpositionpair;
	//	cout << " freePrimerPositionTotal = " << freePrimerPositionTotal << ", randomPrimerPosition = " << randomPrimerPosition << endl;
		primerpositionpair = findPrimerPositionInSet(randomPrimerPosition);
	//	writeSeq(primerpositionpair.first);
	//	cout << ", (" << primerpositionpair.second.fragmentNo1 << ", " << primerpositionpair.second.pos << ")" << endl;
	//	string str;		
	//	getline(cin, str);
		Position fragPos = primerpositionpair.second;
		Primer primerAttached = primerpositionpair.first;
		Coord primerLength = primerAttached.size();
//		cout << "fragPos = (" << fragPos.fragmentNo1 << "," << fragPos.pos << ")" << endl;
		//writeFragmentList();
		/*%length of Fragment should be larger than a primer */
		if (fragPos.pos < primerLength)
			exitMsg((char *) "Error in attachPhi29. There is a primer position in the list which will generate a fragment of size 0 or less.", INTERNAL_WOW_ERROR);
		else 
		{
			bool fragmentOccupancy = false;
			for ( Coord k =  fragPos.pos;  k >= fragPos.pos - primerLength + 1; k--)
			{
				//cout << "pos = " << k << endl;
				fragmentOccupancy = fragmentOccupancy  || (fragmentList.at(fragPos.fragmentNo1 - 1).dna.at(k).occupancy.fragmentNo1 != 0);
			};
			//cout << "fragmentOccupancy = " << fragmentOccupancy << endl;
			if (fragmentOccupancy == true) /* The primer position is in the list of free positions but it is not free in real */
				exitMsg((char *) "Error in attachPhi29. The primer position is in the list of free positions but it is not free on the fragment.", INTERNAL_WOW_ERROR);
			else /* the fragment is free at this position for this primer */
			{
				primerIndex ++;
				FragmentID fragmentIndex = FragmentID(fragmentList.size());
				Phi29 newPhi29;
				newPhi29.fragmentNo1 = fragmentIndex + 1;
				if ((frgAveLength / 10) == 0) 
					newPhi29.expectedLength = frgAveLength;
				else
					newPhi29.expectedLength = frgAveLength + (frgAveLength / 20) - (rand() % (frgAveLength / 10)); 					Coord currentCoord = fragPos.pos - primerLength + 1;
				/*length of the new fragment can not be more than remaining part of the copied fragment */
				if (newPhi29.expectedLength > currentCoord)
					newPhi29.expectedLength = currentCoord;
				newPhi29.currentPosition.fragmentNo1 = fragPos.fragmentNo1;
				newPhi29.currentPosition.pos = currentCoord;
				phi29List.push_back(newPhi29);
				Fragment newFragment;
			//	cout <<"before change : "<<endl;
			//	writeFragmentList();
				for (Coord k = 0 ; k < primerLength ; k ++)
				{
					FragmentBase baseTemp;
					baseTemp.base = primerAttached.at(k);
					baseTemp.occupancy = fragPos;
					newFragment.dna.push_back(baseTemp);
					fragmentList.at(fragPos.fragmentNo1 - 1).dna.at(fragPos.pos).occupancy.fragmentNo1 = -(fragmentIndex + 1);
					fragmentList.at(fragPos.fragmentNo1 - 1).dna.at(fragPos.pos).occupancy.pos = k;
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

void OneStepAheadPhi29() //labs
{
	for (int k = 0; k < phi29StepSize; k++)
	{
		int phi29Counter = 0;
		for (list<Phi29>::iterator it = phi29List.begin(); it != phi29List.end(); it++)
		{
			Position positionOnNewFrag, positionOnOriginalFrag;
			positionOnOriginalFrag.fragmentNo1 = it->currentPosition.fragmentNo1;	
			positionOnOriginalFrag.pos = it->currentPosition.pos; // ??
			positionOnOriginalFrag.pos--;
			if (positionOnOriginalFrag.pos < 0)
				exitMsg((char *)"Error: error in OneStepAheadPhi29. positionOnOriginalFrag.pos is negative.", INTERNAL_WOW_ERROR);
			positionOnNewFrag.fragmentNo1 = it->fragmentNo1;
			positionOnNewFrag.pos = fragmentList.at(it->fragmentNo1 - 1).dna.size();
			FragmentBase newBase, originalBase;
			originalBase = fragmentList.at(positionOnOriginalFrag.fragmentNo1 - 1).dna.at(positionOnOriginalFrag.pos);
			newBase.base = reverseComplement(originalBase.base);
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
				/*%detach fragment occupancyTemp.fragmentNo1 from fragment positionOnOriginalFrag.fragmentNo1 */
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
//		cout <<" phi29 before removing " << endl;
//		writephi29List();
//		cout.flush();
	
		phi29List.remove_if(zeroExpectedLength);	
//		cout <<" phi29 after removing " << endl;
//		writephi29List();
//		cout.flush();
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
		//oldPhiListLength
		Coverage phi29AttachmentNumCurrent = phi29AttachmentNum * singleStrandCoverage / dnaLength * primerCurrentNoTotal / avePrimerNum;
		if (phi29AttachmentNumCurrent < 1)
			phi29AttachmentNumCurrent = 1;
		//cout << counter << " > phi29AttachmentNumCurrent = " << phi29AttachmentNumCurrent << endl;
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
			cout << "counter = " << counter << " fragments counts = " << fragmentList.size() << ", wholeCoverage = " << wholeCoverageRate << ", singleStrandCoverage = " << singleCoverageRate << ", phi29AttachmentNumCurrent = " << phi29AttachmentNumCurrent << ", active current phi29 numbers = " << phi29List.size() << ", primerCurrentNoTotal = " << primerCurrentNoTotal << endl;
		}
	}
	cout << "saving fragments ends" << endl;
	if (Writefragmentasoutput)
		saveFragmentList(outputFragmentsFile);
	double wholeCoverageRate = (double) wholeCoverage / (double) dnaLength;
	double singleCoverageRate = (double) singleStrandCoverage / (double) dnaLength;
	wholeCoverageRate = (double) wholeCoverage / (double) dnaLength;
	singleCoverageRate = (double) singleStrandCoverage / (double) dnaLength;
	cout << "counter = " << counter << " fragments counts = " << fragmentList.size() << ", wholeCoverage = " << wholeCoverageRate << ", singleStrandCoverage = " << singleCoverageRate << ", active current phi29 numbers = " << phi29List.size() << ", primerCurrentNoTotal = " << primerCurrentNoTotal << endl;

}


void cleaveFragments(string filename, double averageReadLength, FragmentID readNumber)
{
	FILE * readFile = NULL;
	readFile = open_file(filename.c_str(), "wt");
	(averageReadLength) = 0;
	(readNumber) = 0;
	Coverage singleStrandCounterAtTheEnd = 0;
	for (FragmentID fragIndex = 0; fragIndex < fragmentList.size(); fragIndex++)
	{	
//		cout << "fragIndex = " << fragIndex << endl;
//		cout.flush();
		bool occupied = (fragmentList[fragIndex].dna[0].occupancy.fragmentNo1 !=0);
		for (Coord posIndex = 0; posIndex < fragmentList[fragIndex].dna.size();)
		{
//			cout << "posIndex = " << posIndex << endl;
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
//				cout <<">" << readTemp.size() <<" | name = R" << readNumber << " | fragment = "<< fragIndex <<" | position = "<< posIndex << endl;
//				writeSeq(readTemp);
//				cout << endl;
//				cout.flush();
				readsList.push_back(readTemp);
				readNumber ++;
				averageReadLength = averageReadLength + readTemp.size();
				fprintf(readFile, ">%ld | name = R%ld | fragment = %ld | position = %ld \n", readTemp.size(), readNumber, fragIndex, posIndex);
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
//				cout <<">" << readTemp2.size() <<" | name = R" << readNumber << " | fragment = "<< fragIndex <<" | position = "<< posIndex << endl;
//				writeSeq(readTemp2);
//				cout << endl;
//				cout.flush();	
				readsList.push_back(readTemp2);
//				cout << "averageReadLength = " << averageReadLength << endl;
				cout.flush();	
				readNumber ++;
				averageReadLength = averageReadLength + readTemp2.size();
				fprintf(readFile, ">%ld | name = R%ld | fragment = %ld | position = %ld \n", readTemp2.size(), readNumber, fragIndex, posIndex);
				fprintfSeq(readFile,"%c",readTemp2);
				fprintf(readFile, "\n");
			}
		}
	}
	(averageReadLength) = (averageReadLength) / double (readNumber);
	cout << "averageReadLength = " << averageReadLength << ", singleStrandCounterAtTheEnd = " << singleStrandCounterAtTheEnd << endl;
	if(readFile)
		close_file(readFile);
}

int main(int argc, char *argv[])
{
	GetOpt opts(argc, argv, OPTIONS);
	
	int files = 0;
	string outputName ("out");
	string inputFileName ("dna.fasta");
	string inputPrimerFileName ("primerList.fasta");
	string outputFragmentsFile ("outFragments.txt");
	string outputReadsName ("outReads.fasta");
	string outputPrimerPositionsFile ("outPrimerPositions.txt");

	while (opts.hasNext())
	{
		Option *current = opts.next();
		char count = current->getShortForm();
		cout << count;
//		if (count == FREE_ARG)
//		{
//			if(files == 0)
//			{
//				inputFilename = current->getArg();
//				files++;
//			}
//			else
//				cerr << "Warning: ignoring additional argument " <<  current->getArg() << endl;
//		}
      		if (count == 'V')
			{} //version("V1.0"); //FILE_STRING); // ?? correct this
      		else if (count == 'h')
		{
			printf("Usage: ");
			printf(FILE_STRING);
		//	printf(" [options] reads-file\n");
		//	printf("       Output files prefix is identified by -O option.\n");
			printf("%s\n", opts.help());
		//	exitMsg(NULL, NO_ERROR);
		}
		else if (count == 'f')
			FragmentasInput = true;
		else if (count == 'I')
			inputFileName = current->getArg();
		else if (count == 'O')
		{
			outputName = current->getArg(); // concatenate with "Fragments"	
			outputFragmentsFile = outputName;
			outputReadsName = outputName; // concatenate with "Reads"
			outputFragmentsFile += "Fragments.txt";
			outputReadsName += "Reads.fasta";
			outputPrimerPositionsFile = outputName;
			outputPrimerPositionsFile += "PrimerPositions.txt";
		}
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
			phi29AttachmentNum = strtoul(current->getArg(), NULL, 0);
	}
//	if(outputName == "")
//	{
//		fputs("Error: output name not specified, try '", stderr);
//		fputs(FILE_STRING, stderr);
//		fputs(" -h' for help.\n", stderr);
//		exitMsg(NULL, OUTPUT_ARG_ERROR);
//	}
	

	//string logFilename = outputName+".log"; //add
//	cout << "FragmentasInput = " << FragmentasInput << "\n inputFileName = " << inputFileName << endl;
	cout << "outputFragmentsFile = " << outputFragmentsFile << ", outputReadsName = " << outputReadsName << ", outputPrimerPositionsFile = " << outputPrimerPositionsFile << ", inputPrimerFileName = " << inputPrimerFileName << endl;
	cout << " fragment Average Length = " << frgAveLength << ", average coverage = " << aveCoverage << ",   avePrimerNum = " << avePrimerNum << ", minimumReadLength = " << minimumReadLength << ", phi29StepSize = " << phi29StepSize << ",  phi29AttachmentNum = " << phi29AttachmentNum << endl;	
	cout << "start to Run" << endl;
	cout.flush();
	//time_t now;
	//time(&now);
	//string str;
	//getline(cin,str);
	if (!FragmentasInput) // DNA sequence as input 
	{
		cout << "loading original dna sequence" << endl;
		if(inputFileName != "")
			loadOriginalSequence(inputFileName);
		//writeDNA();
		cout << "inizializing fragment list: " << endl;
		initializeFragmentList();
		if (Writefragmentasoutput)
			saveFragmentList(outputFragmentsFile);
	
	}
	else
	{
		cout << "loading Fragments" << endl;
		if(inputFileName != "")
			loadFragmentList(inputFileName); 
		cout << "Fragments loaded" << endl;
	};

	if (avePrimerNum == 0) // put it after reading the input dna
		avePrimerNum = dnaLength * aveCoverage / frgAveLength * CoverageToPrimerConst;
	cout << "avePrimerNum = " << avePrimerNum << endl;	

	cout << "reading primers" << endl;
	if (inputPrimerFileName !="")
		loadPrimers(inputPrimerFileName);
	cout << "Primer list loaded" << endl;
//	writePrimers();
	initializePrimerPosition();
	if (Writefragmentasoutput)
		savePrimerPositionSets(outputPrimerPositionsFile);
	cout << "PrimerPositionList has been initialized" << endl;
	cout << " Amplification has been started " << endl;
	cout << " RAND_MAX = " << RAND_MAX << endl;
	srand(time(NULL));
	sequenceFragments(outputFragmentsFile, outputPrimerPositionsFile);
	cout << " Amplification has been finished " << endl;
	if (Writefragmentasoutput)
		savePrimerPositionSets(outputPrimerPositionsFile);
	double averageReadLength;
	FragmentID readNumber;
	cout << "cleavage starts." << endl; 
	cleaveFragments(outputReadsName, averageReadLength, readNumber);
	cout << "cleavage ends." << endl; 
	return 0;
}

