/*
Copyright 2011- Hamidreza Chitsaz (chitsaz@wayne.edu)

    This file is part of HyDA.

    HyDA is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    HyDA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HyDA; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/


/***************************************************************************
 * Title:          sequence.C 
 * Author:         Hamidreza Chitsaz
 * Created:        2011
 * Last modified:  10/10/2011
 *
 * Copyright (c) 2011- Hamidreza Chitsaz
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include <string.h>
#include <ctype.h>
#include <iostream>
#include "sequence.h"


Sequence::Sequence(char *fn)
{
	std::cout << "here" << std::endl;
	name = NULL;
	string = NULL;
	FILE *f = fopen(fn, "rt");
	if(!f)
		perror("Sequence: error opening sequence file: ");

	loaded = (load(f) != -1);
	if(loaded) 
		fclose(f);
}

Sequence::Sequence(FILE *f)
{
	name = NULL;
	string = NULL;
	loaded = (load(f) != -1);
}


Sequence::~Sequence()
{
	if(loaded)
	{
		free(name);
	}
	free(string);
}

int Sequence::load(FILE *f)
{
	return input(f);
}

char* Sequence::getName()
{
	return name;
}

void Sequence::setName(const char *n)
{
	name = (char *)alloc.xrealloc(name, strlen(n)+1);
	strcpy(name, n);
}

char* Sequence::getString()
{
	return string;
}

unsigned int Sequence::getLen()
{
	return strlen(string);
}

int Sequence::input(FILE* file)
{
  int current, last, state;
  unsigned int availableSeq, usedSeq, availableName, usedName;

  string = (char *)alloc.xmalloc(1024);
  availableSeq = 1024;
  usedSeq = 0;
  name = (char *)alloc.xmalloc(80);
  availableName = 80;
  usedName = 0;
  state = 0;
  last = '\n';

  while ((current = fgetc(file)) != EOF)        //while the current char is not EOF
    {
      if (last == '\n' && current == '>')       //if last char read is \n and current is > (*)
	{
          if (usedSeq)                          //if usedSeq != 0
	    {
	      ungetc('>', file);                //ungets the read char >
              break;                            //break here, because a '>' in the middle of the file has been read that does not start a sequence name!!
	    }
	  else
            state = 1;                          //if not usedSeq, set state to 1, which means,  we're reading a sequence NAME atm
	}
      else if (current == '\n')                 //if the current char is a newline \n
        state = 0;                              //set the state to 0, because we don't know what we'll get next
      else if (state == 0)                      //if the state is IN_SEQUENCE and we did not read the start of a sequence name in (*)
	{
          if (('a' <= current && current <= 'z') || ('A' <= current && current <= 'Z') || ('0' <= current && current <= '9') || current == '&') //if current char is out of a-zA-Z0-9&
            string[usedSeq++] = toupper(current);       //increase usedSeq++, the string at usedSeq is the CAPITAL char of current char
          else if (current == ';')                      //else if current is a ;
            break;                                      //break, something went wrong!
	}
      else if (state == 1)                      //if the state is 1, which means we're IN_HEADER
        name[usedName++] = current;             //the name at usedName++ is the current char

      checkArray(&name, &availableName, usedName, 80);      //reallocate memory for name if the array is full
      checkArray(&string, &availableSeq, usedSeq, 1024);    //reallocate memory for sequence if the array is full
      last = current;                           //set the last read char last to the current one
    }                                           //repeat till end of file is reached
  name[usedName] = 0;                 //set name after usedName to 0
  string[usedSeq] = 0;                //same with the sequence
  
  if (strlen(name) == 0)            //if the name length is 0
    {
      free(name);                   //free the memory
      name = NULL;                  //set name to NULL
    }
  else
    name = (char *)alloc.xrealloc(name, strlen(name) + 1);  //else, allocate exactly as much memory as you need for the name

  if (strlen(string) == 0)          //do the same stuff for the sequence!
    {
      free(string);
      string = NULL;
      return -1;                    //however, if there is no sequence, the file is faulty
    }
  else
    string = (char *)alloc.xrealloc(string, strlen(string) + 1);
  return 1;
}

void Sequence::checkArray(char** array, unsigned int* available, unsigned int used, unsigned int increment)
{
  if (used == *available)
    {
      *array = (char *) alloc.xrealloc(*array, *available + increment);
      *available += increment;
    }
}

int Sequence::countNs()
{
	int ret = 0;
	unsigned int i = 0;

	while(string[i])
	{
		char c = toupper(string[i]);
		if(c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'U') ret++;
		i++;
	}

	return ret;
}


void Sequence::print(FILE *out)
{
	if(!out)
		return;

	fprintf(out, ">%s\n", name);
	unsigned int i = 0;	
	while(string[i])
	{
		fputc(string[i], out);
		if(!((i+1) % 80))
			fputc('\n', out);
		i++;
	}
	fputc('\n', out);
}

void Sequence::print(std::ostream &out)
{
	out << ">" << name << std::endl;
	printSequence(out);
}

void Sequence::printSequence(std::ostream &out)
{
	unsigned int i = 0;
	while(string[i])
	{
		out << string[i];
		if(!((i+1) % 80))
			out << std::endl;
		i++;
	}
	out << std::endl;
}


void Sequence::subsequence(size_t pos, size_t n)
{
	int len = getLen();

	if(pos > len)
		pos = len;
	if(pos+n > len)
		n = len-pos;

	char *newstr = (char *) alloc.xmalloc((n+1)*sizeof(char));
	for(size_t i = 0; i < n; i++)
		newstr[i] = string[i+pos];
	newstr[n] = '\0';
	if(string) alloc.xfree((void *&)string);
	string = newstr;
}



