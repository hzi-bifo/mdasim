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


/**
 * Open file, load sequence name and sequence string and close file.
 *
 * @param fn   filename as string of chars
 */
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

/**
 * Load sequence name and sequence string.
 *
 * Wrapper to load sequence name and sequence string from a FILE handle. Opening
 * and closing of the file need to be handled externally.
 *
 * @param f   FILE handle of a file open for reading.
 */
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

/**
 * Load one sequence and its name from a fasta input file.
 *
 * Main function to load one sequence and its name from a FILE handle
 * of a FASTA input file.
 *
 * @param file   open FILE handle
 * @return   `-1` for error, `1` for success
 */
int Sequence::load(FILE* file)
{
  int current, last;
  unsigned int availableSeq, availableName;

  string = (char *)alloc.xmalloc(1024);
  availableSeq = 1024;
  // counter of nucleotides from the current sequence we already read
  int nucleotides_read = 0;
  name = (char *)alloc.xmalloc(80);
  availableName = 80;
  // counter of chars from the current sequence name we already read
  unsigned int name_chars_read = 0;
  // indicator whether we are currently reading a header (1) or not (0)
  unsigned int reading_header = 0;
  last = '\n';

  // read in currecnt char while it is not EOF
  while ((current = fgetc(file)) != EOF)
  {
    // if last char read is \n and current is '>' (*): we have read file and encounter a read name header
    if (last == '\n' && current == '>')
    {
      if (nucleotides_read)
      {
        //ungets the read name starting char '>'
        ungetc('>', file);
        //break here, because a '>' in the middle of the file has been read that does not start a sequence name!!
        break;
      }
      else
      {
        //if no nucleotides_read, we are now reading_header
        reading_header = 1;
      }
    }
    // encounter a newline
    else if (current == '\n')
    {
      // set the reading_header to 0, because we don't know what we'll get next
      reading_header = 0;
    }
    // if we are not reading_header
    else if (reading_header == 0)
    {
      // if current char is from the set [a-zA-Z0-9&]
      if (('a' <= current && current <= 'z') || ('A' <= current && current <= 'Z') || ('0' <= current && current <= '9') || current == '&')
      {
        // add capitalised current char to end of string and post-increment nucleotides_read
        string[nucleotides_read++] = toupper(current);
      }
      // TODO: this is a strange thing to test for in a FASTA file...
      else if (current == ';')
      {
        // break, something went wrong!
        break;
      }
    }
    // we previously determined that we are now reading a header
    else if (reading_header == 1)
    {
      // add the current chard to the end of the name string and post-increment name_chars_read
      name[name_chars_read++] = current;
    }

    //reallocate memory for name if the array is full
    checkArray(&name, &availableName, name_chars_read, 80);
    //reallocate memory for sequence if the array is full
    checkArray(&string, &availableSeq, nucleotides_read, 1024);
    // remember the current char as `last` before reading in a new one
    last = current;
    // repeat till end of file is reached
  }
  // set last character of name string to 0
  name[name_chars_read] = 0;
  // same with the sequence
  string[nucleotides_read] = 0;

  // if the name length is 0
  if (strlen(name) == 0)
  {
    // free the memory
    free(name);
    // set name to NULL
    name = NULL;
  }
  else
  {
    //allocate exactly as much memory as you need for the name
    name = (char *)alloc.xrealloc(name, strlen(name) + 1);
  }

  //do the same stuff for the sequence!
  if (strlen(string) == 0)
  {
    free(string);
    string = NULL;
    //however, if there is no sequence, the file is faulty
    return -1;
  }
  else
  {
    string = (char *)alloc.xrealloc(string, strlen(string) + 1);
  }
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
