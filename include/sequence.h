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
 * Title:          sequence.h 
 * Author:         Hamidreza Chitsaz
 * Created:        2011
 * Last modified:  12/27/2011
 *
 * Copyright (c) 2011- Hamidreza Chitsaz
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#ifndef SEQ_H
#define SEQ_H

#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "alloc.h"

class Sequence {
private:
	Alloc alloc; 
	char *name;
	char *string;
	bool loaded;
	int readSequence(FILE *);
	void checkArray(char** array, unsigned int* available, unsigned int used, unsigned int increment);

public:
	Sequence(char *fn);
	Sequence(FILE*);
	~Sequence();
	int load(FILE *f);
	int input(FILE* file);
	bool isLoaded() {return loaded;};
	char *getName();
	void setName(const char *);
	char *getString();
	unsigned int getLen();
	void print(FILE *);
	void print(std::ostream &out);
	void subsequence(size_t pos, size_t n);
	void printSequence(std::ostream &out);
	int countNs();
};

#endif
