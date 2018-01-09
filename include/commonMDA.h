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
 * Title:          commonMDA.h  
 * Author:         (modified by) Zeinab Taghavi
 * Created:        2012
 * Last modified:  05/27/2014
 *
 * Copyright (c) 2012- Zeinab Taghavi
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/


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
 * Title:          commonMDA.h 
 * Author:         Hamidreza Chitsaz
 * Created:        2011
 * Last modified:  10/10/2011
 *
 * Copyright (c) 2011- Hamidreza Chitsaz
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/


#ifndef COMMONMDA_H
#define COMMONMDA_H

//works only for gcc
#define PACK_MEMORY 1

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string>

#ifndef NULL
#define NULL 0
#endif

#ifndef true
#define true 1
#define false 0
#endif

#define PACKAGE_STRING "MDAsim 1.2"

//#define MIN(A, B)  (((A) < (B)) ? (A) : (B))
//#define MAX(A, B)  (((A) < (B)) ? (B) : (A))

void version(char* prog);

typedef enum {
	NO_ERROR = 0,
	FILE_OPEN_ERROR = 1, 
	UNITIGS_FILE_ERROR = 2, 
	CONFIG_FILE_ERROR = 3, 
	INTERNAL_DIJKSTRA_ERROR = 4, 
	INTERNAL_RADIXHEAP_ERROR = 5,
	INTERNAL_SPLAYTREE_ERROR = 6,
	INPUT_ARG_ERROR = 7,
	OUTPUT_ARG_ERROR = 8,
	INTERNAL_WOW_ERROR = 9
} Error;


FILE *open_file(char *fname, const char *mode);
FILE *open_file(std::string fname, const char *mode);
void exitMsg(char *, Error);


#endif

