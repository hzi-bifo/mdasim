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
 * Title:          alloc.h 
 * Author:         Hamidreza Chitsaz
 * Created:        2011
 * Last modified:  03/14/2012
 *
 * Copyright (c) 2011- Hamidreza Chitsaz
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#ifndef ALLOC_H
#define ALLOC_H

#include <stdlib.h>
#include <stdio.h>

class Alloc {
public:
	

	inline	void* xcalloc(size_t m, size_t n)
	{
		void* ptr;

		if(m == 0 || n == 0)
			return NULL;

		if (!(ptr = calloc(m, n)))
		{
			fputs("Error in calloc()\n", stderr);
			exit(EXIT_FAILURE);
		}

		return ptr;
	}

	inline void* xmalloc(size_t n)
	{
		void* ptr;

		if (!(ptr = malloc(n)))
		{
			fputs("Error in malloc()\n", stderr);
			exit(EXIT_FAILURE);
		}
		return ptr;
	}

	inline void* xrealloc(void* ptr, size_t n)
	{
		if (!(ptr = realloc(ptr, n)) && n != 0)
		{
			printf("%ld\n", n);
			fputs("Error in realloc()\n", stderr);
			perror("");
			exit(EXIT_FAILURE);
		}
		return ptr;
	}

	inline void xfree(void* &ptr)
	{
		if(ptr) free(ptr);
			ptr = NULL;
	}

};

#endif
