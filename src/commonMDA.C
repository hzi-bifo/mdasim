/*
Copyright 2018- Victoria Sack (victoria.sack@helmholtz-hzi.de)

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
 * Title:          commonMDA.C
 * Author:         (modified by) Victoria Sack
 * Created:        2018
 * Last modified:  01/10/2018
 *
 * Copyright (c) 2018- Victoria Sack
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
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
 * Title:          commonMDA.C
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
 * Title:          commonMDA.C
 * Author:         Hamidreza Chitsaz
 * Created:        2011
 * Last modified:  10/10/2011
 *
 * Copyright (c) 2011- Hamidreza Chitsaz
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <cctype>
#include <errno.h>

#include "commonMDA.h"

void version(char* prog)
{
        //******************************************************
        // altered to #2.0

        printf("%s (%s) \n", prog, PACKAGE_STRING);

        puts("Author: D.Lähnemann and V.Sack");
        puts("Copyright (C) 2018");
        puts("Helmholtz Centre for Infection Reseach");
        puts("Braunschweig, Germany");
        puts("based on:");

        //******************************************************

        printf("mdasim 1.2 (MDAsim 1.2) \n");
	puts("Author: Zeinab Taghavi");
	puts("Copyright (C) 2012-2013");
	puts("Wayne State University");
	puts("Detroit, MI");
}

FILE *open_file(char *fname, const char *mode)
{
	FILE *ret = fopen(fname, mode);
	if(!ret)
	{
		char buf[10000] = "Error: failed to open file ";
		strncat(buf, fname, 9500);
		strcat(buf, " in mode ");
		strncat(buf, (char *)mode, 100);
		strcat(buf, ".\n");
		exitMsg(buf, FILE_OPEN_ERROR);
	}
	return ret;
}

FILE *open_file(std::string fname, const char *mode)
{
	return open_file((char*)fname.c_str(), mode);
}

void exitMsg(char *msg, Error c)
{
	if(msg)
	{
		fputs(msg, stderr);
		fputc('\n', stderr);
	}
	exit(c);
}
