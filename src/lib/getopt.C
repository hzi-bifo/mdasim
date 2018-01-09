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
 * Title:          getopt.C 
 * Author:         Hamidreza Chitsaz
 * Created:        2011
 * Last modified:  10/10/2011
 *
 * Copyright (c) 2011- Hamidreza Chitsaz
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include <stdio.h>
#include "getopt.h"


void GetOpt::parseArgs()
{
	optionsNum = -1;
	while(options[++optionsNum].getShortForm());

	int i = 0;
	char buffer[MAX_BUF_SIZE];
	
		
	while (i++ < argc - 1)
	{
		if(argv[i][0] == '-')
		{
			if(argv[i][1] == '-')
			{
				int j = 2;
				bool hasArg = false;
				char a[MAX_BUF_SIZE];
				
				while(argv[i][j] && (j < MAX_BUF_SIZE))
				{
					if(argv[i][j] == '=')
					{
						hasArg = true;
						break;
					}
					else
					{
						buffer[j-2] = argv[i][j];
						j++;
					}
				}
				buffer[j-2] = 0;
				
				j++;
				int k = 0;

				while(argv[i][j] && (j < MAX_BUF_SIZE))
				{
					a[k++] = argv[i][j++];
				}
				a[k] = 0;

				for(j = 0; j < optionsNum; j++)
					if(!strcmp(buffer, options[j].getLongForm()))
						break;

				if(j == optionsNum)
				{
					fprintf(stderr, "Undefined option --%s\n", buffer);
					exit(1);
				}
				if (options[j].getArgRequirement() == NO_ARG && hasArg)
				{
					fprintf(stderr, "Option --%s does not accept arguments.\n", buffer);
					exit(1);
				}

				if (options[j].getArgRequirement() == NEEDS_ARG && !hasArg)
				{
					fprintf(stderr, "Option --%s requires argument.\n", buffer);
					exit(1);
				}
				
				parsedOptions[parsedOptionsNum] = new Option(options[j].getShortForm(), options[j].getLongForm(), options[j].getArgRequirement(), options[j].getDesc());
				if(hasArg)
					parsedOptions[parsedOptionsNum]->setArg(a);
				parsedOptionsNum++;
			} 
			else if(argv[i][1])
			{
				int j;
				for(j = 0; j < optionsNum; j++)
					if(argv[i][1] == options[j].getShortForm())
						break;
				
				if(j == optionsNum)
				{
					fprintf(stderr, "Undefined option -%c\n", argv[i][1]);
					exit(1);
				}
				
				parsedOptions[parsedOptionsNum] = new Option(options[j].getShortForm(), options[j].getLongForm(), options[j].getArgRequirement(), options[j].getDesc());
				
				char a[MAX_BUF_SIZE];
				j = 2;
				bool start = false;
				int k = 0;

				while(argv[i][j] && (j < MAX_BUF_SIZE))
				{
					if(start)
						a[k++] = argv[i][j++];
					else
						if(argv[i][j++] == '=')
							start = true;
				}
				a[k] = 0;

				if (parsedOptions[parsedOptionsNum]->getArgRequirement() == NO_ARG && start)
				{
					fprintf(stderr, "Option -%c does not accept arguments.\n", argv[i][1]);
					exit(1);
				}

				if (parsedOptions[parsedOptionsNum]->getArgRequirement() == NEEDS_ARG && !start)
				{
					fprintf(stderr, "Option -%c requires argument.\n", argv[i][1]);
					exit(1);
				}
				
				if(start)
					parsedOptions[parsedOptionsNum]->setArg(a);

				parsedOptionsNum++;
			}
		}
		else
		{
			parsedOptions[parsedOptionsNum] = new Option(FREE_ARG, NULL, 0, NULL);
			parsedOptions[parsedOptionsNum++]->setArg(argv[i]);
		}

	}
}

char *GetOpt::help()
{
	char buf[MAX_BUF_SIZE];

	strcpy(helpstring, "\n\n");
	for(int j = 0; j < optionsNum; j++)
	{
		buf[0] = 0;
		if((options[j].getShortForm() <= 'Z' && options[j].getShortForm() >= 'A') || (options[j].getShortForm() <= 'z' && options[j].getShortForm() >= 'a'))
		{
			sprintf(buf, "\t-%c,--%s      %s\n", options[j].getShortForm(), options[j].getLongForm(), options[j].getDesc());
			strcat(helpstring, buf);
		}
	}

	for(int j = 0; j < optionsNum; j++)
	{
		buf[0] = 0;
		if(!((options[j].getShortForm() <= 'Z' && options[j].getShortForm() >= 'A') || (options[j].getShortForm() <= 'z' && options[j].getShortForm() >= 'a')))
		{
			sprintf(buf, "\t   --%s      %s\n", options[j].getLongForm(), options[j].getDesc());
			strcat(helpstring, buf);
		}
	}

	return helpstring;
}

