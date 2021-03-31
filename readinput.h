//MD: Input Parameters reading Class. Contains the class to read input parameters via a command line parser (Revision Date: 31st March 2021)
#ifndef _READINPUT_H
#define _READINPUT_H
#include "system.h"
#include "header.h"

class ReadInput
{
	public:
//Import from other classes
		System initialize;
//Functions
	void printHelpMessage();			//Help message in case of -h flag	
	int ReadVariables(int, vector<string>);		//To set the parameters from the command line
};
#endif
