//Contains the class to read input parameters via a command line parser
#ifndef _READINPUT_H
#define _READINPUT_H
//Headers needed in readinput.cpp
#include "system.h"
#include "header.h"
class ReadInput
{
	public:
		System initialize;
//Help message in case of -h flag
	void printHelpMessage();				
//TO set the parameters from the command line
	int ReadVariables(int, vector<string>);
};
#endif
