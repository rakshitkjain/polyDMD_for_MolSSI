//Contains all the data rading files (PDB, PSF, Mass and radii, parameter file on 28th October 2019)
#ifndef _READFILES_H
#define _READFILES_H
//Headers needed in writefiles.cpp
#include "particle.h"
#include "bond.h"
#include "nonbond.h"

class ReadFiles
{
	public:
		vector<std::string> atomtype_vector;

	void ReadPDBFile(std::string, vector<Particle>&, vector<double>&);		
	void ReadPSFFile(std::string, vector<Particle>&, vector<BOND>&, vector<NONBOND>&);		
	void ReadParameterFile(std::string, vector<Particle>&, vector<BOND>&, vector<NONBOND>&);
};
#endif

