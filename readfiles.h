//MD: File Reading Class. Contains all the data rading files (PDB, PSF, Mass and radii, parameter file) (Revision Date: 31st March 2021)
#ifndef _READFILES_H
#define _READFILES_H
#include "particle.h"
#include "bond.h"
#include "nonbond.h"

class ReadFiles
{
	public:
//Class variables
	vector<std::string> atomtype_vector;
//Functions
	void ReadPDBFile(std::string, vector<Particle>&, vector<double>&);		
	void ReadPSFFile(std::string, vector<Particle>&, vector<BOND>&, vector<NONBOND>&);		
	void ReadParameterFile(std::string, vector<Particle>&, vector<BOND>&, vector<NONBOND>&);
};
#endif

