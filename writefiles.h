//Contains all the data writing files (DCD, PDB, DUMP, Restart, TimeVarFile on 8th October 2019)
#ifndef _WRITEFILES_H
#define _WRITEFILES_H
//Headers needed in writefiles.cpp
#include "system.h"

class WriteFiles
{
	public:
		System S;				//To get all the variables for naming files
		int dcd_sets = 0;			//Number of dcd sets written

//For writing headers of DCD Files		
	void DCDHeader(double, int, std::string);					
//For writing each relevant timestep
	void DCDWriteStep(vector<Particle>&, int, std::string);				
//For combining header and temp dcd file and remove them
	void DCDFile(std::string);								
//For PDB initial system, to go along with dcd
	void PDBFile(vector<Particle>&, std::string);						
//Initialising the file to store time variation parameters
	void TimeVarFileIni(std::string);							
//Storing in the above file
	void TimeVarFile(double, double, double, double, VEL, std::string);			
//Writing restart file
	void RestartFile(double, vector<Particle>&, int, std::string);				
//For movies (LAMMPS Traj format)
	void WriteDump(double, vector<Particle>&, int, double, std::string);			
};
#endif
	
