//MD: File writing class. Contains all the data writing files (DCD, PDB, DUMP, Restart, TimeVarFile) (Revision Date: 31st March 2021)
#ifndef _WRITEFILES_H
#define _WRITEFILES_H
#include "system.h"

class WriteFiles
{
	public:
//Import from other classes
		System S;			//To get all the variables for naming files
//Class variables	
	int dcd_sets = 0;			//Number of dcd sets written
//Functions
	void DCDHeader(double, int, std::string);				//For writing headers of DCD Files					
	void DCDWriteStep(vector<Particle>&, int, std::string);			//For writing each relevant timestep	
	void DCDFile(std::string);						//For combining header and temp dcd file and remove them		
	void PDBFile(vector<Particle>&, std::string);				//For PDB initial system, to go along with dcd
	void TimeVarFileIni(std::string);					//Initialising the file to store time variation parameters		
	void TimeVarFile(double, double, double, double, VEL, std::string);	//Storing in the above file		
	void RestartFile(double, vector<Particle>&, int, std::string);		//Writing restart file		
	void WriteDump(double, vector<Particle>&, int, double, std::string);	//For movies (LAMMPS Traj .dump format)		
};
#endif
	
