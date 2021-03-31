//MD: Radial Distribution Function calculation module. Calculates RDF from particle positions during the run (Revision Date: 31st March 2021)
#include "system.h" 
#include "timecalc.h" 
#include "particle.h"

class RDF
{
	public:
//Import from other classes
		TimeCalc mi;
//Class Variables
	vector<int>hist;		//Used to tally the radial distribution function
	vector<double>gr;		//Exact radial distribution function
	vector<double>r;		//Distance between two particles
	int rdf_steps = 0;
	double delta_r;
//Functions
	void RDF_ini(double, int);									//Initialising the required vectors and lists
	void RDF_r(double, int, const vector<Particle> , int);						//Calculating RDF
	void RDF_write(double, int, int, std::string);							//Writing in the file
};
