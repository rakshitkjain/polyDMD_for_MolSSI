//MD: particle.h Particle Class (Revision Date: 4th Nov 2019)
//Store the list of particles 
#ifndef _PARTICLE_H
#define _PARTICLE_H
#include "header.h"
#include "xyz.h"
#include "vel.h"
class Particle
{
public:
	XYZ coordinate;				//coordinates
	VEL velocity;					//velocity
	double velocity2;				//net speed of the particle
	int cell_ID;					//ID of cell list location
	std::string type;				//Type of atom: ATOM/HETATOM etc.
	std::string name;				//Atom name which incorporates the geometry (N2, C3 etc)
	std::string resname;				//Residue name for the particle (CT2, ACY etc)
	char chaintype;				//Chaintype, to identify multiples of same residue (C, A, W etc)
	int chainnumber;				//Each chain has a different number and different residues have a different number
	double occupancy;				//Still NOT CLEAR ON THESE TWO
	double tempfactor;				//Still NOT CLEAR ON THESE TWO
	std::string moltype;				//Name of molecule, helps to join residues which are joined
	std::string atomtype;				//Special name for the atom, useful when reading parameters fron .par file. Initialized when reading a psf
	std::string symbol;				//Symbol of the element concerning the molecule
	double mass;					//Mass of the particle
	double diameter;				//Diamter of the particle
	double charge;					//Charge on the particle (might be useful for Ewald Summation later)
	double bcalc(){return (coordinate.x*velocity.vx + coordinate.y*velocity.vy + coordinate.z*velocity.vz) ;}
	Particle(){coordinate.x=0.0; coordinate.y=0.0; coordinate.z=0.0; velocity.vx=0.0; velocity.vy=0.0; velocity.vz=0.0; velocity2=0.0; cell_ID=0; type="ATOM"; name="C1"; resname="CHI"; chaintype='A'; chainnumber=1; occupancy=0.0; tempfactor=0.0; moltype="CHI"; atomtype = "CG331"; symbol="C"; diameter = 1.0; mass = 1.0; charge = 0.0;}
};
#endif

