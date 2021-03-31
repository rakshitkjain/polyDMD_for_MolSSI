//Here we define the Coll class function declarations, function definitions are in collision.cpp
//This construct ensures you do not re-define headers
#ifndef _COLLISION_H
#define _COLLISION_H

//Headers needed in collision.h
#include "timecalc.h"

class Collision
{
	public:
		TimeCalc TC;
		//Virial accumulator
		double virial = 0;
		//Tally for bond collisions and square well collisions
		int bondcol_counter = 0, swcol_counter = 0, cellchange_counter = 0, nbrmove_counter = 0;
		//Checking if a cellchange or nbrchange move is happening or not
		bool didcellchange = false, didnbrchange = false;
		double	nbrmove_distance = 0.0;

//Recalculating "b" after particle movements
	double recalc_b(Particle&, double, int, double, double);				
//Moving the particles and momentum transfer
	void Bump(double, int);									
//Thermostat function, returning the particle on which applying andersen
	int AndersenThermostat();		
//Function to update particle position when using false position
	void AllParticlePositionUpdater(int, vector<Particle>&, double, double);
};
#endif
