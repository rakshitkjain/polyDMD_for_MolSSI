//Here we define the System class function declarations, function definitions are in system.cpp

//This construct ensures you do not re-define headers
#ifndef _SYSTEM_H
#define _SYSTEM_H

//Headers needed in system.h
#include "header.h"
#include "xyz.h"
#include "particle.h"
#include "grid.h"
#include "bond.h"
#include "nonbond.h"
#include "readfiles.h"

class System
{
	public:
//N=no. of particles, nsweep = no. of sweeps (number of collisions = nsweep*N), maxbin = no. of bins for g(r) calculation
		int N, nsweep, maxbin;	 		
//L= Box size
		double L;				
//temperature = absolute temperature of the system, mass = mass of each particle in amu
		double temperature/*, mass*/;		
//Potential Parameters "Square well" well_depth= energy parameter (depth = +ve, WILL NOT WORK FOR SHOULDER POTENTIALS!), sigma1 = particle diameter, sigma2 = square well diameter
//		double well_depth, sigma1, sigma2;	
//bondlength of the polymer chain, bond_delta = difference in bond length allowed
		double /*bondlength, */bond_delta; 	
//chainlength = Number of particles in chain, nchains = number of chains in simulation box
//		int chainlength, nchains; 		
//One out of every andersen_freq collisions will be a ghost collision
		int andersen_freq;			
//Size of each cell for cell-list
		double cellwidth;			
//Multiplicative factor for the largest potential width giving the range for the neighbor list calculation
		double neighborlist_factor;		
//Maximum velocity allowed, to calculate nbrlist update frequency
		double maxvel;				
//Filenames for pdb, psf and parameter file
		std::string filename_pdb, filename_psf;	
//Partial filename for output files
		std::string filename_out, filename_par;	
//Printing frequency, false positioning update frequency
		int fpupdate_freq;
		double print_freq;				
//Fraction of nsweep to be used as burning collisions
		double burn_percent;			
//Particle positions
		vector<Particle> P;			
//Cell parameters like a, b, c, alpha, beta, gamma
		vector<double> cellparameters;		
//Grid variables
		vector<Grid> G;				
//List of bonded particles (NEED TO INCLUDE PSEUDOBONDS AS WELL)
		vector<BOND> bondlist;			
		vector<NONBOND> nonbondlist;
		vector<int> LinkList;
		vector<int> HeadList;
		vector<int> NeighborLinkList;		
//Contains the position of the last neighbor of particle i in the vector NeighborLinkList
		vector<int> NeighborHeadList;		
//square of square well diameter, needed for time calculations
//		double sigma2sq;			
//Potential energy of the system
		double potential_energy = 0;		
//This variable determines if we are using celllist or not, depending on the initial state of the system
		bool celllist_counter = false;		
		bool neighborlist_counter = false;	
//Distance of neighborlist = neighborlist_factor*largestpotentialwidth
		double neighborlist_range;		
//Maximum number of cells allowed in the system
		int ncell_max;				
//Length of the sides of the box
//		double a, b, c;				
//Angles for the box
//		double alpha, beta, gamma;		
//time passed in the system, time of last false position update
		double TIME, fpupdate_TIME;				
//No. of chains per side of the lattice
//		int Nchain_perside;			
//Maximum of all the well depths in the system (Will be important when multiple potentials etc)
		double max_well_width;			

//Additional functions
//Single coordinate pbc
	void PBC_s(double&);							
//3-D PBC
	void PBC(double&, double&, double&);					
//Minimum image convention
	void min_img(XYZ&, double);		
//Generate a gaussian number for a gaussian velocity profile because the inbuilt normal dist. not working well
	double RandomGaussianNumber();				
//send the square of the distance beyond which it counts as an overlap, Checking particle overlaps
	bool CheckOverlap(vector<Particle>&); 					
//Checking initialising chain overlap
	bool CheckOverlap_bead(int, vector<Particle>&);				
//Potential energy calculating function
	double pe(vector<Particle>&, int);					
//Kinetic Energy calculating function
	double ke(vector<Particle>&, int);					
//Net simulation box momentum
	VEL net_momentum(vector<Particle>&, int);				
//Creating the initial system within the simulation
	void Create();								
//Initialising chain velocities			
	void velocity_initialization();						
//Initialising chains at the start
	void coordinate_initialization(double, int); 				
//Cell list functions								
	//Iniitialisation for headlist and linklist
	void CellList_initialization();						
	//Makes cellID, header and link-list arrays
	void CellList();							
//Neighbor list functions
	void NeighborList_initialization();
	void NeighborList();
//False position updater (forward and backward)
	XYZ OneParticlePositionUpdater(Particle&, double, double);
	XYZ OneParticlePositionBackwarder(Particle&, double, double);
};
#endif
