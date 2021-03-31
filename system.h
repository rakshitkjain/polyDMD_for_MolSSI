//Here we define the System class function declarations, function definitions are in system.cpp (Revision Date: 31st March 2021)
#ifndef _SYSTEM_H
#define _SYSTEM_H
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
//Class variables
	int N, nsweep, maxbin;	 		//No. of particles, No. of sweeps (number of collisions = nsweep*N), No. of bins for g(r) calculation
	double L;				//Box size
	double temperature;			//Absolute temperature of the system
	double bond_delta; 			//Difference in bond length allowed
	int andersen_freq;			//One out of every andersen_freq collisions will be a ghost collision
	double cellwidth;			//Size of each cell for cell-list
	double neighborlist_factor;		//Multiplicative factor for the largest potential width giving the range for the neighbor list calculation
	double maxvel;				//Maximum velocity allowed, to calculate nbrlist update frequency
	std::string filename_pdb, filename_psf;	//Filenames for pdb and psf file
	std::string filename_out, filename_par;	//Partial filename for output files, Filename for par file
	int fpupdate_freq;			//False positioning update frequency
	double print_freq;			//Printing frequency	
	double burn_percent;			//Fraction of nsweep to be used as burning collisions
	vector<Particle> P;			//Particle information vector
	vector<double> cellparameters;		//Cell parameters like a, b, c, alpha, beta, gamma
	vector<Grid> G;				//Grid variables
	vector<BOND> bondlist;			//List of bonded particles (NEED TO INCLUDE DIHEDRAL PSEUDOBONDS AS WELL)
	vector<NONBOND> nonbondlist;		//List of nonbonded particles 
	vector<int> LinkList;			//Storing Cell List
	vector<int> HeadList;			
	vector<int> NeighborLinkList;		//Storing Neighbor List
	vector<int> NeighborHeadList;		//Contains the position of the last neighbor of particle i in the vector NeighborLinkList
	double potential_energy = 0;		//Potential energy of the system
	bool celllist_counter = false;		//This variable determines if we are using celllist or not, depending on the initial state of the system
	bool neighborlist_counter = false;	//This variable determines if we are using nbrlist or not, depending on the initial state of the system
	double neighborlist_range;		//Distance of neighborlist = neighborlist_factor*largestpotentialwidth
	int ncell_max;				//Maximum number of cells allowed in the system
	double TIME, fpupdate_TIME;		//Time passed in the system, Time of last false position update		
	double max_well_width;			//Maximum of all the well depths in the system (Will be important when multiple potentials etc)

//Functions
	void PBC_s(double&);						//Single coordinate pbc	
	void PBC(double&, double&, double&);				//3-D PBC	
	void min_img(XYZ&, double);					//Minimum image convention
	double RandomGaussianNumber();					//Generate a gaussian number for a gaussian velocity profile
	bool CheckOverlap(vector<Particle>&); 				//Checking particle overlaps	
	bool CheckOverlap_bead(int, vector<Particle>&);			//Checking initialising chain overlap	
	double pe(vector<Particle>&, int);				//Potential energy calculating function	
	double ke(vector<Particle>&, int);				//Kinetic Energy calculating function	
	VEL net_momentum(vector<Particle>&, int);			//Net simulation box momentum	
	void Create();							//Creating the initial system within the simulation	
	void velocity_initialization();					//Initialising chain velocities				
	void CellList_initialization();					//Initialisation for celllist headlist and linklist	
	void CellList();						//Makes cellID, header and link-list arrays	
	void NeighborList_initialization();				//Neighbor list functions
	void NeighborList();
	XYZ OneParticlePositionUpdater(Particle&, double, double);	//False position updater (forward and backward)
	XYZ OneParticlePositionBackwarder(Particle&, double, double);
};
#endif
