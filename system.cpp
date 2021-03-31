//Here we define function definitions of system class
#include "system.h"
//Single direction pbc, Complicated functions because we are moving large distances for 
//low densities using "&" to actually edit the values, not just calculate the new positions
void System::PBC_s(double &a)
{
	double adec;
	int aint;
	if(a>L)
	{
		aint = floor(a);
		adec = a-double(aint);
		a = double(aint%int(L)) + adec;
	}   
	else
	{
		if(a< 0.0)
		{
			aint = floor(a); 
			adec = -a+double(aint);
			a = L + double(aint%int(L)) - adec;
		}
	}
}
//3-D PBC, uses the 1-D PBC function
void System::PBC(double &q, double &w, double &e)
{
	PBC_s(q);
	PBC_s(w);
	PBC_s(e);
}
//Minimumum image convention
void System::min_img(XYZ &r, double L)
{
	r.x = r.x - L*round(r.x/L);					
	r.y = r.y - L*round(r.y/L);
	r.z = r.z - L*round(r.z/L);
}

double System::RandomGaussianNumber()
{
	std::random_device rd{};					
   	std::mt19937 gen{rd()};
//	std::mt19937 gen{0};

	std::uniform_real_distribution<> d{-1,1};
	double ran1,ran2,ransq;

	do
	{
		ran1=d(gen);
		ran2=d(gen);
		ransq=ran1*ran1+ran2*ran2;
	}
	while((ransq>1.0)||(ransq<0.0));
	return ran1*sqrt(-2.0*log(ransq)/ransq);
}
//To make sure that the particles do not overlap at any time
bool System::CheckOverlap(vector<Particle> &P)
{
	double tolerance = 1.0e-4;
	int i, j, k;
	bool overlap = false, particle_bonded = false, particle_nonbonded = false;
	XYZ counter;
	double rij, rij1, l1, l2;

	for(i=0; i<N-1; i++)
	{
		for(j=i+1; j<N; j++)
		{
			counter = P[i].coordinate - P[j].coordinate;
			min_img(counter, L);
//			dist2 = (P[i].diameter + P[j].diameter)/2;
			rij = counter.norm();

			for(k = 0; k < bondlist.size(); k++)
			{	//Means the two particles are bonded
				if(bondlist[k].partner1 == i && bondlist[k].partner2 == j)
				{
					particle_bonded = true;
					l1 = (1-bond_delta)*bondlist[k].bondlength;
					l2 = (1+bond_delta)*bondlist[k].bondlength;
					//This means everything is fine and dandy
					if(l1*(1-tolerance)<rij && l2*(1+tolerance)>rij)	
						{break;}
					else
					{
						overlap = true;
						cout<<"Particles overlapping which are bonded i ="<<i;
						cout<<"\t j = "<<j<<"\t l1 = "<<l1<<"\t l2 = "<<l2;
						cout<<"\t rij = "<<rij<<endl;
						cout<<"i = "<<i<<"\t x = "<<P[i].coordinate.x;
						cout<<"\t y =" <<P[i].coordinate.y<<"\t z = ";
						cout<<P[i].coordinate.z<<"\tcell_ID= "<<P[i].cell_ID<<"\tvx=";
						cout<<P[i].velocity.vx<<"\tvy="<<P[i].velocity.vy<<"\tvz=";
						cout<<P[i].velocity.vz<<endl;
						cout<<"j = "<<j<<"\t x = "<<P[j].coordinate.x;
						cout<<"\t y =" <<P[j].coordinate.y<<"\t z = ";
						cout<<P[j].coordinate.z<<"\tcell_ID= "<<P[j].cell_ID<<"\tvx=";
						cout<<P[j].velocity.vx<<"\tvy="<<P[j].velocity.vy<<"\tvz=";
						cout<<P[j].velocity.vz<<endl;
						cout<<"TIME="<<TIME<<"\tfp_TIME="<<fpupdate_TIME<<endl;
						return overlap;
					}
				}	
				else 
					{particle_bonded = false;}	
			}
			//The two particles are not bonded, looked through the whole bondlist
			if(!particle_bonded)			
			{
				for(k = 0; k<nonbondlist.size(); k++)
				{
					//Find the two particles in nonbondlist
					if((nonbondlist[k].partner1 == P[i].atomtype && nonbondlist[k].partner2 == P[j].atomtype) || (nonbondlist[k].partner1 == P[j].atomtype && nonbondlist[k].partner2 == P[i].atomtype) ) 
					{
						particle_nonbonded = true;
						rij1 = rij/nonbondlist[k].L1[0];
//overlapping, just trying to find their times so that we can figure out which event is meing missed out
						if(rij1 + tolerance < 1)
						{
							overlap = true;
							cout<<"Particles overlapping which are nonbonded i="<<i;
						   	cout<<"\t j= "<<j<<"\t nonbondlist[k].L1[0]=";
							cout<<nonbondlist[k].L1[0]<<"\t rij = "<<rij<<endl;
							cout<<"i = "<<i<<"\tx = "<<P[i].coordinate.x<<"\t y =";
							cout<<P[i].coordinate.y<<"\t z = "<<P[i].coordinate.z;
							cout<<"\t cell_ID = "<<P[i].cell_ID<<"\tvx=";
							cout<<P[i].velocity.vx<<"\tvy="<<P[i].velocity.vy;
							cout<<"\tvz="<<P[i].velocity.vz<<endl;
							cout<<"j = "<<j<<"\tx = "<<P[j].coordinate.x<<"\t y =";
							cout<<P[j].coordinate.y<<"\t z = "<<P[j].coordinate.z;
							cout<<"\t cell_ID = "<<P[j].cell_ID<<"\tvx=";
							cout<<P[j].velocity.vx<<"\tvy="<<P[j].velocity.vy;
							cout<<"\tvz="<<P[j].velocity.vz<<endl;
							cout<<"TIME="<<TIME<<"\tfp_TIME="<<fpupdate_TIME<<endl;
							//If once returning this, then its fine
							return overlap;	
						}
						break; 
					}
					else
						{particle_nonbonded = false;}
				}
				if(!particle_nonbonded)		//Particle pairing not found in nonbondedlist
				{
					cout<<"Error, particle pairing not found in nonbondlist while checking";
					cout<<" overlaps, weird stuff bro"<<endl;
					exit(1);
				}
			}
		}
	} 
	return overlap;
}
//Overlap checking function for initialising the system, because at t=0, not all coordinates are initialised
//Wont be needing this after all now
/*bool System::CheckOverlap_bead(int a, vector<Particle> &P)
{
	int beadno, chainno, i, j;
	bool overlap = false;
	double tolerance = 1.0e-4;
	XYZ counter;
	double rij, dist2;

	beadno = a%chainlength;
	chainno = int(double(a)/double(chainlength));
	for(i = 0; i< chainno; i++)//Checking if overlapping with a particle from another chain
	{
		for(j = 0; j < chainlength; j++)
		{
			counter = P[i*chainlength + j].coordinate - P[a].coordinate;
			min_img(counter, L);
			dist2 = (P[i].diameter + P[j].diameter)/2;
			rij = counter.norm2()/dist2;
			if(rij + tolerance < 1)
			{	
				overlap = true;
				cout<<"Particles overlapping (different chains) are i = "<<a<<"\t j = "<<i*chainlength+j<<endl;
				return overlap;	//If once returning this, then its fine
			}
		}
	}
	for(j = 0; j<beadno-1; j++)//Checking if overlapping with a particle from the same chain, except the neighbor
	{
		counter = P[chainno*chainlength + j].coordinate - P[a].coordinate;
		min_img(counter, L);
		dist2 = (P[i].diameter + P[j].diameter)/2;
		rij = counter.norm2()/dist2;
		if(rij + tolerance < 1)
		{
			overlap = true;
			cout<<"Particles overlapping (same chain) are i = "<<a<<"\t j = "<<chainno*chainlength+j<<endl;
			return overlap;	//If once returning this, then its fine
		}
	}
	//checking overlap between particle and its neighbor
	counter = P[chainno*chainlength + beadno -1].coordinate - P[a].coordinate;
	min_img(counter, L);
	dist2 = (P[i].diameter + P[j].diameter)/2;
	rij = counter.norm2()/dist2;
	if(rij + tolerance < 1)
	{
		overlap = true;
		cout<<"Particles overlapping (same chain) are i = "<<a<<"\t j = "<<chainno*chainlength+j<<endl;
		return overlap;	//If once returning this, then its fine
	}
	return overlap;
}*/
//Function to calculate the potential energy //Need to add multiple potential possibilities
double System::pe(vector<Particle> &P, int N)
{
	double pote = 0, r;
	int i, j, k, l;
	XYZ counter;
	cout<<"Size of nonbondlist = "<<nonbondlist.size()<<endl;
	for(i=0; i<N-1; i++)
	{
		for(j=i+1; j<N; j++)
		{
			for(k = 0; k < nonbondlist.size(); k++)
			{//This means we have found the particle in the nonbondlist
				if((nonbondlist[k].partner1 == P[i].atomtype && nonbondlist[k].partner2 == P[j].atomtype) || (nonbondlist[k].partner2 == P[i].atomtype && nonbondlist[k].partner1 == P[j].atomtype))
				{
//					cout<<"Entered for particle number i ="<<i<<"\t j = "<<j<<"\t atomtype[i] = "<<P[i].atomtype<<"\t atomtype[j] = "<<P[j].atomtype<<endl;
					counter = P[i].coordinate - P[j].coordinate;
					min_img(counter, L);
					r = counter.norm();
					l = 0;
//					cout<<"Ghuse hain for l = "<<l<<"\t nonbondlist[k].L1.size() = "<<nonbondlist[k].L1.size()<<endl;
					while(l < nonbondlist[k].L1.size())
					{
						if(r>nonbondlist[k].L1[l] && r<nonbondlist[k].L2[l]) 
						{
							pote = pote + nonbondlist[k].epsilon[l];
							break;
						}
//						else
//							{continue;}
						l = l+1;
					}
				}
//				else
//					{continue;}
			}
		}
	}
	return pote;
}
//Function to calculate the kinetic energy
double System::ke(vector<Particle> &P, int N)
{
	double kinetic_energy=0;
	for(int i=0; i<N; i++)
	{
		kinetic_energy = kinetic_energy + 0.5*P[i].mass*P[i].velocity2;
	}
	return kinetic_energy;
}
//Net momentum calculation function
VEL System::net_momentum(vector<Particle> &P, int N)
{
	VEL momentum;
	for(int i=0; i<N; i++)
	{
		momentum.vx = momentum.vx + P[i].mass*P[i].velocity.vx;
		momentum.vy = momentum.vy + P[i].mass*P[i].velocity.vy;
		momentum.vz = momentum.vz + P[i].mass*P[i].velocity.vz;
	}
	return momentum;
}
//to initalize the system when no run has been performed, or when a dcd file or lammps 
//file is present, pick up the last time and the last configuration of particles
void System::Create()	
{
//	N = nchains*chainlength; 	
//IF the filename_out is not added, then use the same name as pdb file
	if(filename_out == "test")
	{
		filename_out = filename_pdb;
		//Remove .pdb from the name
		filename_out.erase(filename_out.end()-4,filename_out.end());
//Inserts _ in front of the name from pdb
		filename_out.insert(0,"_");
	}
	std::string filename_restart;
//	std::ostringstream fnd;					//FileName dummy
//	fnd<<"_restart_nc_"<<nchains<<"_cl_"<<chainlength<<"_T_"<<temperature<<"_L_"<<L<<"_m_ "<<mass<<"_s1_"<<sigma1<<"_s2_"<<sigma2<<"_bl_"<<bondlength<<"_del_"<<bond_delta<<"_d_"<<well_depth<<".restart";
	filename_restart = filename_out;
	filename_restart.append(".restart");
//	std::string FileName = fnd.str();
//Defining a,b,c,theta, phi, psi
	cellparameters.push_back(L);
	cellparameters.push_back(L);
	cellparameters.push_back(L);
	cellparameters.push_back(90.0);
	cellparameters.push_back(90.0);
	cellparameters.push_back(90.0);
	
	ncell_max = int(L/cellwidth) - 1;
	cout<<"ncell_max = "<<ncell_max<<endl;
	if(ncell_max > 9)
	{
		cout<<"Exiting because ncell_max > 9"<<endl;
		exit(10);	
	}

	ReadFiles filereader;
	filereader.ReadPDBFile(filename_pdb, P, cellparameters);
	cout<<"Initialised the coordinates from pdb"<<endl;	
	
	filereader.ReadPSFFile(filename_psf, P, bondlist, nonbondlist);
	cout<<"Initialised the info from psf"<<endl;		

	filereader.ReadParameterFile(filename_par, P, bondlist, nonbondlist);
	cout<<"Initialised the info from par"<<endl;		

	N = P.size();
//	Checking if a previous particle readin exists. If yes, then start the simulation from that point
	if(std::ifstream(filename_restart))
	{
		std::cout<<"Restart file exists, will take values from the last file"<<endl;
		std::string dumps;
		int Npar;
		double time, length;
		double dumpi, dumpd;
		double x, y, z, vx, vy, vz;

		ifstream in;
		in.open(filename_restart);						
			//Checking if the data is for the same N and L
  			std::getline(in,dumps);
    			in>>time;
  			std::getline(in,dumps);
		   	in>>Npar;
  			std::getline(in,dumps);
  			std::getline(in,dumps);
    			in>>dumpd>>length;
		in.close();

		if(Npar == N && length == L)
		{
			cout<<"Restart file for same N and L, continue simulations, N= "<<N<<"\t L = "<<L<<endl;
		}
		TIME = time;
		Particle p;

		in.open(filename_restart);
			//reading x, y and z coordinates as well as velocities
  			std::getline(in,dumps);
  			std::getline(in,dumps);
    			std::getline(in,dumps);					
  			std::getline(in,dumps);
  			std::getline(in,dumps);
  			std::getline(in,dumps);
  			std::getline(in,dumps);
  			std::getline(in,dumps);
  			std::getline(in,dumps);
			while (!in.eof())
			{
				in>>dumpi>>x>>y>>z>>vx>>vy>>vz;
				P[int(dumpi)].coordinate.x = x;
				P[int(dumpi)].coordinate.y = y;
				P[int(dumpi)].coordinate.z = z;
				P[int(dumpi)].velocity.vx = vx;
				P[int(dumpi)].velocity.vy = vy;
				P[int(dumpi)].velocity.vz = vz;
				P[int(dumpi)].velocity2 = vx*vx + vy*vy + vz*vz;
			//	P.push_back(p);			
			}
		in.close();
		for(int i = 0; i<N; i++)
		{
			cout<<"Restart X = "<<P[i].coordinate.x<<"\t Y = "<<P[i].coordinate.y<<"\t Z = ";
			cout<<P[i].coordinate.z<<"\t VX = "<<P[i].velocity.vx<<"\t VY = "<<P[i].velocity.vy;
			cout<<"\t VZ = "<<P[i].velocity.vz<<endl;
		}
	}
	else//Means file doesn't exist
	{
 		cout<<"Creating system during simulation"<<endl;
//		double gridsize;//width of each grid 
//		cout<<"N= "<<N<<"\t L= "<<L<<"\t nsweep= "<<nsweep<<endl;	
//Allocate Grid
//		Nchain_perside = int(pow(nchains,0.33333333))+1;	
//		gridsize = double(L)/double(Nchain_perside);
//		cout<<"Grid width="<<gridsize<<"\t No. of chains per side= "<<Nchain_perside<<endl;
//		coordinate_initialization(gridsize, Nchain_perside);//Initialising time

		velocity_initialization();//Initialising velocities
		cout<<"Initialised the velocities"<<endl;
		TIME = 0.0;							//initialising time
		fpupdate_TIME=0.0;
	}
//Couting for testing the correct reading of the system
	cout<<"The read bondlist is here: "<<endl;
	for (int i=0; i<bondlist.size(); i++)
	{
		cout<<"partner1= "<<bondlist[i].partner1<<"\t partner2= "<<bondlist[i].partner2;
		cout<<"\t middle= "<<bondlist[i].middle<<"\t bondlength= "<<bondlist[i].bondlength<<endl;
	}

	cout<<"The read nonbondlist is here: "<<endl;
	for (int i=0; i<nonbondlist.size(); i++)
	{
		cout<<"partner1= "<<nonbondlist[i].partner1<<"\t partner2= "<<nonbondlist[i].partner2;
		for(int j=0; j<nonbondlist[i].L1.size(); j++)
		{
			cout<<"\t L1= "<<nonbondlist[i].L1[j]<<"\t L2= "<<nonbondlist[i].L2[j];
			cout<<"\t epsilon= "<<nonbondlist[i].epsilon[j];
		}
		cout<<endl;
	}

	if(CheckOverlap(P))
	{
		cout<<"Overlap of particles at initialisation, so ending simulation"<<endl;
		exit(1);
	}
//	cout<<"overlap checking complete"<<endl;
	//Using this is useful when different pairs have different well width
	//No need to consider bonds because usually the rnage of nonbonded is longer anyway
	max_well_width = 0;
//	for(int i = 0; i < bondlist.size(); i++)
//	{
//		if(bondlist[i].bondlength*(1+bond_delta) > max_well_width)
//			{max_well_width = bondlist[i].bondlength*(1+bond_delta);}
//	}
	//Since we are ordering them in increasing order of lengths, this makes sense
	for(int i = 0; i< nonbondlist.size(); i++)
	{
		if(nonbondlist[i].L2[nonbondlist[i].L2.size()-1] > max_well_width)
			{max_well_width = nonbondlist[i].L2[nonbondlist[i].L2.size()-1];}
	}
	cout<<"Max well width = "<<max_well_width<<endl;
	//Minimum number of cells on each side so that the calculation of cell-list is beneficial = 4
	if(4*cellwidth > L && cellwidth > max_well_width)
	{
		cout<<"Will not be using celllist, still running simulation"<<endl;
		celllist_counter = false;
	}
	else
	{
		cout<<"Using CellList"<<endl;
		celllist_counter = true;
		CellList_initialization();
		CellList();

		neighborlist_range = neighborlist_factor*max_well_width;
		if(neighborlist_range < cellwidth)
		{
			cout<<"Will use neighborlist with celllist"<<endl;
			neighborlist_counter = true;
			NeighborList_initialization();
			NeighborList();
		}
		else
		{
			cout<<"Using celllist but not neighborlist"<<endl;
			neighborlist_counter = false;
		}
	}
	potential_energy = pe(P, N);
	cout<<"Potential energy = "<<potential_energy<<endl;
}
//Initializing the coordinate, when initializing within the code
//Not needed anymore, because using pdb
/*void System::coordinate_initialization(double size_of_grid, int chains_per_side)
{
	Particle p;
	int i,j,k,n, current_chain, bead_overlap[chainlength], chain_overlap[nchains];
	XYZ rand;
	std::random_device rd;  						//Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); 						//Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dist1(-0.5, 0.5);
	std::uniform_real_distribution<> dist2(-0.5, 0.5);			//For randomly placing the particles
	std::uniform_real_distribution<> dist3(-0.5, 0.5);

	for(i = 0; i<N; i++)
	{
		p.coordinate.x = 0;
		p.coordinate.y = 0;
		p.coordinate.z = 0;						//Initialising vectors for position and velocity
		p.velocity.vx = 0;
		p.velocity.vy = 0;
		p.velocity.vz = 0;
		p.velocity2 = 0;
		P.push_back(p);
	}
	//initialising the first bead of all the chains
	Grid g;
	for(k=0; k<chains_per_side; k++)
	{
		for(j=0; j<chains_per_side; j++)
		{
			for(i=0; i<chains_per_side; i++)
       			{
       				g.com.x=(double(i)+0.5)*size_of_grid;
       				g.com.y=(double(j)+0.5)*size_of_grid;
        			g.com.z=(double(k)+0.5)*size_of_grid;
				g.velocity_com.vx = 0;
				g.velocity_com.vy = 0;
				g.velocity_com.vz = 0;
	       			G.push_back(g);
	       		}
	       	}
	}
	//Initialising first beads of the chain
	for(i = 0; i<nchains; i++)
	{
		P[i*chainlength].coordinate = G[i].com;
		chain_overlap[i] = 0;
	}
	//Now, to initialize all the other beads
	for(current_chain = 0; current_chain<nchains; current_chain++)
	{
		n = current_chain*chainlength;
		for(i = 0; i<chainlength; i++)
		{
			bead_overlap[i] = 0;				//Counter for overlap in the current chain, first particle of each chain already initialised
		}
		for(i = 1; i<chainlength; i++)
		{	
			rand.x = dist1(gen);
			rand.y = dist2(gen);
			rand.z = dist3(gen);
			rand = rand*bondlength/rand.norm();			//Centering the new particle at the distance of bondlength

			P[n + i].coordinate = P[n + i - 1].coordinate + rand;
			PBC(P[n + i].coordinate.x, P[n + i].coordinate.y, P[n + i].coordinate.z);

		 	if(CheckOverlap_bead(n + i, P))	//overlap occured when we put the next bead
			{
				bead_overlap[i]  = bead_overlap[i] + 1;
				cout<<"current_chain = "<<current_chain<<"\t i = "<<i<<"bead_overlap[i] = "<<bead_overlap[i]<<endl;
				if(bead_overlap[i] > 1000)		//Configurational space around the particle is filled, better to restart chain
				{
					chain_overlap[current_chain] = chain_overlap[current_chain] + 1;
					cout<<"current_chain = "<<current_chain<<"chain_overlap[current_chain] = "<<chain_overlap[current_chain]<<endl;
					if(chain_overlap[current_chain] > 500)
					{
						sigma2 = 0.999*sigma2;
						bondlength = 0.999*bondlength;
						for(j = 0; j<N; j++)
						{
							P[j].diameter = P[j].diameter*0.999;
							P[j].coordinate = P[j].coordinate*0.999;
							PBC(P[j].coordinate.x, P[j].coordinate.y, P[j].coordinate.z);
						}
						chain_overlap[current_chain] = 0;
					}
					current_chain = current_chain - 1;
					break;					
				}
				i = i-1;				//if bead_overlap <1000, then it starts to fill it again, if not, then it doesnt matter, because we are reinitialising the chain
			}
		}
	}
}*/
//Velocity initialization after the coordinates are fixed
void System::velocity_initialization()
{
//c.o.m. velocity in each direction, so that com velocity can be zero
	VEL avg_center_of_mass_momentum, avg_com_momentum_new, momentum_counter;				
	double avg_com_mv2 = 0.0, net_mass = 0.0;
	double max_particle_vel = 0.0;
	cout<<"N = "<<N<<endl;
//For maxwell boltzmann distribution of molecules
//	std::random_device rd{};					
//   	std::mt19937 gen{rd()};

//	std::mt19937 gen{0};
//	std::normal_distribution<> d1{0,1};				//0 is mean, 1 is SD
//	std::normal_distribution<> d2{0,1};
//	std::normal_distribution<> d3{0,1};
//	std::uniform_real_distribution<> d1{-1,1};
//	std::uniform_real_distribution<> d2{-1,1};
//	std::uniform_real_distribution<> d3{-1,1};
	
	for (int i = 0; i<N; i++)
	{
		P[i].velocity.vx = RandomGaussianNumber();
		P[i].velocity.vy = RandomGaussianNumber();
		P[i].velocity.vz = RandomGaussianNumber();
	
//	for (int i = 0; i<N; i++)
//	{
//		P[i].velocity.vx = (d1(gen)-0.5)*sqrt(temperature/P[i].mass);
//		P[i].velocity.vy = (d2(gen)-0.5)*sqrt(temperature/P[i].mass);
//		P[i].velocity.vz = (d3(gen)-0.5)*sqrt(temperature/P[i].mass);
//		cout<<"P[i].velocity.vx = "<<P[i].velocity.vx<<"\t P[i].velocity.vy = "<<P[i].velocity.vy<<"\t P[i].velocity.vz = "<<P[i].velocity.vz<<endl;
    
		momentum_counter=P[i].velocity*P[i].mass;
		avg_center_of_mass_momentum = avg_center_of_mass_momentum + momentum_counter;
		avg_com_mv2 = avg_com_mv2 + P[i].velocity.norm2()*P[i].mass;
		net_mass = net_mass + P[i].mass;
		P[i].velocity2 = P[i].velocity.norm2();
	}
//This makes it average velocity
	avg_center_of_mass_momentum = avg_center_of_mass_momentum/net_mass;
	avg_com_mv2 = avg_com_mv2/net_mass;
	cout<<"COM_vx = "<<avg_center_of_mass_momentum.vx<<"\tCOM_vy = "<<avg_center_of_mass_momentum.vy<<"\tCOM_vz = "<<avg_center_of_mass_momentum.vz<<"\tCOM_v2 = "<<avg_com_mv2<<endl;
	for(int i=0; i<N; i++)
	{
		//Making com velocity zero and scaling to give required thermal velocities
		momentum_counter=avg_center_of_mass_momentum*P[i].mass;

		P[i].velocity = (P[i].velocity-avg_center_of_mass_momentum)*sqrt(temperature/P[i].mass);	
		cout<<"P[i].velocity.vx = "<<P[i].velocity.vx<<"\t P[i].velocity.vy = "<<P[i].velocity.vy<<"\t P[i].velocity.vz = "<<P[i].velocity.vz<<endl;
		P[i].velocity2 = P[i].velocity.norm2();
		avg_com_mv2 = avg_com_mv2 + P[i].velocity2*P[i].mass;
//		P[i].velocity = (P[i].velocity-momentum_counter);//*sqrt(temperature/P[i].mass);	
//		cout<<"P[i].velocity.vx = "<<P[i].velocity.vx<<"\t P[i].velocity.vy = "<<P[i].velocity.vy<<"\t P[i].velocity.vz = "<<P[i].velocity.vz<<endl;
//		P[i].velocity2 = P[i].velocity.norm2();
//		avg_com_mv2 = avg_com_mv2 + P[i].velocity2*P[i].mass;

//New momentum after updating particle velocity
		momentum_counter=P[i].velocity*P[i].mass;
		avg_com_momentum_new = avg_com_momentum_new + momentum_counter; 
	}

	avg_com_momentum_new = avg_com_momentum_new/net_mass;
	avg_com_mv2 = avg_com_mv2/net_mass;
	cout<<"COM_vx = "<<avg_com_momentum_new.vx<<"\tCOM_vy = "<<avg_com_momentum_new.vy<<"\tCOM_vz = "<<avg_com_momentum_new.vz<<"\tCOM_v2 = "<<avg_com_mv2<<endl;
	//Defining maxvel for neighborlist calculation
	for(int i = 0; i<N; i++)
	{
		if(P[i].velocity2 > max_particle_vel)
		{
			max_particle_vel = P[i].velocity2;
		}
	}
	if(10*sqrt(max_particle_vel)> maxvel)	
		{maxvel = 10*sqrt(max_particle_vel);}
}
//Cell list vectors initialization and length fixedness
void System::CellList_initialization()
{
	for(int i =0; i<N; i++)
		{LinkList.push_back(-1);} 
	for(int i =0; i<1000; i++)
		{HeadList.push_back(-1);}
}
//Cell calculation for various particles
void System::CellList()
{
	int cell_x, cell_y, cell_z;
	double remainder_x, remainder_y, remainder_z;
	int i;
	XYZ particler;
	//Reinitialization of linklist and headlist
	for(i =0; i<N; i++)
	{
		LinkList[i] = -1;
	}
	for(i =0; i<1000; i++)
	{
		HeadList[i] = -1;
	}
//Cell Assignment
	for(i = 0; i < N; i++)
	{
//		PBC(P[i].coordinate.x, P[i].coordinate.y, P[i].coordinate.z);	
		particler = OneParticlePositionUpdater(P[i], TIME, fpupdate_TIME);	
//		PBC(P[i].coordinate.x, P[i].coordinate.y, P[i].coordinate.z);
		cell_x = int(particler.x/cellwidth);
		cell_y = int(particler.y/cellwidth);
		cell_z = int(particler.z/cellwidth);

		if(cell_x < 0)
			{cell_x = ncell_max;}
		if(cell_y < 0)
			{cell_y = ncell_max;}
		if(cell_z < 0)
			{cell_z = ncell_max;}
		if(cell_x > ncell_max)
			{cell_x = 0;}
		if(cell_y > ncell_max)
			{cell_y = 0;}
		if(cell_z > ncell_max)
			{cell_z = 0;}

		P[i].cell_ID = 100*cell_z + 10*cell_y + cell_x;
		LinkList[i] = HeadList[P[i].cell_ID];
		HeadList[P[i].cell_ID] = i;
	}
}
//Neighbor list vector initialization and length fixedness
void System::NeighborList_initialization()
{
	for(int i =0; i<100*N; i++)
	{
		NeighborLinkList.push_back(-1);
	} 
	for(int i =0; i<N; i++)
	{
		NeighborHeadList.push_back(-1);
	}	
}
//Use the already formed celllist to find the neighbors for each particle
void System::NeighborList()
{
	int j, cell_x, cell_y, cell_z, cell_xcalc, cell_ycalc, cell_zcalc;
	int m1, m2, m3, cell_ID;
	int x_l, x_u, y_l, y_u, z_l, z_u;
	int nbrlinklist_counter = 0;
	Particle difference;
	double r;
//Reinitialization of the neighbor headlist and linklist	
	for(int i = 0; i< N; i++)
	{
		NeighborHeadList[i] = -1;
	}
	for(int i = 0; i<100*N; i++)
	{
		NeighborLinkList[i] = -1;
	}
	for(int i = 0; i < N; i++)
	{
		//Extracting cell_ID
		cell_x = P[i].cell_ID%10;
		cell_y = int(double(P[i].cell_ID%100)/10);
		cell_z = int(double(P[i].cell_ID)/100);
		//lower and upper neighboring cells
		x_l = cell_x-1;
		y_l = cell_y-1;
		z_l = cell_z-1;
		x_u = cell_x+1;
		y_u = cell_y+1;
		z_u = cell_z+1;
		for(m1 = x_l; m1<=x_u; m1++)
		{	
			for(m2 = y_l; m2<=y_u; m2++)
			{
				for(m3 = z_l; m3<=z_u; m3++)
				{
					cell_xcalc = m1;
					cell_ycalc = m2;
					cell_zcalc = m3;
					//PBC to find neighboring cells
					if(cell_xcalc < 0)
						{cell_xcalc = ncell_max;}
					if(cell_ycalc < 0)
						{cell_ycalc = ncell_max;}
					if(cell_zcalc < 0)
						{cell_zcalc = ncell_max;}
					if(cell_xcalc > ncell_max)
						{cell_xcalc = 0;}
					if(cell_ycalc > ncell_max)
						{cell_ycalc = 0;}
					if(cell_zcalc > ncell_max)
						{cell_zcalc = 0;}
					cell_ID = 100*cell_zcalc + 10*cell_ycalc + cell_xcalc;
					j = HeadList[cell_ID];
					while(j>=0)
					{
						if(j != i)
						{
							difference.coordinate = P[i].coordinate-P[j].coordinate;
							difference.velocity = P[i].velocity-P[j].velocity;
							difference.coordinate=OneParticlePositionUpdater(difference,TIME,fpupdate_TIME);
							min_img(difference.coordinate, L);
							r = difference.coordinate.norm();
							if(r <= neighborlist_range)
							{
								NeighborLinkList[nbrlinklist_counter] = j;
								NeighborHeadList[i] = nbrlinklist_counter;
								nbrlinklist_counter = nbrlinklist_counter + 1;
							}
						}
						j = LinkList[j];
					}

				}
			}
		}
		NeighborLinkList[nbrlinklist_counter] = -1;
		nbrlinklist_counter = nbrlinklist_counter + 1;
	}
}
XYZ System::OneParticlePositionUpdater(Particle &A, double TIME_current, double TIME_fpupdated)
{
	XYZ returner;
	returner.x=A.coordinate.x + A.velocity.vx*(TIME_current-TIME_fpupdated);
	returner.y=A.coordinate.y + A.velocity.vy*(TIME_current-TIME_fpupdated);
	returner.z=A.coordinate.z + A.velocity.vz*(TIME_current-TIME_fpupdated);
//	min_img(returner, L);
	PBC(returner.x, returner.y, returner.z);
	return returner;
}
XYZ System::OneParticlePositionBackwarder(Particle &A, double TIME_current, double TIME_fpupdated)
{
	XYZ returner;
	returner.x=A.coordinate.x - A.velocity.vx*(TIME_current-TIME_fpupdated);
	returner.y=A.coordinate.y - A.velocity.vy*(TIME_current-TIME_fpupdated);
	returner.z=A.coordinate.z - A.velocity.vz*(TIME_current-TIME_fpupdated);
//	min_img(returner, L);
	PBC(returner.x, returner.y, returner.z);
	return returner;
}
