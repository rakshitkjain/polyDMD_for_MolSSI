//Function definitions of system class (Revision Date: 31st March 2021)
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
//Generate gaussian number for velocity initialization and Andersen Thermostat
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
			rij = counter.norm();

			for(k = 0; k < bondlist.size(); k++)
			{	//Means the two particles are bonded
				if(bondlist[k].partner1 == i && bondlist[k].partner2 == j)
				{
					particle_bonded = true;
					l1 = (1-bond_delta)*bondlist[k].bondlength;
					l2 = (1+bond_delta)*bondlist[k].bondlength;
					//This means everything is fine
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
						//Overlapping, just trying to find their times so that we can figure out which event is being missed out
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
							return overlap;	
						}
						break; 
					}
					else
						{particle_nonbonded = false;}
				}
				//Particle pairing not found in nonbondedlist
				if(!particle_nonbonded)		
				{
					cout<<"Error, particle pairing not found in nonbondlist while checking";
					cout<<" overlaps, very weird"<<endl;
					exit(1);
				}
			}
		}
	} 
	return overlap;
}
//Function to calculate the potential energy 
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
					counter = P[i].coordinate - P[j].coordinate;
					min_img(counter, L);
					r = counter.norm();
					l = 0;
					while(l < nonbondlist[k].L1.size())
					{
						if(r>nonbondlist[k].L1[l] && r<nonbondlist[k].L2[l]) 
						{
							pote = pote + nonbondlist[k].epsilon[l];
							break;
						}
						l = l+1;
					}
				}
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
//To initalize the system when no run has been performed, when a dcd file or lammps 
//file is present, pick up the last time and the last configuration of particles
void System::Create()	
{
//If the filename_out is not added, then use the same name as pdb file
	if(filename_out == "test")
	{
		filename_out = filename_pdb;
		//Remove .pdb from the name
		filename_out.erase(filename_out.end()-4,filename_out.end());
		//Inserts _ in front of the name from pdb
		filename_out.insert(0,"_");
	}
	std::string filename_restart;
	filename_restart = filename_out;
	filename_restart.append(".restart");
//Defining a,b,c,theta, phi, psi
	cellparameters.push_back(L);
	cellparameters.push_back(L);
	cellparameters.push_back(L);
	cellparameters.push_back(90.0);
	cellparameters.push_back(90.0);
	cellparameters.push_back(90.0);
//Maximum number of cells in the system for cell list calculation
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
//Checking if a previous particle readin exists. If yes, then start the simulation from that point
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
 		cout<<"Starting new simulation here"<<endl;

		velocity_initialization();			//Initialising velocities
		cout<<"Initialised the velocities"<<endl;
		TIME = 0.0;					//Initialising time
		fpupdate_TIME=0.0;
	}
//Printing to make sure that the input is being read correctly
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
	max_well_width = 0;
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

//Velocity initialization after the coordinates are fixed
void System::velocity_initialization()
{
//c.o.m. velocity in each direction, so that com velocity can be zero
	VEL avg_center_of_mass_momentum, avg_com_momentum_new, momentum_counter;				
	double avg_com_mv2 = 0.0, net_mass = 0.0;
	double max_particle_vel = 0.0;
	cout<<"N = "<<N<<endl;
	
	for (int i = 0; i<N; i++)
	{
		P[i].velocity.vx = RandomGaussianNumber();
		P[i].velocity.vy = RandomGaussianNumber();
		P[i].velocity.vz = RandomGaussianNumber();
	    
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
//Cell list vectors initialization and fixing lengths
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
		particler = OneParticlePositionUpdater(P[i], TIME, fpupdate_TIME);	
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
