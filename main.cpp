//DMD: MAIN CODE (Revision Date: 7th Jan 2020)
#include "header.h"
#include "system.h"
#include "timecalc.h"
#include "collision.h"
#include "rdf.h"
#include "writefiles.h"
#include "readinput.h"
#include "readfiles.h"

int main(int argc, char *argv[])
{
	int a,i, f, col_counter, ghost_counter = 0, fpupdate_counter = 0;
	double tmin, tmin1, start_time;
//Initial and Final kinetic and total energy
	double eni_ke = 0, enf_ke = 0, eni_te = 0, enf_te = 0;			
//Compressibility factor, rate of collisions, volume fractions
	double comp_fac = 0, rate = 0, vol_frac = 0;				
//Variables for calculations during run
	double kineticenergy, totalenergy, burn_time;				
	int thermostat_particle;
	VEL momentum;
	time_t initialization_time = time(NULL);
	time_t percent_time;

	ReadInput reader;
	WriteFiles writer;
	vector<string> arguments;
	for(int i = 0; i< argc; i++)
		{arguments.push_back(string(argv[i]));}
//Read Input using command line parsing
	reader.ReadVariables(argc, arguments);
//Create system
	reader.initialize.Create();	

	Collision C;
	RDF rdf;
	C.TC.S = reader.initialize;
	writer.S = reader.initialize;
//Kinetic energy and total energy functions (Potential energy is updated within the collision segment)
	eni_ke = C.TC.S.ke(C.TC.S.P, C.TC.S.N);
	eni_te = eni_ke + C.TC.S.potential_energy;								
	momentum = C.TC.S.net_momentum(C.TC.S.P, C.TC.S.N);
//Parameter time variation file initialization
	writer.TimeVarFileIni(C.TC.S.filename_out);
//Parameters time variation file printing
	writer.TimeVarFile(C.TC.S.TIME, eni_ke, C.TC.S.potential_energy, eni_te, momentum, C.TC.S.filename_out);
//initialising RDF
	rdf.RDF_ini(C.TC.S.L, C.TC.S.maxbin);
//Write dump, PDB and binary file for particle coords
	writer.WriteDump(0.0, C.TC.S.P, C.TC.S.N, C.TC.S.L, C.TC.S.filename_out);
	writer.PDBFile(C.TC.S.P, C.TC.S.filename_out);	
	writer.DCDWriteStep(C.TC.S.P, C.TC.S.N, C.TC.S.filename_out);

//Useful to check if the file already exists or not
	ofstream out;
                std::ostringstream fnd;//FileName dummy
                fnd<<C.TC.S.filename_out<<".txt";
                std::string FileName = fnd.str();                                                      
//Checking initialisation overlap
	if(C.TC.S.CheckOverlap(C.TC.S.P))
	{
		cout<<"Overlap of particles inside the simulation run, so ending simulation at TIME= "<<C.TC.S.TIME<<endl;
		goto exit;
	}
//Initialising time array
	C.TC.timearray_initialization(C.TC.S.N);

//If using celllist or neighborlist or not, then making the first timelist
	if(C.TC.S.celllist_counter)
	{
		if(C.TC.S.neighborlist_counter)
			{C.TC.Collision_time_neighborlist(C.TC.S.L, C.TC.S.N);}
		else
			{C.TC.Collision_time_celllist(C.TC.S.L, C.TC.S.N);}
	}
	else
		{C.TC.Collision_time(C.TC.S.L, C.TC.S.N);}

	start_time = C.TC.S.TIME;
//Main collision loop starts
	for (a=1; int(double(a)/double(C.TC.S.N))<C.TC.S.nsweep; a++)
	{
		cout<<"Entered loop"<<endl;
		//The system has a thermostat applied
		if(C.TC.S.andersen_freq != 0)
		{
			//Apply the thermostat in this step
			if(a%C.TC.S.andersen_freq==0 && int(double(a)/double(C.TC.S.N))>=C.TC.S.burn_percent*C.TC.S.nsweep)
			{
				//Thermostat, returns the particle on which it is applied
				thermostat_particle = C.AndersenThermostat();
				cout<<"The thermostat_particle is: "<<thermostat_particle<<endl;
				//When thermostat, count as ghost collision					
				ghost_counter = ghost_counter +1;
				//Now you have to re-form the timelist of the concerned particle(s)
				//Using celllist but not neighborlist
				if(C.TC.S.celllist_counter && (!C.TC.S.neighborlist_counter))
				{
					//Updating timelist after collision
					for (int i=0; i<C.TC.S.N; i++)
					{
						if(i==thermostat_particle || C.TC.Partner[i]==thermostat_particle)
						{				
							C.TC.UpdateUpCellList(i, C.TC.S.N);
						}
					}	
					C.TC.UpdateDownCellList(thermostat_particle, C.TC.S.N);
				}
				//Using both celllist and neighborlist, but just updating neighborlist
				else if(C.TC.S.celllist_counter && C.TC.S.neighborlist_counter)
				{
					//Updating timelist after collision
					for (int i=0; i<C.TC.S.N; i++)
					{
						if(i==thermostat_particle||C.TC.Partner[i]==thermostat_particle)
						{				
							C.TC.UpdateUpNeighborList(i, C.TC.S.N);
						}
					}	
					C.TC.UpdateDownNeighborList(thermostat_particle, C.TC.S.N);
				}
				//Neither celllist nor neighborlist used
				else if(!(C.TC.S.celllist_counter) && !(C.TC.S.neighborlist_counter))
				{
/*					for (i=0; i<C.TC.S.N; i++)
					{
						if(i==C.TC.collider)
							{cout<<"i for C.TC.collider="<<C.TC.collider<<endl;}
						if(C.TC.Partner[i]==C.TC.collider)
							{cout<<"i for C.TC.Partner[i]="<<i<<endl;}
						if(i==C.TC.col_partner)
							{cout<<"i for C.TC.col_partner="<<i<<endl;}
						if(C.TC.Partner[i]==C.TC.col_partner)
							{cout<<"i for C.TC.Partner[i]=C.TC.col_partner="<<i<<endl;}
					}
					//Updating timelist after collision*/
					for (int i=0; i<C.TC.S.N; i++)
					{
						if(i==thermostat_particle||C.TC.Partner[i]==thermostat_particle)
						{				
							C.TC.UpdateUpList(i, C.TC.S.N);
						}
					}	
					C.TC.UpdateDownList(thermostat_particle, C.TC.S.N);
				}
				//Finding smallest collision time after updating timelist
				tmin = C.TC.timbig;
				col_counter = C.TC.S.N;
				for(i = 0; i<C.TC.S.N; i++)
				{
					if(C.TC.TimeList[i]<tmin)
					{
						tmin = C.TC.TimeList[i];	
						col_counter = i;
					}
				}
				C.TC.minimum_time = tmin;
				C.TC.collider = col_counter;
			}
		}
		//Collision happens here
		C.Bump(C.TC.S.L, C.TC.S.N);
		//Checking overlap after every update of positions of all the particles
		if(a%C.TC.S.fpupdate_freq==0)
		{
			C.AllParticlePositionUpdater(C.TC.S.N, C.TC.S.P, C.TC.S.TIME, C.TC.S.fpupdate_TIME);
			if(C.TC.S.CheckOverlap(C.TC.S.P))					
			{
				cout<<"Overlap of particles inside the simulation run, so ending simulation at run number ="<<a<<"and TIME= "<<C.TC.S.TIME<<endl;
				goto exit;
			}

		//Rewriting fpupdate_TIME to the last time when system was updated
			C.TC.S.fpupdate_TIME=C.TC.S.TIME;
			fpupdate_counter=fpupdate_counter+1;
		//Writing restart, datafiles and updating RDF everytime when particles at actual positions
			writer.WriteDump(C.TC.S.TIME, C.TC.S.P, C.TC.S.N, C.TC.S.L, C.TC.S.filename_out);
			writer.RestartFile(C.TC.S.TIME, C.TC.S.P, C.TC.S.N, C.TC.S.filename_out);
			writer.DCDWriteStep(C.TC.S.P, C.TC.S.N, C.TC.S.filename_out);
			rdf.RDF_r(C.TC.S.L, C.TC.S.N, C.TC.S.P, C.TC.S.maxbin);
		}
		//If only cell list, use cell list functions

		if(C.TC.S.celllist_counter && (!C.TC.S.neighborlist_counter))
		{
			if(C.didcellchange)
			{
				//Updating timelist after cell change
				for (i=0; i<C.TC.S.N; i++)
				{
					if(i == C.TC.collider || C.TC.Partner[i] == C.TC.collider)
					{
						C.TC.UpdateUpCellList(i, C.TC.S.N);
					}
				}
				C.TC.UpdateDownCellList(C.TC.collider, C.TC.S.N);
			}
			else
			{
				//Updating timelist after collision
				for (i=0; i<C.TC.S.N; i++)
				{//Tells to update the list for collider and their listed partners
				if(i == C.TC.collider || C.TC.Partner[i] == C.TC.collider || i == C.TC.col_partner || C.TC.Partner[i] == C.TC.col_partner)
					{
						C.TC.UpdateUpCellList(i, C.TC.S.N);
					}
				}
				C.TC.UpdateDownCellList(C.TC.collider, C.TC.S.N);
				C.TC.UpdateDownCellList(C.TC.col_partner, C.TC.S.N);
			}
			//Finding smallest collision time after updating timelist
			tmin = C.TC.timbig;
			//Reassign colcounter
			col_counter = C.TC.S.N;
			//Finding the min time from the new timelist
			for(i = 0; i<C.TC.S.N; i++)
			{
				if(C.TC.TimeList[i]<tmin)
				{
					tmin = C.TC.TimeList[i];	
					col_counter = i;
				}
			}
			C.TC.minimum_time = tmin;
			C.TC.collider = col_counter;				
		}
		//IF both cell list and neighbor list, use neighborlist functions
		else if(C.TC.S.celllist_counter && C.TC.S.neighborlist_counter)
		{
			if(C.didnbrchange)
			{
				C.TC.Collision_time_neighborlist(C.TC.S.L, C.TC.S.N);
			}
			else 
			{
				if(C.didcellchange)
				{
					//Updating timelist after cell change
					for (i=0; i<C.TC.S.N; i++)
					{
						if(i == C.TC.collider || C.TC.Partner[i] == C.TC.collider)
						{
							C.TC.UpdateUpNeighborList(i, C.TC.S.N);
						}
					}
					C.TC.UpdateDownNeighborList(C.TC.collider, C.TC.S.N);
				}
				else
				{
					//Updating timelist after collision
					for (i=0; i<C.TC.S.N; i++)
					{
						if(i == C.TC.collider || C.TC.Partner[i] == C.TC.collider || i == C.TC.col_partner || C.TC.Partner[i] == C.TC.col_partner)
						{
							C.TC.UpdateUpNeighborList(i, C.TC.S.N);
						}
					}
					C.TC.UpdateDownNeighborList(C.TC.collider, C.TC.S.N);
					C.TC.UpdateDownNeighborList(C.TC.col_partner, C.TC.S.N);
				}	
				//Finding smallest collision time after updating timelist
				tmin = C.TC.timbig;
				col_counter = C.TC.S.N;
				for(i = 0; i<C.TC.S.N; i++)
				{
					if(C.TC.TimeList[i]<tmin)
					{
						tmin = C.TC.TimeList[i];	
						col_counter = i;
					}
				}
				C.TC.minimum_time = tmin;
				C.TC.collider = col_counter;
			}
		}
		//neither celllist not neighborlist used
		else if(!(C.TC.S.celllist_counter) && !(C.TC.S.neighborlist_counter))
		{
/*			for (i=0; i<C.TC.S.N; i++)
			{
				if(i==C.TC.collider)
				{cout<<"i for C.TC.collider="<<C.TC.collider<<endl;}
				if(C.TC.Partner[i]==C.TC.collider)
				{cout<<"i for C.TC.Partner[i]="<<i<<endl;}
				if(i==C.TC.col_partner)
				{cout<<"i for C.TC.col_partner="<<i<<endl;}
				if(C.TC.Partner[i]==C.TC.col_partner)
				{cout<<"i for C.TC.Partner[i]=C.TC.col_partner="<<i<<endl;}
			}
			//Updating timelist after collision*/
			for (i=0; i<C.TC.S.N; i++)
			{
				if(i == C.TC.collider || C.TC.Partner[i] == C.TC.collider || i == C.TC.col_partner || C.TC.Partner[i] == C.TC.col_partner)
				{
					C.TC.UpdateUpList(i, C.TC.S.N);
				}
			}
			C.TC.UpdateDownList(C.TC.collider, C.TC.S.N);
			C.TC.UpdateDownList(C.TC.col_partner, C.TC.S.N);
			//Finding smallest collision time after updating timelist
			tmin = C.TC.timbig;
			col_counter = C.TC.S.N;
			for(i = 0; i<C.TC.S.N; i++)
			{
				if(C.TC.TimeList[i]<tmin)
				{
					tmin = C.TC.TimeList[i];	
					col_counter = i;
				}
			}
			C.TC.minimum_time = tmin;
			C.TC.collider = col_counter;
		}

		if(C.TC.collider == C.TC.S.N)
		{
			cout<<"Crashing due to wrong assigment of the particle being moved"<<endl;
			exit(1);
		}

		//Data printed at every update frequency
		if(a%int(C.TC.S.print_freq*C.TC.S.fpupdate_freq) == 0)
		{
		//Calculation of KE, PE, TE, momentum
			kineticenergy = C.TC.S.ke(C.TC.S.P, C.TC.S.N);
			totalenergy = kineticenergy + C.TC.S.potential_energy;		
			momentum = C.TC.S.net_momentum(C.TC.S.P, C.TC.S.N);						
		//Time variation file printing
			writer.TimeVarFile(C.TC.S.TIME, kineticenergy, C.TC.S.potential_energy, totalenergy, momentum, C.TC.S.filename_out);
		}


//		if(a%50 == 0 && !(C.TC.S.celllist_counter) && !(C.TC.S.neighborlist_counter))
//		{			
//				{C.TC.Collision_time(C.TC.S.L, C.TC.S.N);}				
//		}
		cout<<"Collider= "<<C.TC.collider<<"\t Partner= "<<C.TC.Partner[C.TC.collider];
		cout<<"\t Time for this= "<<C.TC.S.TIME<<"Collision number= "<<a<<endl;
		cout<<"Collider x= "<<C.TC.S.P[C.TC.collider].coordinate.x<<"\t y =";
		cout<<C.TC.S.P[C.TC.collider].coordinate.y<<"\t z = ";
		cout<<C.TC.S.P[C.TC.collider].coordinate.z<<endl;
		cout<<"Partner x = "<<C.TC.S.P[C.TC.Partner[C.TC.collider]].coordinate.x<<"\t y ="; 
		cout<<C.TC.S.P[C.TC.Partner[C.TC.collider]].coordinate.y<<"\t z = ";
		cout<<C.TC.S.P[C.TC.Partner[C.TC.collider]].coordinate.z<<endl;
		cout<<"Event type = "<<C.TC.Event_type[C.TC.collider]<<"\t Collision type=";
		cout<<C.TC.Collision_type[C.TC.collider]<<endl;
		//Makes sure to start counting the virial after the burning time
		if(int(double(a)/double(C.TC.S.N)) == C.TC.S.burn_percent*C.TC.S.nsweep)	
			{C.virial = 0; burn_time = C.TC.S.TIME;}

		if(a == int(0.1*double(C.TC.S.N*C.TC.S.nsweep)))		
		{
			percent_time = time(NULL);
			cout<<"10 percent done"<<endl;
			out.open(FileName, ios::app);
                                out<<"No. of collisions = "<<a<<endl;
                                out<<"No. of ghost collisions: "<<ghost_counter<<endl;
                                out<<"Time in DMD units: "<<C.TC.S.TIME<<endl;
                                out<<"No. of bond_cols: "<<C.bondcol_counter<<"\t No. of sqw_cols: ";
				out<<C.swcol_counter<<"\t No. of cellchanges: "<<C.cellchange_counter;
				out<<"\t No. of nbrlistchanges: "<<C.nbrmove_counter;
				out<<"\t No. of FP updates"<<fpupdate_counter<<endl;
                                out<<"Start time: "<<ctime(&initialization_time)<<"Completion time: ";
				out<<ctime(&percent_time)<<endl;
                        out.close();

		}
		if(a == int(0.2*double(C.TC.S.N*C.TC.S.nsweep)))
		{
			percent_time = time(NULL);
			cout<<"20 percent done"<<endl;
			out.open(FileName, ios::app);
                                out<<"No. of collisions = "<<a<<endl;
                                out<<"No. of ghost collisions: "<<ghost_counter<<endl;
                                out<<"Time in DMD units: "<<C.TC.S.TIME<<endl;
                                out<<"No. of bond_cols: "<<C.bondcol_counter<<"\t No. of sqw_cols: ";
				out<<C.swcol_counter<<"\t No. of cellchanges: "<<C.cellchange_counter;
				out<<"\t No. of nbrlistchanges: "<<C.nbrmove_counter;
				out<<"\t No. of FP updates"<<fpupdate_counter<<endl;
                                out<<"Start time: "<<ctime(&initialization_time)<<"Completion time: ";
				out<<ctime(&percent_time)<<endl;
                        out.close();
		}
		if(a == int(0.3*double(C.TC.S.N*C.TC.S.nsweep)))
		{
			percent_time = time(NULL);
			cout<<"30 percent done"<<endl;
			out.open(FileName, ios::app);
                                out<<"No. of collisions = "<<a<<endl;
                                out<<"No. of ghost collisions: "<<ghost_counter<<endl;
                                out<<"Time in DMD units: "<<C.TC.S.TIME<<endl;
                                out<<"No. of bond_cols: "<<C.bondcol_counter<<"\t No. of sqw_cols: ";
				out<<C.swcol_counter<<"\t No. of cellchanges: "<<C.cellchange_counter;
				out<<"\t No. of nbrlistchanges: "<<C.nbrmove_counter;
				out<<"\t No. of FP updates"<<fpupdate_counter<<endl;
                                out<<"Start time: "<<ctime(&initialization_time)<<"Completion time: ";
				out<<ctime(&percent_time)<<endl;
                        out.close();
		}
		if(a == int(0.4*double(C.TC.S.N*C.TC.S.nsweep)))
		{
			percent_time = time(NULL);
			cout<<"40 percent done"<<endl;
			out.open(FileName, ios::app);
                                out<<"No. of collisions = "<<a<<endl;
                                out<<"No. of ghost collisions: "<<ghost_counter<<endl;
                                out<<"Time in DMD units: "<<C.TC.S.TIME<<endl;
                                out<<"No. of bond_cols: "<<C.bondcol_counter<<"\t No. of sqw_cols: ";
				out<<C.swcol_counter<<"\t No. of cellchanges: "<<C.cellchange_counter;
				out<<"\t No. of nbrlistchanges: "<<C.nbrmove_counter;
				out<<"\t No. of FP updates"<<fpupdate_counter<<endl;
                                out<<"Start time: "<<ctime(&initialization_time)<<"Completion time: ";
				out<<ctime(&percent_time)<<endl;
                        out.close();
		}
		if(a == int(0.5*double(C.TC.S.N*C.TC.S.nsweep)))
		{
			percent_time = time(NULL);
			cout<<"50 percent done"<<endl;
			out.open(FileName, ios::app);
                                out<<"No. of collisions = "<<a<<endl;
                                out<<"No. of ghost collisions: "<<ghost_counter<<endl;
                                out<<"Time in DMD units: "<<C.TC.S.TIME<<endl;
                                out<<"No. of bond_cols: "<<C.bondcol_counter<<"\t No. of sqw_cols: ";
				out<<C.swcol_counter<<"\t No. of cellchanges: "<<C.cellchange_counter;
				out<<"\t No. of nbrlistchanges: "<<C.nbrmove_counter;
				out<<"\t No. of FP updates"<<fpupdate_counter<<endl;
                                out<<"Start time: "<<ctime(&initialization_time)<<"Completion time: ";
				out<<ctime(&percent_time)<<endl;
                        out.close();
		}
		if(a == int(0.6*double(C.TC.S.N*C.TC.S.nsweep)))
		{
			percent_time = time(NULL);
			cout<<"60 percent done"<<endl;
			out.open(FileName, ios::app);
                                out<<"No. of collisions = "<<a<<endl;
                                out<<"No. of ghost collisions: "<<ghost_counter<<endl;
                                out<<"Time in DMD units: "<<C.TC.S.TIME<<endl;
                                out<<"No. of bond_cols: "<<C.bondcol_counter<<"\t No. of sqw_cols: ";
				out<<C.swcol_counter<<"\t No. of cellchanges: "<<C.cellchange_counter;
				out<<"\t No. of nbrlistchanges: "<<C.nbrmove_counter;
				out<<"\t No. of FP updates"<<fpupdate_counter<<endl;
                                out<<"Start time: "<<ctime(&initialization_time)<<"Completion time: ";
				out<<ctime(&percent_time)<<endl;
                        out.close();
		}
		if(a == int(0.7*double(C.TC.S.N*C.TC.S.nsweep)))
		{
			percent_time = time(NULL);
			cout<<"70 percent done"<<endl;
			out.open(FileName, ios::app);
                                out<<"No. of collisions = "<<a<<endl;
                                out<<"No. of ghost collisions: "<<ghost_counter<<endl;
                                out<<"Time in DMD units: "<<C.TC.S.TIME<<endl;
                                out<<"No. of bond_cols: "<<C.bondcol_counter<<"\t No. of sqw_cols: ";
				out<<C.swcol_counter<<"\t No. of cellchanges: "<<C.cellchange_counter;
				out<<"\t No. of nbrlistchanges: "<<C.nbrmove_counter;
				out<<"\t No. of FP updates"<<fpupdate_counter<<endl;
                                out<<"Start time: "<<ctime(&initialization_time)<<"Completion time: ";
				out<<ctime(&percent_time)<<endl;
                        out.close();
		}
		if(a == int(0.8*double(C.TC.S.N*C.TC.S.nsweep)))
		{
			percent_time = time(NULL);
			cout<<"80 percent done"<<endl;
			out.open(FileName, ios::app);
                                out<<"No. of collisions = "<<a<<endl;
                                out<<"No. of ghost collisions: "<<ghost_counter<<endl;
                                out<<"Time in DMD units: "<<C.TC.S.TIME<<endl;
                                out<<"No. of bond_cols: "<<C.bondcol_counter<<"\t No. of sqw_cols: ";
				out<<C.swcol_counter<<"\t No. of cellchanges: "<<C.cellchange_counter;
				out<<"\t No. of nbrlistchanges: "<<C.nbrmove_counter;
				out<<"\t No. of FP updates"<<fpupdate_counter<<endl;
                                out<<"Start time: "<<ctime(&initialization_time)<<"Completion time: ";
				out<<ctime(&percent_time)<<endl;
                        out.close();
		}
		if(a == int(0.9*double(C.TC.S.N*C.TC.S.nsweep)))
		{
			percent_time = time(NULL);
			cout<<"90 percent done"<<endl;
			out.open(FileName, ios::app);
                                out<<"No. of collisions = "<<a<<endl;
                                out<<"No. of ghost collisions: "<<ghost_counter<<endl;
                                out<<"Time in DMD units: "<<C.TC.S.TIME<<endl;
                                out<<"No. of bond_cols: "<<C.bondcol_counter<<"\t No. of sqw_cols: ";
				out<<C.swcol_counter<<"\t No. of cellchanges: "<<C.cellchange_counter;
				out<<"\t No. of nbrlistchanges: "<<C.nbrmove_counter;
				out<<"\t No. of FP updates"<<fpupdate_counter<<endl;
                                out<<"Start time: "<<ctime(&initialization_time)<<"Completion time: ";
				out<<ctime(&percent_time)<<endl;
                        out.close();
		}
	}
 	exit:
//Writing the final dump and DCD files
//	C.AllParticlePositionUpdater(C.TC.S.N, C.TC.S.P, C.TC.S.TIME, C.TC.S.fpupdate_TIME);
	writer.WriteDump(C.TC.S.TIME, C.TC.S.P, C.TC.S.N, C.TC.S.L, C.TC.S.filename_out);
	writer.DCDWriteStep(C.TC.S.P, C.TC.S.N, C.TC.S.filename_out);
	writer.DCDHeader(start_time, C.TC.S.N, C.TC.S.filename_out);
	writer.DCDFile(C.TC.S.filename_out);

	cout<<"Main simulation loop finished"<<endl;
//Final KE and TE
	enf_ke = C.TC.S.ke(C.TC.S.P, C.TC.S.N);	
	enf_te = enf_ke + C.TC.S.potential_energy;
	momentum = C.TC.S.net_momentum(C.TC.S.P, C.TC.S.N);
//Time variation file printing
	writer.TimeVarFile(C.TC.S.TIME,enf_ke,C.TC.S.potential_energy,enf_te,momentum,C.TC.S.filename_out);
//Calculation of compressibility factor, rate per col and final vol_frac
	comp_fac = C.virial/(3*double(C.TC.S.N)*(C.TC.S.TIME-burn_time)*C.TC.S.temperature);
	rate = a/C.TC.S.TIME;
	for(int i = 0; i< C.TC.S.N; i++)
	{
		vol_frac = vol_frac + PI*(C.TC.S.P[i].diameter/double(C.TC.S.L))*(C.TC.S.P[i].diameter/double(C.TC.S.L))*(C.TC.S.P[i].diameter/double(C.TC.S.L))/6;
	}
	cout<<"Writing file for final data collected from the simultions, txt format"<<endl;
	cout<<"Number of bond collisions:"<<C.bondcol_counter<<"\tNumber of square well collisions: ";
	cout<<C.swcol_counter<<endl;
  
	time_t completion_time = time(NULL);

	out.open(FileName, ios::app);
		out<<"Output file for L=\t"<<C.TC.S.L<<"\t & N=\t"<<C.TC.S.N<<endl;
		out<<"No. of collisions = "<<a<<endl;
		out<<"Initial TE= "<<eni_te<<"\t Per particle basis "<<eni_te/double(C.TC.S.N)<<"Initial KE= ";
		out<<eni_ke<<"\t Per particle basis "<<eni_ke/double(C.TC.S.N)<<endl;
		out<<"Final TE = "<<enf_te<<"\t Per particle basis "<<enf_te/double(C.TC.S.N)<<"Final KE= ";
		out<<enf_ke<<"\t Per particle basis "<<enf_ke/double(C.TC.S.N)<<endl;
		out<<"Compressibility factor = "<<comp_fac<<"Volume Fraction = "<<vol_frac<<endl;
		out<<"No. of ghost collisions: "<<ghost_counter<<endl;
		out<<"Time in system units: "<<C.TC.S.TIME<<"\t Burn time"<<burn_time<<endl;
		out<<"No. of bond_cols: "<<C.bondcol_counter<<"\t No. of sqw_cols: "<<C.swcol_counter;
		out<<"\t Cellchanges:"<<C.cellchange_counter<<"\t Nbrlistchanges:"<<C.nbrmove_counter<<endl;
		out<<"Start time: "<<ctime(&initialization_time)<<"Completion time: ";
		out<<ctime(&completion_time)<<endl;
	out.close();
//Finally writing finished RDF
	rdf.RDF_write(C.TC.S.L, C.TC.S.N, C.TC.S.maxbin, C.TC.S.filename_out); 
}
