//Here we define functions for timecalc class	
#include "timecalc.h" 
//Minimum image convention to find the nearest neighbor

//Initialising the time array and to store the successive collisions
void TimeCalc::timearray_initialization(int N)
{
	for(int i = 0; i < N; i++) 
	{
		TimeList.push_back(timbig);
		Partner.push_back(N);
		Collision_type.push_back(0);
		Event_type.push_back(0);
	}
	for(int i =0; i<5*N; i++)
	{
		ColPars.push_back(0.0);
	}
}
//Time until collision to select from a vector of L's and epsilons
void TimeCalc::CollisionTime_ij(int i, int j, int event_type, vector<double> r_inner, vector<double> r_outer, vector<double> epsilon)
{
	double r;
	bool list_updated = false;
	Particle counter;
	counter.coordinate = S.P[i].coordinate - S.P[j].coordinate;
	counter.velocity = S.P[i].velocity - S.P[j].velocity;
//Updating to the actual positions to calculate collision times
	counter.coordinate=S.OneParticlePositionUpdater(counter, S.TIME, S.fpupdate_TIME);
	S.min_img(counter.coordinate, S.L);

	r = counter.coordinate.norm();

	if(event_type != 2)
	{
		cout<<"Wrong event_type assigned"<<endl;
		exit(1);
	}
	for(int k = 0; k<=r_inner.size()-1; k++)
	{
//This takes care of the multiple well cases: 1.0e08 represents hard wall, 0 represents being outside the well
		if(r < r_outer[k] && r > r_inner[k])// Particle is in that well
		{
			list_updated = true;
			if(k == 0 && r_inner.size()==1)	
				{CollisionTime_ij(i, j, 2, r_inner[k], r_outer[k], epsilon[k], 1.0e08, 0);}
			else if(k==0 && r_inner.size() > 1)
				{CollisionTime_ij(i, j, 2, r_inner[k], r_outer[k], epsilon[k], 1.0e08, epsilon[k+1]);}
			else if(k!=0 && k==r_inner.size()-1)
				{CollisionTime_ij(i, j, 2, r_inner[k], r_outer[k], epsilon[k], epsilon[k-1], 0);}
			else
				{CollisionTime_ij(i, j, 2, r_inner[k], r_outer[k], epsilon[k], epsilon[k-1], epsilon[k+1]);}
			break;
		}
		else
			{list_updated = false;}
	}
//This means list has not been updated and particle is outside all the wells
	if(!list_updated)	
	{
//This means that the particle is outside the potential range (so the current epsilon should be zero)
		if(r > r_outer[r_outer.size()-1])	
		{
//zero epsilon shows that the current epsilon is zero and the epsilon to the right is also zero
//Using S.L as the outer limit is a counter that tells us that the particle is outside the well
			CollisionTime_ij(i, j, 2, r_outer[r_outer.size()-1], S.L, 0.0, epsilon[r_outer.size()-1], 0.0);
		}
	}	
}
//General time until collision occurs// This is valid for negative epsilon  
//If epsilon_inner is 1e08, then it means that we are looking at a hard sphere collision
void TimeCalc::CollisionTime_ij(int i, int j, int event_type, double r_inner, double r_outer, double epsilon_now, double epsilon_inner, double epsilon_outer)
{
	double b=0, r2=0, t=0, D1=0, D2=0, c1=0, c2=0;
	Particle counter;
	int col_type =0;
//POTENTIAL < 0 MEANS SHOULDERS ARE INVOLVED, ADD CODE FOR THIS! 
	double potential_inner, potential_outer, r_inner2, r_outer2, reduced_mass;			

	counter.coordinate = S.P[i].coordinate - S.P[j].coordinate;
	counter.velocity = S.P[i].velocity - S.P[j].velocity;
//Updating to the actual positions to calculate collision times
	counter.coordinate=S.OneParticlePositionUpdater(counter, S.TIME, S.fpupdate_TIME);
	S.min_img(counter.coordinate, S.L);

	b = counter.bcalc();
	counter.velocity2 = counter.velocity.norm2();
	r2 = counter.coordinate.norm2();
	reduced_mass = S.P[i].mass*S.P[j].mass/(S.P[i].mass + S.P[j].mass);

	r_inner2 = r_inner*r_inner;
	r_outer2 = r_outer*r_outer;

//If potential_inner and potential_outer are both positive and any one of them is not equal to 1.0e8
// then it is in a square well currently
//For this kind of collision, the concerned distances are average of radii and sigma2
	if(event_type == 2)							
	{
		potential_inner = -epsilon_now + epsilon_inner;
		potential_outer = -epsilon_now + epsilon_outer;
	}
	if(event_type == 1)
	{
		//because bond potential is a square well with zero depth
		potential_inner = 1.0e08;							
		potential_outer = 1.0e08;
	}
//	cout<<"Coltype eval for i="<<i<<"\t j="<<j<<"\t event_type="<<event_type<<"\t  r_inner="<<r_inner<<"\t r_outer="<<r_outer<<"\t epsilon_now"<<epsilon_now<<"\t epsilon_inner"<<epsilon_inner<<"\t epsilon_outer="<<epsilon_outer<<endl;
//Checking if inside square well or bond well, depending on the type of collision
	c1 = r2 - r_inner2;
	c2 = r2 - r_outer2;
//Corresponding discriminants
	D1 = b*b - counter.velocity2*c1;					
	D2 = b*b - counter.velocity2*c2;
//Centers coming towards one another
	if(b < 0)								
	{
//particles inside square well (BY THE CURRENT FORMULATION, THIS IS ALWAYS TRUE, BUT STILL KEEP FOR SOME TIME)
//Taking the small value because it is the margin of error
		if(c2/r_outer2 <= 1.0e-4 && r_outer != S.L)						
		{
			if(D1 > 0)						//Cores Collision
			{
				t = (-b - sqrt(D1))/counter.velocity2;
				if(t < 0)
				{
					t = (-b + sqrt(D1))/counter.velocity2;
				}
				col_type = 1;
//				col_type = 1;
//Only bounce on left side possible //The left side has a hard wall
				if(epsilon_inner == 1.0e08)
					{col_type=1;}
//Particle is in a sq well right now, no hard wall on left side
				else if(potential_inner > 0.0)
				{
//Energy enough to get out of the well, or collision with outer wall in case of bond event
					if(b*b > 2*r_inner2*potential_inner/reduced_mass)
						{col_type = 2;}				//Transfer on left side
					else						//Energy not enough
						{col_type = 1;}				//Bounce on left side
				}
//Left side has a sq well, bounce never possible
				else if(potential_inner < 0.0)			
					{col_type=2;}				//Transfer on left side
			}
//Cores will not collide, attractive collision will occur, D1<0//no concept of attractive collision for bonded particles
			else 							
			{
				t = (-b + sqrt(D2))/counter.velocity2;
				col_type = 3;
//This means that it has an infinite wall at the outside
				if(epsilon_outer == 1.0e08)			
				{
					col_type = 3;
				}
//Particle in square well and not outside wells
				else if(potential_outer > 0.0)		
				{					
//Energy enough to get out of the well, or collision with outer wall in case of bond event. 
					if(b*b > 2*r_outer2*potential_outer/reduced_mass)
						//Transfer on the right side
						{col_type = 4;}				
					else						//Energy not enough
						{col_type = 3;}				//Bounce on right side
				}
				//Bounce can never happen on right side now
				else if(potential_outer < 0.0)			
					{col_type=4;}				//Transfer on right side
			}
		}
//This means that the particles are outside square well, only possible collisions are 1, 2, 0 (which is 4 for no boundaries)
		else if(c2/r_outer2 <= 1.0e-4 && r_outer == S.L)	
		{
			if(D1 > 0)						//Cores Collision
			{
				t = (-b - sqrt(D1))/counter.velocity2;
				if(t < 0)
				{
					t = (-b + sqrt(D1))/counter.velocity2;
				}
				col_type = 2;					//Making transfer the default
//Only bounce on left side possible //The left side has a hard wall
				if(epsilon_inner == 1.0e08)			 
					{col_type=1;}				
//Particle is in a sq well right now, no hard wall on left side
				else if(potential_inner > 0.0)
				{
//Energy enough to get out of the well, or collision with outer wall in case of bond event
					if(b*b > 2*r_inner2*potential_inner/reduced_mass)
						{col_type = 2;}				//Transfer on left side
					else						//Energy not enough
						{col_type = 1;}				//Bounce on left side
				}
//Left side has a sq well, bounce never possible
				else if(potential_inner < 0.0)			
					{col_type=2;}				//Transfer on left side
			}
//Particles will not collide at all, will just be going away from outside the well boundary
			else 							
			{
				t = timbig;
				col_type = 0;
			}
		}
		else
		{//This shouldn't happen anymore (IDEALLY), so let's see
		//(int i, int j, int event_type, double r_inner, double r_outer, double epsilon_now, double epsilon_inner, double epsilon_outer)
			col_type = 5;
			cout<<"Exiting because particles should not be outside wells anymore, coltype = ";
			cout<<col_type<<endl;
			cout<<"Particle1="<<i<<"\t Particle2="<<j<<"\t event_type="<<event_type<<"\t r_outer=";
			cout<<r_outer<<"\t c2="<<c2<<"\t r_outer2="<<r_outer2<<"\t S.L="<<S.L<<"\t r2="<<r2;
			cout<<"\t c2/r_outer2="<<c2/r_outer2<<endl;
			cout<<"P1x="<<S.P[i].coordinate.x<<"\tP1y="<<S.P[i].coordinate.y<<"\tP1z=";
			cout<<S.P[i].coordinate.z<<"\tP1vx="<<S.P[i].velocity.vx<<"\tP1vy=";
			cout<<S.P[i].velocity.vy<<"\tP1vz="<<S.P[i].velocity.vz<<endl;
			cout<<"P2x="<<S.P[j].coordinate.x<<"\tP2y="<<S.P[j].coordinate.y<<"\tP2z=";
			cout<<S.P[j].coordinate.z<<"\tP2vx="<<S.P[j].velocity.vx<<"\tP2vy=";
			cout<<S.P[j].velocity.vy<<"\tP2vz="<<S.P[j].velocity.vz<<endl;
			cout<<"TIME="<<S.TIME<<"\tfp_TIME="<<S.fpupdate_TIME<<endl;

			exit(1);
//			if(D2 > 0) 						//It'll be a capture
//			{
//				t = (-b - sqrt(D2))/counter.velocity2;
//				col_type = 2;
//			}
//			else 							//No collision
//			{
//				t = timbig;
//				col_type = 0;
//			}
		}
	}
//Centers going away from each other, b>0
	else									
	{
		if(c2/r_outer2 <= 1.0e-4 && r_outer != S.L)				//Inside square well
		{
			t = (-b + sqrt(D2))/counter.velocity2;
			col_type = 3;
//Means that this is a bond collision
			if(epsilon_outer ==  1.0e08)			
			{
				col_type = 3;
			}
//Particle in square well and not outside wells
			else if(potential_outer > 0.0)			
			{					
//Energy enough to get out of the well, or collision with outer wall in case of bond event
				if(b*b > 2*r_outer2*potential_outer/reduced_mass)	
					{col_type = 4;}				//Transfer on the right side
				else						//Energy not enough
					{col_type = 3;}				//Bounce on right side
			}
			else if(potential_outer < 0.0)			//Bounce can never happen on left side now
				{col_type=4;}				//Transfer on right side
		}
//Centers are going away and it is outside the square well, so can only have collision type 0
		else if(c2/r_outer2 <= 1.0e-4 && r_outer == S.L)
		{
			t=timbig;
			col_type=0;
		}
//Means centers going away from each other and it is outside the range, ideally shouldn't happen
		else // no collision 
		{
			col_type = 6;
			cout<<"Exiting,particles should not be outside wells anymore,coltype= "<<col_type<<endl;
			cout<<"Particle1="<<i<<"\t Particle2="<<j<<"\t event_type="<<event_type<<"\t r_outer=";
			cout<<r_outer<<"\t c2="<<c2<<"\t r_outer2="<<r_outer2<<"\t S.L="<<S.L<<"\t r2="<<r2;
			cout<<"\t c2/r_outer2="<<c2/r_outer2<<endl;
			cout<<"P1x="<<S.P[i].coordinate.x<<"\tP1y="<<S.P[i].coordinate.y<<"\tP1z=";
			cout<<S.P[i].coordinate.z<<"\tP1vx="<<S.P[i].velocity.vx<<"\tP1vy=";
			cout<<S.P[i].velocity.vy<<"\tP1vz="<<S.P[i].velocity.vz<<endl;
			cout<<"P2x="<<S.P[j].coordinate.x<<"\tP2y="<<S.P[j].coordinate.y<<"\tP2z=";
			cout<<S.P[j].coordinate.z<<"\tP2vx="<<S.P[j].velocity.vx<<"\tP2vy=";
			cout<<S.P[j].velocity.vy<<"\tP2vz="<<S.P[j].velocity.vz<<endl;
			cout<<"TIME="<<S.TIME<<"\tfp_TIME="<<S.fpupdate_TIME<<endl;
//			exit(1);
		}
	}
	
//Will this actually be needed? Because if we believe the accuracy, then it is already taken care of in the skipping wells portion//Wait to see if this actually happens
/*	if(t <= 1.0e-11 || (event_type == 1 && (col_type == 2 || col_type == 4)))
	{
		if(event_type == 2)					//Repitition of a square well collision
		{	//Transfer on right
			if(col_type == 4)			//4 -----> 2 or 0, checking criteria for 2
			{
				if(b<0)
				{
					col_type = 2;
					t = (-b - sqrt(D2))/counter.velocity2;
				}
				if(b>0)
				{ 
					if(r_outer == S.L)
					{
						col_type = 0;
						t = timbig;
					}
				}
			}
			else if(col_type == 2)			//2 -----> 1 or 3 or 4, checking criteria for them
			{
				if(D1>=0)
				{
					col_type = 1;
					t = (-b - sqrt(D1))/counter.velocity2;
					if(t < 0)
					{ t = (-b + sqrt(D1))/counter.velocity2;}
				}
				else				
				{
					t = (-b+sqrt(D2))/counter.velocity2;
					if(b*b > 2*r_outer2*potential/reduced_mass)	
						{col_type = 3;}
					else				
						{col_type = 4;}
				}
			}//bounce on left
			else if(col_type == 1)			//1 -----> 3 or 4, checking criteria for them
			{
				t = (-b+sqrt(D2))/counter.velocity2;
				if(b*b > 2*r_outer2*potential/reduced_mass)	
					{col_type = 3;}
				else				
					{col_type = 4;}
			}//Bounce on right
			else if(col_type == 4)
			{
				if(b<=0)
				{
					col_type = 1;
					t = (-b - sqrt(D1))/counter.velocity2;
					if(t < 0)
					{ t = (-b + sqrt(D1))/counter.velocity2;}
				}
				else				
				{
					t = (-b+sqrt(D2))/counter.velocity2;
					if(b*b > 2*r_outer2*potential/reduced_mass)	
						{col_type = 3;}
					else				
						{col_type = 4;}
				}
			}
		}
//NOT SURE ABOUT THIS, CHECK WHILE EXECUTING
		else if(event_type == 1)
		{
			if(col_type == 4 || col_type == 2)		//4 or 2 -----> 1, checking criteria
			{
				col_type = 1;
				t = (-b - sqrt(D1))/counter.velocity2;
				if(t < 0)
					{t = (-b + sqrt(D1))/counter.velocity2;}
			}
			else if(col_type == 3 || col_type == 1)		//3 or 1 -----> 3 , checking criteria
			{
				t = (-b+sqrt(D2))/counter.velocity2;
				col_type = 3;
			}	
		}
	}*/
	if(t < TimeList[i] && t > 1.0e-11)
	{
//		cout<<"i="<<i<<"\t j="<<j<<"\t t="<<t<<"\t r2="<<r2<<"\t c2="<<c2;
//		cout<<"\tevent_type="<<event_type<<"\tcol_type"<<col_type<<endl;
//		cout<<"\t P1x="<<S.P[i].coordinate.x<<"\t P1y="<<S.P[i].coordinate.y;
//		cout<<"\t P1z="<<S.P[i].coordinate.z<<endl;
//		cout<<"\t P2x="<<S.P[j].coordinate.x<<"\t P2y="<<S.P[j].coordinate.y;
//		cout<<"\t P2z="<<S.P[j].coordinate.z<<endl;
//		cout<<"\t P1vx="<<S.P[i].velocity.vx<<"\t P1vy="<<S.P[i].velocity.vy;
//		cout<<"\t P1vz="<<S.P[i].velocity.vz<<endl;
//		cout<<"\t P2vx="<<S.P[j].velocity.vx<<"\t P2vy="<<S.P[j].velocity.vy;
//		cout<<"\t P2vz="<<S.P[j].velocity.vz<<endl;
		TimeList[i] = t;
		Partner[i] = j;
		Collision_type[i] = col_type;
		Event_type[i] = event_type;
		ColPars[5*i] = r_inner;
		ColPars[5*i + 1] = epsilon_now;
		ColPars[5*i + 2] = r_outer;
		ColPars[5*i + 3] = epsilon_inner;
		ColPars[5*i + 4] = epsilon_outer;
	}
}
//For updating the time array
void TimeCalc::UpdateUpList(int a, int N)
{
	double l1, l2, epsilon;
	int event_type;
	bool particle_bonded = false, particle_nonbonded = false;
//no need to edit, no uplist molecules present, a=N means its an andersen thermostat
	if(a == N-1)	
		{return;}
	TimeList[a] = timbig;
	
	for(int j=a+1; j<N; j++)
	{
		for(int k=0; k < S.bondlist.size(); k++)
		{
			//Means the two particles are bonded ( cov or pseudo, doesn't matter)
			if(S.bondlist[k].partner1 == a && S.bondlist[k].partner2 == j) 
			{
				particle_bonded = true;
				l1 = (1-S.bond_delta)*S.bondlist[k].bondlength;
				l2 = (1+S.bond_delta)*S.bondlist[k].bondlength;
				epsilon = 0;
				event_type = 1;
				CollisionTime_ij(a, j, event_type, l1, l2, 0.0, 1.0e08, 1.0e08);
				break;
			}
			else 
				{particle_bonded = false;}	
		}
		//The two particles are not bonded, looked through the whole bondlist
		if(!particle_bonded)			
		{
			for(int k = 0; k<S.nonbondlist.size(); k++)
			{//Find the two particles in nonbondlist
				if((S.nonbondlist[k].partner1 == S.P[a].atomtype && S.nonbondlist[k].partner2 == S.P[j].atomtype) || (S.nonbondlist[k].partner1 == S.P[j].atomtype && S.nonbondlist[k].partner2 == S.P[a].atomtype))	
				{
					particle_nonbonded = true;
					event_type = 2;
					CollisionTime_ij(a, j, event_type, S.nonbondlist[k].L1, S.nonbondlist[k].L2, S.nonbondlist[k].epsilon);
					break; 
				}
				else
					{particle_nonbonded = false;}
			}
			if(!particle_nonbonded)		//Particle pairing not found in nonbondedlist
			{
				cout<<"Error, particle pairing not found in nonbondlist and bondlist, weird weird"<<endl;
				exit(1);
			}
		}
	}
}
//Checking if "a" was the collision partner of any other particle in the list
void TimeCalc::UpdateDownList(int a, int N)
{
	double l1, l2, epsilon;
	int event_type;
	bool particle_bonded = false, particle_nonbonded = false;

	if(a == 0)					//Because no particle below it for updating downlist
		{return;}

	for(int j=0; j<a; j++)
	{
		for(int k=0; k < S.bondlist.size(); k++)
		{//Means the two particles are bonded
			if(S.bondlist[k].partner1 == j && S.bondlist[k].partner2 == a) 
			{
				particle_bonded = true;
				l1 = (1-S.bond_delta)*S.bondlist[k].bondlength;
				l2 = (1+S.bond_delta)*S.bondlist[k].bondlength;
				epsilon = 0;
				event_type = 1;
				CollisionTime_ij(j, a, event_type, l1, l2, 0.0, 1.0e08, 1.0e08);
				break;
			}
			else
				{particle_bonded = false;}	
		}//The two particles are not bonded, looked through the whole bondlist
		if(!particle_bonded)			
		{
			for(int k = 0; k<S.nonbondlist.size(); k++)
			{//Find the two particles in nonbondlist
				if((S.nonbondlist[k].partner1 == S.P[a].atomtype && S.nonbondlist[k].partner2 == S.P[j].atomtype) || (S.nonbondlist[k].partner1 == S.P[j].atomtype && S.nonbondlist[k].partner2 == S.P[a].atomtype))
				{
					particle_nonbonded = true;
					event_type = 2;
					CollisionTime_ij(j, a, event_type, S.nonbondlist[k].L1, S.nonbondlist[k].L2, S.nonbondlist[k].epsilon);
					break; 
				}
				else
					{particle_nonbonded = false;}
			}
			if(!particle_nonbonded)		//Particle pairing not found in nonbondedlist
			{
				cout<<"Error, particle pairing not found in nonbondlist and bondlist, weird stuff bro"<<endl;
				exit(1);
			}
		}
	}
}
//Initialise the time array, also for reseting it after some particle collisions
void TimeCalc::Collision_time(double L, int N)
{
	int i;
	minimum_time = timbig;
//Reinitialising the partner and the time array everytime this function is called, is included in the updateuplist part
	for(i=0; i<N; i++)
	{
		UpdateUpList(i, N);
		if(TimeList[i]<minimum_time)
		{
			minimum_time = TimeList[i];
			collider = i;
		}
	}
}
//Updating uplist in the case of cell list usage
void TimeCalc::UpdateUpCellList(int i, int N)
{
	double l1, l2, epsilon;
	int event_type;
	bool particle_bonded = false, particle_nonbonded = false;

	int j, cell_x, cell_y, cell_z, cell_xcalc, cell_ycalc, cell_zcalc;
	int m1, m2, m3, cell_ID;
	int x_l, x_u, y_l, y_u, z_l, z_u;
	double time_cellchange;
	TimeList[i] = timbig;
	//Extracting cell_ID
	cell_x = S.P[i].cell_ID%10;
	cell_y = int(double(S.P[i].cell_ID%100)/10);
	cell_z = int(double(S.P[i].cell_ID)/100);
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
					{cell_xcalc = S.ncell_max;}
				if(cell_ycalc < 0)
					{cell_ycalc = S.ncell_max;}
				if(cell_zcalc < 0)
					{cell_zcalc = S.ncell_max;}
				if(cell_xcalc > S.ncell_max)
					{cell_xcalc = 0;}
				if(cell_ycalc > S.ncell_max)
					{cell_ycalc = 0;}
				if(cell_zcalc > S.ncell_max)
					{cell_zcalc = 0;}
				
				cell_ID = 100*cell_zcalc + 10*cell_ycalc + cell_xcalc;
				j = S.HeadList[cell_ID];
				while(j >= 0)
				{
					if (j > i)
					{
						for(int k=0; k < S.bondlist.size(); k++)
						{//Means the two particles are bonded
							if(S.bondlist[k].partner1==i && S.bondlist[k].partner2 ==j)		
							{
								particle_bonded = true;
								l1 = (1-S.bond_delta)*S.bondlist[k].bondlength;
								l2 = (1+S.bond_delta)*S.bondlist[k].bondlength;
								epsilon = 0;
								event_type = 1;
								CollisionTime_ij(i, j, event_type, l1, l2, 0.0, 1.0e08, 1.0e08);
								break;
							}
							else 
								{particle_bonded = false;}	
						}//Particles are not bonded, looked through the whole bondlist
						if(!particle_bonded)
						{
							for(int k = 0; k<S.nonbondlist.size(); k++)
							{//Find the two particles in nonbondlist
								if((S.nonbondlist[k].partner1 == S.P[i].atomtype && S.nonbondlist[k].partner2 == S.P[j].atomtype) || (S.nonbondlist[k].partner1 == S.P[j].atomtype && S.nonbondlist[k].partner2 == S.P[i].atomtype)) 	
								{
									particle_nonbonded = true;
									event_type = 2;
									CollisionTime_ij(i, j, event_type, S.nonbondlist[k].L1, S.nonbondlist[k].L2, S.nonbondlist[k].epsilon);
									break; 
								}
								else
									{particle_nonbonded = false;}
							}//Particle pairing not found in nonbondedlist
							if(!particle_nonbonded)
							{
								cout<<"Error, particle pairing not found in nonbondlist and bondlist, weird weird"<<endl;
								exit(1);
							}
						}
						j = S.LinkList[j];
					}
					else
					{
						j = S.LinkList[j];
						continue;
					}
				}
			}
		}
	}
	time_cellchange = CellChangeTime(i, S.cellwidth, cell_x, cell_y, cell_z);
	if(time_cellchange < TimeList[i] && time_cellchange != 0)
	{
		TimeList[i] = time_cellchange;
		Event_type[i] = 3;
		Collision_type[i] = 0;
		Partner[i] = N; 
	}
}
//updating the list for all the particles which are below in numbering than this particular particle for cell list
void TimeCalc::UpdateDownCellList(int i, int N)
{
	double l1, l2, epsilon;
	int event_type;
	bool particle_bonded = false, particle_nonbonded = false;

	int j, cell_x, cell_y, cell_z, cell_xcalc, cell_ycalc, cell_zcalc;
	int m1, m2, m3, cell_ID;
	int x_l, x_u, y_l, y_u, z_l, z_u;
	double time_cellchange;

	//Extracting cell_ID
	cell_x = S.P[i].cell_ID%10;
	cell_y = int(double(S.P[i].cell_ID%100)/10);
	cell_z = int(double(S.P[i].cell_ID)/100);
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
					{cell_xcalc = S.ncell_max;}
				if(cell_ycalc < 0)
					{cell_ycalc = S.ncell_max;}
				if(cell_zcalc < 0)
					{cell_zcalc = S.ncell_max;}
				if(cell_xcalc > S.ncell_max)
					{cell_xcalc = 0;}
				if(cell_ycalc > S.ncell_max)
					{cell_ycalc = 0;}
				if(cell_zcalc > S.ncell_max)
					{cell_zcalc = 0;}
				
				cell_ID = 100*cell_zcalc + 10*cell_ycalc + cell_xcalc;
				j = S.HeadList[cell_ID];
				while(j >= 0)
				{
					if (j < i)
					{
						for(int k=0; k < S.bondlist.size(); k++)
						{//Means the two particles are bonded
							if(S.bondlist[k].partner1==j && S.bondlist[k].partner2 ==i)
							{
								particle_bonded = true;
								l1 = (1-S.bond_delta)*S.bondlist[k].bondlength;
								l2 = (1+S.bond_delta)*S.bondlist[k].bondlength;
								epsilon = 0;
								event_type = 1;
								CollisionTime_ij(j, i, event_type, l1, l2, 0.0, 1.0e08, 1.0e08);
								break;
							}
							else 
								{particle_bonded = false;}	
						}//Particles are not bonded, looked through the whole bondlist
						if(!particle_bonded)
						{
							for(int k = 0; k<S.nonbondlist.size(); k++)
							{//Find the two particles in nonbondlist
								if((S.nonbondlist[k].partner1 == S.P[i].atomtype && S.nonbondlist[k].partner2 == S.P[j].atomtype) || (S.nonbondlist[k].partner1 == S.P[j].atomtype && S.nonbondlist[k].partner2 == S.P[i].atomtype))
								{
									particle_nonbonded = true;
									event_type = 2;
									CollisionTime_ij(j, i, event_type, S.nonbondlist[k].L1, S.nonbondlist[k].L2, S.nonbondlist[k].epsilon);
									break; 
								}
								else
									{particle_nonbonded = false;}
							}//Particle pairing not found in nonbondedlist
							if(!particle_nonbonded)
							{
								cout<<"Error, particle pairing not found in nonbondlist, weird stuff bro"<<endl;
								exit(1);
							}
						}
						j = S.LinkList[j];
					}
					else
					{	
						j = S.LinkList[j];
						continue;
					}
				}
			}
		}
	}
}
double TimeCalc::CellChangeTime(int i, double cellsize, int cell_x, int cell_y, int cell_z)
{
	double time_x, time_y, time_z, time_min;
	XYZ coord;
	coord=S.OneParticlePositionUpdater(S.P[i], S.TIME, S.fpupdate_TIME);

	if(S.P[i].velocity.vx > 0)
		{time_x = (double(cell_x+1)*cellsize-coord.x)/S.P[i].velocity.vx;}
	//To take care of the sign
	else
		{time_x = (double(cell_x)*cellsize - coord.x)/S.P[i].velocity.vx;}	

	if(S.P[i].velocity.vy > 0)
		{time_y = (double(cell_y+1)*cellsize-coord.y)/S.P[i].velocity.vy;}
	//To take care of the sign
	else
		{time_y = (double(cell_y)*cellsize - coord.y)/S.P[i].velocity.vy;}	

	if(S.P[i].velocity.vz > 0)
		{time_z = (double(cell_z+1)*cellsize-coord.z)/S.P[i].velocity.vz;}
	//To take care of the sign
	else
		{time_z = (double(cell_z)*cellsize - coord.z)/S.P[i].velocity.vz;}

	//Sorting the shortest time for cell crossings
	if(time_x < time_y)
		{time_min = time_x;}
	else
		{time_min = time_y;}
	if(time_z < time_min)
	{	
		time_min = time_z;
	}

	if(time_min < 1.0e-11)
		{return timbig;}
	else
		{return time_min;}
}
//Initialize time list in case of using cell-list
void TimeCalc::Collision_time_celllist(double L, int N)
{
	int i;	

	minimum_time = timbig;
//Reinitialising the partner and the time array everytime this function is called, is included in the updateuplist part
	for(i=0; i<N; i++)
	{
		UpdateUpCellList(i, N);
		if(TimeList[i]<minimum_time)
		{
			minimum_time = TimeList[i];
			collider = i;
		}		
	}
}
//Updating timelist for j's > i using neighborlist
void TimeCalc::UpdateUpNeighborList(int i, int N)
{
	double l1, l2, epsilon;
	int event_type;
	bool particle_bonded = false, particle_nonbonded = false;

	int j, nbr_num;
	int cell_x, cell_y, cell_z;
	double time_cellchange;
	if(i >= N)
	{
		cout<<"Wrong assignment for particle time calculation in neighboruplist"<<endl;
		return;
	}

	TimeList[i] = timbig;
	//Extracting cell_ID
	cell_x = S.P[i].cell_ID%10;
	cell_y = int(double(S.P[i].cell_ID%100)/10);
	cell_z = int(double(S.P[i].cell_ID)/100);

	if(i != 0)
		{nbr_num = S.NeighborHeadList[i] - S.NeighborHeadList[i-1] - 1;}
	else
		{nbr_num = S.NeighborHeadList[i] + 1;}
	for(int count = 0; count < nbr_num; count++)
	{
		j = S.NeighborLinkList[S.NeighborHeadList[i] - count];
		if(j > i)
		{
			for(int k=0; k < S.bondlist.size(); k++)
			{//Means the two particles are bonded
				if(S.bondlist[k].partner1 == i && S.bondlist[k].partner2 == j) 
				{
					particle_bonded = true;
					l1 = (1-S.bond_delta)*S.bondlist[k].bondlength;
					l2 = (1+S.bond_delta)*S.bondlist[k].bondlength;
					epsilon = 0;
					event_type = 1;	
					CollisionTime_ij(i, j, event_type, l1, l2, 0.0, 1.0e08, 1.0e08);
					break;
				}	
				else 
					{particle_bonded = false;}	
			}//The two particles are not bonded, looked through the whole bondlist
			if(!particle_bonded)
			{
				for(int k = 0; k<S.nonbondlist.size(); k++)
				{//Find the two particles in nonbondlist
					if((S.nonbondlist[k].partner1 == S.P[i].atomtype && S.nonbondlist[k].partner2 == S.P[j].atomtype) || (S.nonbondlist[k].partner1 == S.P[j].atomtype && S.nonbondlist[k].partner2 == S.P[i].atomtype)) 
					{
						particle_nonbonded = true;
						event_type = 2;
						CollisionTime_ij(i, j, event_type, S.nonbondlist[k].L1, S.nonbondlist[k].L2, S.nonbondlist[k].epsilon);
						break; 
					}
					else
						{particle_nonbonded = false;}
				}
				if(!particle_nonbonded)		//Particle pairing not found in nonbondedlist
				{
					cout<<"Error, particle pairing not found in nonbondlist, weird weird"<<endl;
					exit(1);
				}
			}
		}	
	}
	time_cellchange = CellChangeTime(i, S.cellwidth, cell_x, cell_y, cell_z);
	if(time_cellchange < TimeList[i] && time_cellchange != 0)
	{
		TimeList[i] = time_cellchange;
		Event_type[i] = 3;
		Collision_type[i] = 0;
		Partner[i] = N; 
	}	
}
//Updating timelist for j's < i using neighborlist
void TimeCalc::UpdateDownNeighborList(int i, int N)
{
	double l1, l2, epsilon;
	int event_type;
	bool particle_bonded = false, particle_nonbonded = false;

	int j, nbr_num, chainno_i, chainno_j;
	if(i >= N)
	{
		cout<<"Wrong assignment for particle time calculation in neighbordownlist"<<endl;
		return;
	}

	if(i != 0)
		{nbr_num = S.NeighborHeadList[i] - S.NeighborHeadList[i-1] - 1;}
	else
		{nbr_num = S.NeighborHeadList[i] + 1;}
	for(int count = 0; count < nbr_num; count++)
	{
		j = S.NeighborLinkList[S.NeighborHeadList[i] - count];
		if(j < i)
		{
			for(int k=0; k < S.bondlist.size(); k++)
			{//Means the two particles are bonded
				if(S.bondlist[k].partner1 == j && S.bondlist[k].partner2 == i) 
				{
					particle_bonded = true;
					l1 = (1-S.bond_delta)*S.bondlist[k].bondlength;
					l2 = (1+S.bond_delta)*S.bondlist[k].bondlength;
					epsilon = 0;
					event_type = 1;	
					CollisionTime_ij(j, i, event_type, l1, l2, 0.0, 1.0e08, 1.0e08);
					break;
				}	
				else 
					{particle_bonded = false;}	
			}//The two particles are not bonded, looked through the whole bondlist
			if(!particle_bonded)
			{
				for(int k = 0; k<S.nonbondlist.size(); k++)
				{//Find the two particles in nonbondlist
					if((S.nonbondlist[k].partner1 == S.P[i].atomtype && S.nonbondlist[k].partner2 == S.P[j].atomtype) || (S.nonbondlist[k].partner1 == S.P[j].atomtype && S.nonbondlist[k].partner2 == S.P[i].atomtype)) 	
					{
						particle_nonbonded = true;
						event_type = 2;
						CollisionTime_ij(j, i, event_type, S.nonbondlist[k].L1, S.nonbondlist[k].L2, S.nonbondlist[k].epsilon);
						break; 
					}
					else
						{particle_nonbonded = false;}
				}
				if(!particle_nonbonded)		//Particle pairing not found in nonbondedlist
				{
					cout<<"Error, particle pairing not found in nonbondlist, weird stuff bro"<<endl;
					exit(1);
				}
			}
		}	
	}
}
//Calculating timelist for all using neighborlist
void TimeCalc::Collision_time_neighborlist(double L, int N)
{
	int i;	
	minimum_time = timbig;
//Reinitialising the partner and the time array everytime this function is called, is included in the updateuplist part
	for(i=0; i<N; i++)
	{
		UpdateUpNeighborList(i, N);
		if(TimeList[i]<minimum_time)
		{
			minimum_time = TimeList[i];
			collider = i;
		}		
	}
}
