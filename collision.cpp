//Here we define functions for collision class
#include "collision.h" 

double Collision::recalc_b(Particle &a, double L, int col_type, double r_inner2, double r_outer2)
{
	double b, r2, rij1 = 0, rij2 = 0;

	r2 = a.coordinate.norm2();
	rij1 = sqrt(r2/r_inner2);
	rij2 = sqrt(r2/r_outer2);
	//Means a core collision
	if(col_type == 1)										
	{
		if(1-rij1 >= 0)
		{
			b = a.bcalc();
			return b; 
		}
		else
		{
			TC.S.min_img(a.coordinate, L);
			b = a.bcalc();
			return b;
		}
	}
	else
	{
		if(1-rij2 >= 0)
		{
			b = a.bcalc();
			return b; 
		}
		else
		{
			TC.S.min_img(a.coordinate, L);
			b = a.bcalc();
			return b;
		}
	}
}

int Collision::AndersenThermostat()
{
	int n;

	double vx, vy, vz, v2 = 0, net_mass=0.0;
	VEL thermostat_vel, com_vel, deltav, momentum_counter, after_com_vel;

	Particle n_old, n_new;
	XYZ deltax;
//Will be used to obtain a seed for the random number engine
	std::random_device rd;  									
//Standard mersenne_twister_engine seeded with rd()
	std::mt19937 gen(rd());										
//	std::mt19937 gen{0};
//For selecting particle number
	std::uniform_real_distribution<> par(0.0, 1.0);							
	n = int(double(TC.S.N)*par(gen));
	if (n == TC.collider || n == TC.col_partner)
	{
		n = int(double(TC.S.N)*par(gen));
		if (n == TC.collider || n == TC.col_partner)
		{
			n = int(double(TC.S.N)*par(gen));
		}
	}
//A possible error can be that you aren't reforming timelist after changing the particle velocities, for all the particles, but it might not be necessary as low velocity change
//Moving the particle to the actual position
	TC.S.P[n].coordinate=TC.S.OneParticlePositionUpdater(TC.S.P[n],TC.S.TIME,TC.S.fpupdate_TIME);
	n_old=TC.S.P[n];
//Will be used to obtain a seed for the random number engine
//	std::random_device rd1;
//  	std::mt19937 gen1{rd1()};
//	std::mt19937 gen1{0};
//For getting a gaussian velocity profile
//	std::normal_distribution<> d1{0,1};
//	std::normal_distribution<> d2{0,1};
//	std::normal_distribution<> d3{0,1};	
//	std::uniform_real_distribution<> d1{-1,1};
//	std::uniform_real_distribution<> d2{-1,1};
//	std::uniform_real_distribution<> d3{-1,1};

//Reinitialized non-normalized velocity
	thermostat_vel.vx = TC.S.RandomGaussianNumber()*sqrt(TC.S.temperature/TC.S.P[n].mass);
	thermostat_vel.vy = TC.S.RandomGaussianNumber()*sqrt(TC.S.temperature/TC.S.P[n].mass);
	thermostat_vel.vz = TC.S.RandomGaussianNumber()*sqrt(TC.S.temperature/TC.S.P[n].mass);

	for (int i = 0; i<TC.S.N; i++)
	{
		if(i == n)
			{momentum_counter = thermostat_vel*TC.S.P[i].mass;}
		else
			{momentum_counter = TC.S.P[i].velocity*TC.S.P[i].mass;}
		com_vel=com_vel+momentum_counter;
		net_mass=net_mass+TC.S.P[i].mass;
	}

	com_vel = com_vel/net_mass;
	cout<<"Before thermostat, COMv.vx="<<com_vel.vx<<"\t COMv.vy="<<com_vel.vy<<"\t COMv.vz="<<com_vel.vz<<endl;

	TC.S.P[n].velocity = thermostat_vel;
	TC.S.P[n].velocity2 = TC.S.P[n].velocity.norm2();
	for (int i = 0; i<TC.S.N; i++)
	{
		TC.S.P[i].velocity = TC.S.P[i].velocity - com_vel;
		TC.S.P[i].velocity2 = TC.S.P[i].velocity.norm2();
		momentum_counter=TC.S.P[i].velocity*TC.S.P[i].mass;
		after_com_vel=after_com_vel+momentum_counter;
	}
	after_com_vel=after_com_vel/net_mass;
	cout<<"Using thermostat reformed velocity, COMv.vx="<<after_com_vel.vx<<"\t COMv.vy="<<after_com_vel.vy<<"\t COMv.vz="<<after_com_vel.vz<<endl;

//	TC.S.P[n].velocity2 = TC.S.P[n].velocity.norm2();
//Reinitialized the velocity, updating to make sure it is the true position.
	TC.S.P[n].coordinate=TC.S.OneParticlePositionBackwarder(TC.S.P[n],TC.S.TIME,TC.S.fpupdate_TIME);
//	n_new=TC.S.P[n];
//	n_new.coordinate=TC.S.OneParticlePositionUpdater(n_new,TC.S.TIME,TC.S.fpupdate_TIME);
//	deltax= n_new.coordinate-n_old.coordinate;
//	deltav.vx=deltax.x/(TC.S.TIME-TC.S.fpupdate_TIME);
//	deltav.vy=deltax.y/(TC.S.TIME-TC.S.fpupdate_TIME);
//	deltav.vz=deltax.z/(TC.S.TIME-TC.S.fpupdate_TIME);
//	TC.S.P[n].velocity=TC.S.P[n].velocity+deltav;
//	TC.S.P[n].velocity2 = TC.S.P[n].velocity.norm2();
//	cout<<"Thermostat particle= "<<n<<"\t deltax_x= "<<deltax.x<<"\t deltax_y= "<<deltax.y<<"\t deltax_z="<<deltax.z<<endl;
//	cout<<"deltav.vx= "<<deltav.vx<<"\t deltav.vy= "<<deltav.vy<<"\t deltav.vz= "<<deltav.vz<<endl;

//Moving the particle back in time but with the new velocity
	return n;
}

void Collision::Bump(double L, int N)
{
	int count = 0, vmax_i = 0, n;
	double r_inner2, r_outer2, b = 0, dpx = 0, dpy = 0, dpz = 0, r2 = 0, D = 0, reduced_mass, epsilon_now, epsilon_inner, epsilon_outer, potential_inner, potential_outer;
	Particle calculator;
	for(count=0; count<N; count++)
	{
		TC.TimeList[count] = TC.TimeList[count] - TC.minimum_time;	
	}
	TC.S.TIME = TC.S.TIME + TC.minimum_time;
//Update cell list when you change cells
	if(TC.Event_type[TC.collider] == 3)
	{
		TC.S.CellList();
		cellchange_counter = cellchange_counter + 1;
		didcellchange = true;
	}
	else
	{
//Moving just the colliding particles ahead to their true positions
		TC.col_partner = TC.Partner[TC.collider];	

		TC.S.P[TC.collider].coordinate=TC.S.OneParticlePositionUpdater(TC.S.P[TC.collider], TC.S.TIME, TC.S.fpupdate_TIME);
		TC.S.P[TC.col_partner].coordinate=TC.S.OneParticlePositionUpdater(TC.S.P[TC.col_partner], TC.S.TIME, TC.S.fpupdate_TIME);

		didcellchange = false;

		r_inner2 = TC.ColPars[5*TC.collider]*TC.ColPars[5*TC.collider];
		r_outer2 = TC.ColPars[5*TC.collider + 2]*TC.ColPars[5*TC.collider + 2];
		epsilon_now = TC.ColPars[5*TC.collider + 1];
		epsilon_inner = TC.ColPars[5*TC.collider + 3];
		epsilon_outer = TC.ColPars[5*TC.collider + 4];
		if(TC.Event_type[TC.collider] == 2)							
		{
			potential_inner = -epsilon_now + epsilon_inner;
			potential_outer = -epsilon_now + epsilon_outer;
		}
		if(TC.Event_type[TC.collider] == 1)
		{
			//because bond potential is a square well with zero depth
			potential_inner = 1.0e08;							
			potential_outer = 1.0e08;
		}

		calculator.coordinate = TC.S.P[TC.collider].coordinate - TC.S.P[TC.col_partner].coordinate;
		calculator.velocity = TC.S.P[TC.collider].velocity - TC.S.P[TC.col_partner].velocity;
//		calculator.coordinate=TC.S.OneParticlePositionUpdater(calculator, TC.S.TIME, TC.S.fpupdate_TIME);
		TC.S.min_img(calculator.coordinate, L);
		r2 = calculator.coordinate.norm2();
		reduced_mass = TC.S.P[TC.collider].mass*TC.S.P[TC.col_partner].mass/(TC.S.P[TC.collider].mass + TC.S.P[TC.col_partner].mass);
		//Square well collision
		if(TC.Event_type[TC.collider] == 2)
		{//Hard collision on the left side
			if(TC.Collision_type[TC.collider] == 1)
			{
//Necessary to include this for particles colliding on edges of the box//Reference: Alder, Wainwright 1959
				b = recalc_b(calculator, L, TC.Collision_type[TC.collider], r_inner2, r_outer2);
				dpx = -2*reduced_mass*b*calculator.coordinate.x/r2;
				dpy = -2*reduced_mass*b*calculator.coordinate.y/r2;
				dpz = -2*reduced_mass*b*calculator.coordinate.z/r2;
			}
			//Transfer on the left side
			if(TC.Collision_type[TC.collider] == 2)
			{//Necessary to include this for particles colliding on edges of the box
				b = recalc_b(calculator, L, TC.Collision_type[TC.collider], r_inner2, r_outer2);
				D = b*b - (2*r2*potential_inner/reduced_mass);
				dpx = -reduced_mass*calculator.coordinate.x*(-sqrt(D) + b)/r2;
				dpy = -reduced_mass*calculator.coordinate.y*(-sqrt(D) + b)/r2;
				dpz = -reduced_mass*calculator.coordinate.z*(-sqrt(D) + b)/r2;
				TC.S.potential_energy = TC.S.potential_energy + potential_inner;
			}
			//Hard collision on the right side
			if(TC.Collision_type[TC.collider] == 3)
			{//Necessary to include this for particles colliding on edges of the box
				b = recalc_b(calculator, L, TC.Collision_type[TC.collider], r_inner2, r_outer2);
				dpx = -2*reduced_mass*b*calculator.coordinate.x/r2;
				dpy = -2*reduced_mass*b*calculator.coordinate.y/r2;
				dpz = -2*reduced_mass*b*calculator.coordinate.z/r2;
			}
			//Transfer on the right side
			if(TC.Collision_type[TC.collider] == 4)
			{//Necessary to include this for particles colliding on edges of the box
				b = recalc_b(calculator, L, TC.Collision_type[TC.collider], r_inner2, r_outer2);
				D = b*b - (2*r2*potential_outer/reduced_mass);
				dpx = -reduced_mass*calculator.coordinate.x*(-sqrt(D) + b)/r2;
				dpy = -reduced_mass*calculator.coordinate.y*(-sqrt(D) + b)/r2;
				dpz = -reduced_mass*calculator.coordinate.z*(-sqrt(D) + b)/r2;	
				TC.S.potential_energy = TC.S.potential_energy + potential_outer;
			}
			if(TC.Collision_type[TC.collider] == 0)
			{
				cout<<"Error: Coltype registered as 0 but the collision time is not (should not be) timbig! Collision type not stored!"<<endl;
				exit(1);
			}
			swcol_counter = swcol_counter + 1;
		}
		//Bond collision
		else if(TC.Event_type[TC.collider] == 1)
		{//Reference: Alder, Wainwright 1959
			if(TC.Collision_type[TC.collider] == 1)
			{//Necessary to include this for particles colliding on edges of the box
				b = recalc_b(calculator, L, TC.Collision_type[TC.collider], r_inner2, r_outer2);
				dpx = -2*reduced_mass*b*calculator.coordinate.x/r2;
				dpy = -2*reduced_mass*b*calculator.coordinate.y/r2;
				dpz = -2*reduced_mass*b*calculator.coordinate.z/r2;
			}
			if(TC.Collision_type[TC.collider] == 2 || TC.Collision_type[TC.collider] == 4)
			{
				cout<<"Type 2 collision in bond mechanism, not possible so exiting"<<endl;
				cout<<"Particle 1 x = "<<TC.S.P[TC.collider].coordinate.x<<"\t y = "<<TC.S.P[TC.collider].coordinate.y<<"\t z = "<<TC.S.P[TC.collider].coordinate.z<<"\t i ="<<TC.collider<<endl;
				cout<<"Particle 2 x = "<<TC.S.P[TC.col_partner].coordinate.x<<"\t y = "<<TC.S.P[TC.col_partner].coordinate.y<<"\t z = "<<TC.S.P[TC.col_partner].coordinate.z<<"\t j ="<<TC.col_partner<<endl;
				cout<<"Particle 1 vx = "<<TC.S.P[TC.collider].velocity.vx<<"\t vy = "<<TC.S.P[TC.collider].velocity.vy<<"\t vz = "<<TC.S.P[TC.collider].velocity.vz<<endl;
				cout<<"Particle 2 vx = "<<TC.S.P[TC.col_partner].velocity.vx<<"\t vy = "<<TC.S.P[TC.col_partner].velocity.vy<<"\t vz = "<<TC.S.P[TC.col_partner].velocity.vz<<endl;
				exit(1);
			}
			if(TC.Collision_type[TC.collider] == 3)
			{//Necessary to include this for particles colliding on edges of the box
				b = recalc_b(calculator, L, TC.Collision_type[TC.collider], r_inner2, r_outer2);
				dpx = -2*reduced_mass*b*calculator.coordinate.x/r2;
				dpy = -2*reduced_mass*b*calculator.coordinate.y/r2;
				dpz = -2*reduced_mass*b*calculator.coordinate.z/r2;
			}
			bondcol_counter = bondcol_counter + 1;
		}
		else
		{
			cout<<"Wrong Event type assigned"<<endl;
			exit(1);
		}
		
		virial = virial+(dpx*calculator.coordinate.x+dpy*calculator.coordinate.y+ dpz*calculator.coordinate.z)/reduced_mass;
		//New velocities of the collided particles	
		TC.S.P[TC.collider].velocity.vx = TC.S.P[TC.collider].velocity.vx + (dpx/TC.S.P[TC.collider].mass);
		TC.S.P[TC.collider].velocity.vy = TC.S.P[TC.collider].velocity.vy + (dpy/TC.S.P[TC.collider].mass);	
		TC.S.P[TC.collider].velocity.vz = TC.S.P[TC.collider].velocity.vz + (dpz/TC.S.P[TC.collider].mass);
		TC.S.P[TC.col_partner].velocity.vx = TC.S.P[TC.col_partner].velocity.vx - (dpx/TC.S.P[TC.col_partner].mass);
		TC.S.P[TC.col_partner].velocity.vy = TC.S.P[TC.col_partner].velocity.vy - (dpy/TC.S.P[TC.col_partner].mass);
		TC.S.P[TC.col_partner].velocity.vz = TC.S.P[TC.col_partner].velocity.vz - (dpz/TC.S.P[TC.col_partner].mass);
		TC.S.P[TC.collider].velocity2 = TC.S.P[TC.collider].velocity.norm2();
		TC.S.P[TC.col_partner].velocity2 = TC.S.P[TC.col_partner].velocity.norm2();
		//Moving the particles backwards in time with the new velocities given
		TC.S.P[TC.collider].coordinate=TC.S.OneParticlePositionBackwarder(TC.S.P[TC.collider], TC.S.TIME, TC.S.fpupdate_TIME);
		TC.S.P[TC.col_partner].coordinate=TC.S.OneParticlePositionBackwarder(TC.S.P[TC.col_partner], TC.S.TIME, TC.S.fpupdate_TIME);

//		cout<<"Particle 1 x = "<<TC.S.P[TC.collider].coordinate.x<<"\t y = "<<TC.S.P[TC.collider].coordinate.y<<"\t z = "<<TC.S.P[TC.collider].coordinate.z<<"\t i ="<<TC.collider<<endl;
//		cout<<"Particle 2 x = "<<TC.S.P[TC.col_partner].coordinate.x<<"\t y = "<<TC.S.P[TC.col_partner].coordinate.y<<"\t z = "<<TC.S.P[TC.col_partner].coordinate.z<<"\t j ="<<TC.col_partner<<endl;
//		cout<<"Particle 1 vx = "<<TC.S.P[TC.collider].velocity.vx<<"\t vy = "<<TC.S.P[TC.collider].velocity.vy<<"\t vz = "<<TC.S.P[TC.collider].velocity.vz<<endl;
//		cout<<"Particle 2 vx = "<<TC.S.P[TC.col_partner].velocity.vx<<"\t vy = "<<TC.S.P[TC.col_partner].velocity.vy<<"\t vz = "<<TC.S.P[TC.col_partner].velocity.vz<<endl;
//		cout<<"col_type = "<<TC.Collision_type[TC.collider]<<"\t event_type ="<<TC.Event_type[TC.collider]<<"\t t_min == "<<TC.minimum_time<<endl;
//		cout<<"dvx = "<<dvx<<"\t dvy = "<<dvy<<"\t dvz = "<<dvz<<"\t D = "<<D<<"\t b = "<<b<<endl;
	}
	//Checking if the max velocity has moved by crossing the neighborlist range
	if(TC.S.neighborlist_counter)
	{
		nbrmove_distance = nbrmove_distance + TC.S.maxvel*TC.minimum_time;
		if(nbrmove_distance >= (TC.S.neighborlist_factor - 1.05)*TC.S.max_well_width)
		{
			TC.S.NeighborList();
			didnbrchange = true;
			nbrmove_counter = nbrmove_counter + 1;
			nbrmove_distance = 0.0;
//			cout<<"Doing a neighbor list recalculation due to max distance moved"<<endl;
		}
		else
			{didnbrchange = false;}
	}
}	

void Collision::AllParticlePositionUpdater(int N, vector<Particle> &P, double TIME, double fpupdate_TIME)
{
	//Move all the particles to the new positions
	for(int count=0; count<N; count++)
	{
//		cout<<"i="<<count<<"\tX= "<<P[count].coordinate.x<<"\tY= "<<P[count].coordinate.y<<"\tZ= ";
//		cout<<P[count].coordinate.z<<"\tVX= "<<P[count].velocity.vx<<"\tVY= "<<P[count].velocity.vy;
//		cout<<"\tVZ= "<<P[count].velocity.vz<<endl;

		P[count].coordinate=TC.S.OneParticlePositionUpdater(P[count], TIME, fpupdate_TIME);

//		cout<<"i="<<count<<"\tX= "<<P[count].coordinate.x<<"\tY= "<<P[count].coordinate.y<<"\tZ= ";
//		cout<<P[count].coordinate.z<<"\tVX= "<<P[count].velocity.vx<<"\tVY= "<<P[count].velocity.vy;
//		cout<<"\tVZ= "<<P[count].velocity.vz<<endl;

	}
}
