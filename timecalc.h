//Here we define the TimeCalc class function declarations, function definitions are in timecalc.cpp
//This construct ensures you do not re-define headers
#ifndef _TIMECALC_H
#define _TIMECALC_H

#include "system.h"
class TimeCalc
{
	public:
		System S;
		//For storing collision times
		vector<double> TimeList;			
		//For storing particles involved in collision
		vector<int> Partner;				
		//Stores the type of collision occuring in the system 
		vector<int> Collision_type;			
		//Stores the types of collisions for collision
		vector<int> Event_type;				
		//Collision parameters ( stores L1, L2 and potential difference for the particular collision)
		vector<double> ColPars;				
		//Random huge time
		double timbig = 1.0e5;				
		//For the minimum time of collision
		double minimum_time;				
		//For the collision count, collider and its partner
		int collider, col_partner;			
		//Tolerance to check particle overlap
		double tolerance = 1.0e-04;			
//General calculation and initialization functions
	//Function so that we don't have to write time caluculation again and again
	void CollisionTime_ij(int, int, int, double, double, double, double, double);		
	//Overloading the function to reduce code work
	void CollisionTime_ij(int, int, int, vector<double>, vector<double>, vector<double>);	
	//Initialising time array and associated vectors
	void timearray_initialization(int);							
//For not cell list calculations
	//Updating upward list
	void UpdateUpList(int, int);							
	//Updating downward list
	void UpdateDownList(int, int);							
	//Finding the soonest collision time
	void Collision_time(double, int);						
//Cell List functions
	//Time calculation method if we use cell-list
	void Collision_time_celllist(double, int);					
	//Time for soonest cell change event for this particle
	double CellChangeTime(int, double, int, int, int);				
	//UpdateCellList for the case of cell-list implementation
	void UpdateUpCellList(int, int);						
	//UpdateDownList for cell list implementation
	void UpdateDownCellList(int, int);						
//Neighborlist functions
	void Collision_time_neighborlist(double, int);
	void UpdateUpNeighborList(int, int);
	void UpdateDownNeighborList(int, int);
};
#endif

//OLD__________________________Collision_type: 1= core collision, 2 = Capture (particles outside well, and then come inside) 3 = Dissociation (Opp. of 2) 4 = Bounce, meaning the KE of molecules is not sufficient to get out of the well, so it will just get bounced from the well boundary
//OLD__________________________Event_type: 1=bond collision, 2=Well event, 3= Cell change event

//NEW Collision_type: 1- Bounce on left side, 2- Transfer on left side (inner side), 3-Bounce on right side, 4-Transfer on right side (outer side), 0- Particles are going away 
// LEFT IS INNER, RIGHT IS OUTER
//Event_type:  1 = bond collision 2 = Well event 3 = Cell change event
//For Event_type=1, we should only ever have 1 or 3 Collision_type

