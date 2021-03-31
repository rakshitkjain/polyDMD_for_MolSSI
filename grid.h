//MD: Grid Class (Revision Date: 31st March 2021)
#ifndef _GRID_H
#define _GRID_H
#include "header.h"
#include "xyz.h"
#include "vel.h"

class Grid
{
	public:
//Import from other classes
		XYZ com;			//coordinates
		VEL velocity_com;		//velocity
//Class variables
	Grid(){com.x=0.0; com.y=0.0; com.z=0.0; velocity_com.vx = 0.0; velocity_com.vy = 0.0; velocity_com.vz = 0.0;}
};
#endif
