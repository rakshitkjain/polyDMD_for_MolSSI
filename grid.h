//KMC: grid.h Grid Class (Revision Date: 6th Jan 2019)
#ifndef _GRID_H
#define _GRID_H
#include "header.h"
#include "xyz.h"
#include "vel.h"
class Grid
{
public:
	XYZ com;//coordinates
	VEL velocity_com;//velocity
	Grid(){com.x=0.0; com.y=0.0; com.z=0.0; velocity_com.vx = 0.0; velocity_com.vy = 0.0; velocity_com.vz = 0.0;}
};
#endif
