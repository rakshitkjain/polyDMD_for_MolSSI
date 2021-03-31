//KMC: VEL Class to store and manipulate velocities  (Revision Date: 6th january 2019)
#ifndef _VEL_H
#define _VEL_H
#include "header.h"
class VEL
{
	public:
		double vx,vy,vz;
		VEL(double _vx, double _vy, double _vz){vx=_vx;        vy=_vy;   vz=_vz;}
		VEL(){vx=0.0;      vy=0.0;  vz=0.0;}
		VEL operator + (VEL& other){return VEL(vx+other.vx,vy+other.vy,vz+other.vz);}
		VEL operator - (VEL& other){return VEL(vx-other.vx,vy-other.vy,vz-other.vz);}
		VEL operator / (double s){return VEL(vx/s,vy/s,vz/s);}
		VEL operator * (double s){return VEL(vx*s,vy*s,vz*s);}

		double real_v2(VEL& other){return (vx-other.vx)*(vx-other.vx)+(vy-other.vy)*(vy-other.vy)+(vz-other.vz)*(vz-other.vz);}
		double real_v(VEL& other){return sqrt((vx-other.vx)*(vx-other.vx)+(vy-other.vy)*(vy-other.vy)+(vz-other.vz)*(vz-other.vz));}
		void my_abs(){vx=fabs(vx); vy=fabs(vy); vz=fabs(vz);}
		double norm2(){ return vx*vx+ vy*vy+ vz*vz; }
		double norm(){ return sqrt(vx*vx +vy*vy +vz*vz); }
};
#endif
