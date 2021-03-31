//MD: Nonbond Class to store and manipulate nonbonded atoms along with the potentials and ranges stored  (Revision Date: 4th January 2020)
#ifndef _NONBOND_H
#define _NONBOND_H
#include "header.h"
class NONBOND
{
	public:
		std::string partner1, partner2;
		vector<double> L1;
		vector<double> L2;
//ep<0 means attractive and ep>0 means repulsive. For shoulder, STILL NEED TO ADD FUNCTION FOR THIS
		vector<double> epsilon;				
		
	 	NONBOND(std::string _partner1, std::string _partner2, vector<double> _L1, vector<double> _L2, vector<double> _epsilon){partner1=_partner1;        partner2=_partner2;   L1=_L1;	L2=_L2;	epsilon=_epsilon;}
		NONBOND(){partner1="CG331";      partner2="CG331";}
};
#endif

//Nonbonlist is now defined using the atomtypes
