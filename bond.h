//MD: Bond Class to store and manipulate bonded atoms along with the type of bond stored  (Revision Date: 4th November 2019)
#ifndef _BOND_H
#define _BOND_H
#include "header.h"
class BOND
{
	public:
		int partner1, partner2, middle;
//		std::string type;	
		double bondlength;					
	 	BOND(int _partner1, int _partner2, int _middle, double _bondlength){partner1=_partner1;        partner2=_partner2;   middle=_middle; bondlength=_bondlength;}
		BOND(){partner1=-1;      partner2=-1;  middle=-1;	bondlength=0.0;}
};
#endif

//(NOT USED ANYMORE)type => cov == covalent bond, pseudo == angle pseudobond (to keep angles fixed)
//middle = -1 means normal bond. middle > 0  means angle bond 
//partner1 is the smaller atom number and partner2 is the bigger one
