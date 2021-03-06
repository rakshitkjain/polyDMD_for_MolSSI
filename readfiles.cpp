//MD: Reads the input data from different types of files (Revision Date: 31st March 2021)
#include "readfiles.h"
//Reading PDB files 
void ReadFiles::ReadPDBFile(std::string FileName, vector<Particle> &P, vector<double> &parameters)
{
	Particle p;
	std::string dumps1, dumps2, dumps, line;
	double a, b, c, gamma, alpha, beta;

	int cryst_line = 0, remark_line = 0, test;
	int particle_number;
//T0 check if the first two lines are of Remark or CRYST1
	ifstream in;
	in.open(FileName);
		std::getline(in,dumps1);
		std::getline(in,dumps2);
	in.close();
	std::string remark ="REMARK";
	std::string cryst="CRYST1";
	size_t remark_position= dumps1.find(remark);
//This means that there is no 'REMARK' in 1st line line
	if(remark_position == std::string::npos)			
	{
		remark_position= dumps2.find(remark);
		//This means that there is a 'REMARK' in 2nd line
		if(remark_position != std::string::npos)		
		{
			remark_line = 2;
			size_t cryst_position = dumps1.find(cryst);
			if(cryst_position != std::string::npos)			//1st line is the cryst line
				{cryst_line = 1;}
			else
				{cryst_line = 0;}				//No cryst line
		}
		else
		{
			remark_line = 0;					//No remark line
			size_t cryst_position = dumps1.find(cryst);
			if(cryst_position != std::string::npos)			//1st line is the cryst line
				{cryst_line = 1;}
			else
			{
				size_t cryst_position = dumps2.find(cryst);
				if(cryst_position != std::string::npos)		//2nd line is the cryst line
					{cryst_line = 2;}
				else
					{cryst_line = 0;}			//No cryst line
			}
		}
	}
//Remark is in the 1st line
	else							
	{
		remark_line = 1;
		size_t cryst_position = dumps2.find(cryst);
		if(cryst_position != std::string::npos)				//2nd line is the cryst line
			{cryst_line = 2;}
		else
			{cryst_line = 0;}					//No cryst line
	}
	try 
	{
		if(cryst_line == 1)
		{
			in.open(FileName);
				in>>dumps>>a>>b>>c>>alpha>>beta>>gamma;
			in.close();
			throw 1;
		}
		else if(cryst_line == 2)
		{
			in.open(FileName);
				std::getline(in,dumps);
				in>>dumps>>a>>b>>c>>alpha>>beta>>gamma;
			in.close();
			throw 2;
		}
		else if(cryst_line == 0)
			{throw 0;}
	}
	catch(int x)
	{
		test = x;
		if(x!=0)
		{
			if(a <= parameters[0] && b <= parameters[1] && c <= parameters[2] && alpha == parameters[3] && beta == parameters[4] && gamma == parameters[5])
				{cout<<"Everything is fine, start reading the main file"<<endl;}
			else
				{cout<<"PDB file reading error, You gon done a mistake in the CRYST1 field"<<endl;}
		}
	}
//Start to read the main file now
//Now get the remark/cryst lines from the PDB
	in.open(FileName);
		if((test == 2 && remark_line == 0) || (test==1 && remark_line ==2))
		{
			std::getline(in,dumps);
			std::getline(in,dumps);
		}
		else if((test ==1 && remark_line == 0) || (test==0 && remark_line ==1) )
		{
			std::getline(in,dumps);
		}
		int count = 0;	
		while(!in.eof())
		{
			std::getline(in,line);
			//Now start reading the lines for the data
			if(line.find("ATOM") != std::string::npos )
			{
				p.type = line.substr(0,4);
				std::istringstream particle_numberstream(line.substr(6,5));
				particle_numberstream >> particle_number;
				p.name = line.substr(12,4);
				p.resname = line.substr(17,3);
				std::istringstream chaintypestream(line.substr(21,1));
				chaintypestream >> p.chaintype;
				std::istringstream chainnumberstream(line.substr(22,4));
				chainnumberstream >> p.chainnumber;
				std::istringstream coordinatestream(line.substr(30,24));
				coordinatestream >> p.coordinate.x >> p.coordinate.y >> p.coordinate.z;
				std::istringstream occupancystream(line.substr(54,6));
				occupancystream >> p.occupancy;
				std::istringstream diameterstream(line.substr(60,6));
				diameterstream >> p.diameter;
				p.moltype = line.substr(72,4);
				p.symbol = line.substr(76,2);
				if(count + 1 == particle_number)
				{
					P.push_back(p);
				}
				else	
				{
					cerr<<"You done goofed up in reading particle number and for particle number = "<<particle_number<<endl;
					exit(1);
				}
			}
			else if(line.find("HETATM") != std::string::npos)
			{
				p.type = line.substr(0,6);
				std::istringstream particle_numberstream(line.substr(8,5));
				particle_numberstream >> particle_number;
				p.name = line.substr(14,4);
				p.resname = line.substr(19,3);
				std::istringstream chaintypestream(line.substr(23,1));
				chaintypestream >> p.chaintype;
				std::istringstream chainnumberstream(line.substr(24,4));
				chainnumberstream >> p.chainnumber;
				std::istringstream coordinatestream(line.substr(32,24));
				coordinatestream >> p.coordinate.x >> p.coordinate.y >> p.coordinate.z;
				std::istringstream occupancystream(line.substr(56,6));
				occupancystream >> p.occupancy;
				std::istringstream diameterstream(line.substr(62,6));
				diameterstream >> p.diameter;
				p.moltype = line.substr(74,4);
				p.symbol = line.substr(78,2);
				if(count + 1 == particle_number)
				{
					P.push_back(p);
				}
				else	
				{
					cerr<<"You done goofed up in reading particle number and for particle number = "<<particle_number<<endl;
					exit(1);
				}
			}
			count = count+1;
		}
	in.close();
	cout<<"Read the pdb file completely"<<endl;
}
//Reading PSF files for bonds, pseudobonds, NOTHING FOR DIHEDRAL RIGHT NOW, charge, mass and names for parameter file
void ReadFiles::ReadPSFFile(std::string FileName, vector<Particle> &P, vector<BOND> &bondlist, vector<NONBOND> &nonbondlist)
{
	std::string line;
	bool natom_found=false, nbond_found=false, anglebond_found = false, atomtype_found = false, diditfeed;
	bool swapped;
	int n_residue, n_bonds, n_angles, bond_counter = 0, bond_number, angle_counter = 0, angle_number;
	int bond[8], anglebonds[9];			//Taking in multiple bonds and pseudobonds from angles
	BOND feed_bond, feed_pseudo, bond_orderer;	//feeder for covalent bond and pseudobonds
	NONBOND filler;
//This is needed bcoz vector::insert needs an iterator variable
	auto iterator=bondlist.begin();
	double charge, mass;
	int particle_number = 0, chainnumber, dump;
	std::string name, moltype, resname, atomtype;

	ifstream in;
	in.open(FileName);
		std::getline(in,line);
		while(!in.eof())
		{
			//Atom records not found condition
			while((line.find("!NATOM")) == std::string::npos)			
			{
				std::getline(in,line);
			}
			natom_found = true;
			std::istringstream linestream(line);
			linestream >> n_residue;			//Number of atom records in the file
			while(particle_number < n_residue)
			{
				in >> particle_number >> moltype >> chainnumber >> resname >> name >> atomtype >> charge >> mass >> dump;
				if(P[particle_number - 1].chainnumber != chainnumber)
				{
					cout<<"Error when reading chainnumber for particle number"<<particle_number<<endl;
					exit(1);
				}
				else if(P[particle_number - 1].resname != resname)
				{
					cout<<"Error when reading resname for particle number"<<particle_number<<endl;
					exit(1);
				}
				//This means all the entries are fine in the psf file, matching with the pdb file, so take data from psf
				else		
				{
					P[particle_number-1].charge = charge;
					P[particle_number-1].mass = mass;
					P[particle_number-1].atomtype.replace(P[particle_number-1].atomtype.begin(), P[particle_number-1].atomtype.end(), atomtype);
				}	
			}
			//Now start to collect bond records
			//Bond records not found condition
			while((line.find("!NBOND")) == std::string::npos)			
			{
				std::getline(in,line);
			}
			nbond_found = true;
			std::istringstream linestream1(line);
			linestream1 >> n_bonds;
			bond_number = n_bonds - (n_bonds%4);
			//Now, each line till bond_number is reached, will have 4 bonds 
			while(bond_counter < bond_number)			
			{
				in>>bond[0]>>bond[1]>>bond[2]>>bond[3]>>bond[4]>>bond[5]>>bond[6]>>bond[7];
				for(int i = 0; i<8; i=i+2)
				{
					if(bond[i]<bond[i+1])
					{
						//Sorting them into smaller particle first
						feed_bond.partner1 = bond[i]-1;		
						feed_bond.partner2 = bond[i+1]-1;
						feed_bond.middle = -1;
					}
					else if(bond[i] > bond[i+1])
					{
						feed_bond.partner1 = bond[i + 1]-1;
						feed_bond.partner2 = bond[i]-1;
						feed_bond.middle = -1;
					}
					bondlist.push_back(feed_bond);
				}
				bond_counter = bond_counter + 4;
			}
			//Now, take whatever amount of bonds remain in the file
			//Will not enter here if n_bonds%4==0
			while(bond_counter < n_bonds)			
			{
				//Means only one bond left
				if(n_bonds%4 ==1)
					{in>>bond[0]>>bond[1];}
				//Means only two bonds left
				else if(n_bonds%4 ==2)
					{in>>bond[0]>>bond[1]>>bond[2]>>bond[3];}
				//Means only three bonds left
				else if(n_bonds%4 ==3)
					{in>>bond[0]>>bond[1]>>bond[2]>>bond[3]>>bond[4]>>bond[5];}
				else
				{	
					cout<<"Some erroring in reading the bonds from psf file, counter should NOT be here!!"<<endl;
					exit(1);
				}
				for(int i = 0; i<2*(n_bonds%4); i = i+2)
				{
					if(bond[i]<bond[i+1])
					{
						feed_bond.partner1 = bond[i]-1;
						feed_bond.partner2 = bond[i+1]-1;
						feed_bond.middle = -1;
					}
					else if(bond[i] > bond[i+1])
					{
						feed_bond.partner1 = bond[i + 1]-1;
						feed_bond.partner2 = bond[i]-1;
						feed_bond.middle = -1;
					}
					bondlist.push_back(feed_bond);
					bond_counter = bond_counter + 1;
				}
			}
			//Now read the angles as pseudobonds and use the values for parameters iff you want to keep it fixed
			//Angle records not found condition
			while((line.find("!NTHETA")) == std::string::npos)			
			{
				std::getline(in,line);
			}

			anglebond_found = true;
			std::istringstream linestream2(line);
			linestream2 >> n_angles;
			angle_number = n_angles - (n_angles%3);
			while(angle_counter < angle_number)
			{
				in>>anglebonds[0]>>anglebonds[1]>>anglebonds[2]>>anglebonds[3]>>anglebonds[4]>>anglebonds[5]>>anglebonds[6]>>anglebonds[7]>>anglebonds[8];
				for(int i = 0; i<9; i=i+3)
				{
					//Placing the anglebond near previous bond occurence for that first atom, makes it easy to order the bondlist
					diditfeed = false;
					if(anglebonds[i]<anglebonds[i+2])
					{
						feed_pseudo.partner1 = anglebonds[i]-1;
						feed_pseudo.partner2 = anglebonds[i+2]-1;
						feed_pseudo.middle = anglebonds[i+1]-1;
						for(int q=0; q<bondlist.size(); q++)
						{
							if(bondlist[q].partner1==anglebonds[i]-1)
							{
								iterator=bondlist.begin();
								bondlist.insert(iterator+q, feed_pseudo);
								diditfeed=true;
								break;
							}
						}
					}
					else if(anglebonds[i] > anglebonds[i+2])
					{
						feed_pseudo.partner1 = anglebonds[i + 2]-1;
						feed_pseudo.partner2 = anglebonds[i]-1;
						feed_pseudo.middle = anglebonds[i+1]-1;
						for(int q=0; q<bondlist.size(); q++)
						{
							if(bondlist[q].partner1==anglebonds[i+2]-1)
							{
								iterator=bondlist.begin();
								bondlist.insert(iterator+q, feed_pseudo);
								diditfeed=true;
								break;
							}
						}
					}
					if(!diditfeed)
						{bondlist.push_back(feed_pseudo);}
				}
				angle_counter = angle_counter +3;
			}
			while(angle_counter < n_angles)
			{
				if(n_angles%3 == 1)
					{in>>anglebonds[0]>>anglebonds[1]>>anglebonds[2];}
				else if (n_angles%3 == 2)
					{in>>anglebonds[0]>>anglebonds[1]>>anglebonds[2]>>anglebonds[3]>> anglebonds[4]>>anglebonds[5];}
				else
				{	
					cout<<"Some erroring in reading the angles from psf file, counter should NOT be here!!"<<endl;
					exit(1);
				}
				for(int i = 0; i<3*(n_angles%3); i=i+3)
				{
					//Placing the anglebond near previous bond occurence for that first atom, makes it easy to order the bondlist
					diditfeed = false;
					if(anglebonds[i]<anglebonds[i+2])
					{
						feed_pseudo.partner1 = anglebonds[i]-1;
						feed_pseudo.partner2 = anglebonds[i+2]-1;
						feed_pseudo.middle = anglebonds[i+1]-1;
						for(int q=0; q<bondlist.size(); q++)
						{
							if(bondlist[q].partner1==anglebonds[i]-1)
							{
								iterator=bondlist.begin();
								bondlist.insert(iterator+q, feed_pseudo);
								diditfeed=true;
								break;
							}
						}
					}
					else if(anglebonds[i] > anglebonds[i+2])
					{
						feed_pseudo.partner1 = anglebonds[i + 2]-1;
						feed_pseudo.partner2 = anglebonds[i]-1;
						feed_pseudo.middle = anglebonds[i+1]-1;
						for(int q=0; q<bondlist.size(); q++)
						{
							if(bondlist[q].partner1==anglebonds[i+2]-1)
							{
								iterator=bondlist.begin();
								bondlist.insert(iterator+q, feed_pseudo);
								diditfeed=true;
								break;
							}
						}
					}
					if(!diditfeed)
						{bondlist.push_back(feed_pseudo);}
					angle_counter = angle_counter +1;
				}
			}
			break;
		}
		if(!natom_found || !nbond_found || !anglebond_found)
		{
			cout<<"You did a booboo in natom records in psf"<<endl;
			exit(1);
		}
	in.close();
	cout<<"Reading from psf file complete"<<endl;
//Ordering the bondlist perfectly
	for(int i=0; i<bondlist.size()-1; i++)
	{
		swapped=false;
		for(int j=0; j<bondlist.size()-i-1; j++)
		{//IF this happens, then swap
			if(bondlist[j].partner1>bondlist[j+1].partner1)
			{
				bond_orderer=bondlist[j+1];
				bondlist[j+1]=bondlist[j];
				bondlist[j]=bond_orderer;
				swapped=true;
			}
		}
		if(swapped==false)
			{break;}
	}

//making a list of different types of atomtypes in the system
	for(int i = 0; i < P.size(); i++)
	{
		atomtype.replace(atomtype.begin(),atomtype.end(), P[i].atomtype);
		for(int j = 0; j<atomtype_vector.size(); j++)
		{
			if(atomtype == atomtype_vector[j])
			{
				atomtype_found = true;
				break;
			}
			else
				{atomtype_found = false;}
		}
		if(atomtype_found == false)
		{
			atomtype_vector.push_back(atomtype);
		}
		else
			{continue;}
	}
//Initializing the nonbond list, this prevents repetition of pairs
	for(int i = 0; i < atomtype_vector.size(); i++)
	{
		for(int j = i; j < atomtype_vector.size(); j++)
		{
			filler.partner1.replace(filler.partner1.begin(), filler.partner1.end(), atomtype_vector[i]);
			filler.partner2.replace(filler.partner2.begin(), filler.partner2.end(), atomtype_vector[j]);
			nonbondlist.push_back(filler);
		}
	}
}
//Read parameter file for bondlength for covalents, bondlength for pseudos, and all parameters for nonbonded particles (INCLUDES MULTIPLE POTENTIALS AS WELL)
void ReadFiles::ReadParameterFile(std::string FileName, vector<Particle> &P, vector<BOND> &bondlist, vector<NONBOND> &nonbondlist )
{
	double L1, L2, epsilon;
	NONBOND filler;
	std::string line, dumper[4];
	std:istringstream linestream;
	bool nonbond_found = false;
	double bondlen;
	bool foundit = false;
	int counter_covalent = 0, counter_pseudo = 0;
	
	ifstream in;
	in.open(FileName);
		std::getline(in,line);
		//Bond records not found condition
		while((line.find("BONDS")) == std::string::npos)			
		{
			std::getline(in,line);
		}
		//Means that the pseudobond section has still not started
		while((line.find("PSEUDO")) == std::string::npos)		
		{
			std::istringstream linestream(line);
			linestream>>dumper[0]>>dumper[1];
			linestream>>bondlen;
			for(int i = 0; i< bondlist.size();i++)
			{
				if((P[bondlist[i].partner1].atomtype == dumper[0] && P[bondlist[i].partner2].atomtype == dumper[1]) || (P[bondlist[i].partner1].atomtype == dumper[1] && P[bondlist[i].partner2].atomtype == dumper[0]))
					{foundit = true;}
				else
					{foundit = false;}
				if(foundit == true && bondlist[i].middle == -1)
				{
					bondlist[i].bondlength = bondlen;
					counter_covalent = counter_covalent + 1;
				}
			}				
			std::getline(in,line);
		}
		//Now start with the pseudobonds section
		//Pseudobond records not found condition
		while((line.find("PSEUDO")) == std::string::npos)			
		{
			std::getline(in,line);
		}
		//Means that the nonbond section has still not started
		while((line.find("NONBOND")) == std::string::npos)		
		{
			std::istringstream linestream(line);
			linestream>>dumper[0]>>dumper[1];
			linestream>>bondlen;
			for(int i = 0; i< bondlist.size();i++)
			{
				if((P[bondlist[i].partner1].atomtype == dumper[0] && P[bondlist[i].partner2].atomtype == dumper[1]) || (P[bondlist[i].partner1].atomtype == dumper[1] && P[bondlist[i].partner2].atomtype == dumper[0]))
					{foundit = true;}
				else
					{foundit = false;}
				if(foundit == true && bondlist[i].middle >= 0)
				{
					bondlist[i].bondlength = bondlen;
					counter_pseudo = counter_pseudo + 1;
				}
			}				
			std::getline(in,line);
		}
		//MAKING SURE ALL THE BONDS AND PSEUDOBONDS ARE FILLED
		int j = 0;
		while(j< bondlist.size())
		{
			if(bondlist[j].middle == -1 && bondlist[j].bondlength == 0.0)
			{	
				cout<<"Error reading covalent bonds in parameter file"<<endl;
				cout<<"P1= "<<bondlist[j].partner1<<"P2= "<<bondlist[j].partner2<<endl;
				exit(1);
			}
			if(bondlist[j].middle >= 0 && bondlist[j].bondlength == 0.0)
			{	
				cout<<"Error reading pseudo bonds in parameter file"<<endl;
				cout<<"P1= "<<bondlist[j].partner1<<"P2= "<<bondlist[j].partner2<<endl;
				exit(1);
			}
			j = j+1;
		}
		//Now start with the nonbonded section
		//Nonbond records not found condition
		while((line.find("NONBOND")) == std::string::npos)			
		{
			std::getline(in,line);
		}
		std::getline(in,line);
		//Means that the nonbond section has still not ended
		while(line.find("END") == std::string::npos)			
		{
			std::istringstream linestream(line);
			linestream>>dumper[0]>>dumper[1];
			linestream>>L1>>L2>>epsilon;
			for(int i = 0; i< nonbondlist.size(); i++)
			{
				if((nonbondlist[i].partner1 == dumper[0] && nonbondlist[i].partner2 == dumper[1]) || (nonbondlist[i].partner1 == dumper[1] && nonbondlist[i].partner2 == dumper[0]))
				{
					//So, here, the particle pair is in the parameter file and it is nonbonded pair, so read the data
					//If multiple potentials and this is the first one
					if(nonbondlist[i].L1.size() == 0)
					{
						nonbondlist[i].L1.push_back(L1);
						nonbondlist[i].L2.push_back(L2);
						nonbondlist[i].epsilon.push_back(epsilon);
					}
					else
					{	//This makes sure that the wells are continuous
						if(nonbondlist[i].L2[nonbondlist[i].L2.size()]==L1)
						{
							nonbondlist[i].L1.push_back(L1);	
							nonbondlist[i].L2.push_back(L2);
							nonbondlist[i].epsilon.push_back(epsilon);
						}//Do nothing if repetition 
						else if(nonbondlist[i].L2[nonbondlist[i].L2.size()]==L2 && nonbondlist[i].L1[nonbondlist[i].L1.size()]==L1 && nonbondlist[i].epsilon[nonbondlist[i].epsilon.size()]==epsilon)
						{}
						else//Something weird is happening then
						{
							cout<<"Weird stuff reading from the parameter file"<<endl;
							exit(1);
						}
					}
				}
			}
			std::getline(in,line);
		}
		//MAKING SURE ALL THE NONBONDS ARE FILLED AND OF THE SAME SIZE
		j = 0;
		while(j< nonbondlist.size())
		{
			if((nonbondlist[j].L1.size() != nonbondlist[j].L2.size()) || (nonbondlist[j].L1.size() != nonbondlist[j].epsilon.size()) || (nonbondlist[j].epsilon.size() != nonbondlist[j].L2.size()))
			{
				cout<<"Exiting because L1, L2 and epsilon not of same length. Do par read again, check everything!!"<<endl;
				exit(1);
			}
			else
				{j = j+1;}
		}
	in.close();
}	
