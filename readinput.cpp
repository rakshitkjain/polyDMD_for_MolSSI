#include "readinput.h"

void ReadInput::printHelpMessage()
{
	cout<<endl;
	cout<<"Usage:"<<endl;
	cout<<"DMD <options>"<<endl;
	cout<<"  -h, help---------------->Print usage message"<<endl;
	cout<<"  -l, L------------------->Box length (default 20)"<<endl;
	cout<<"  -e, nsweep-------------->No. of sweeps (default 1000)"<<endl;
	cout<<"  -T, temperature--------->Absolute temperature (default 30.0)"<<endl;
	cout<<"  -i, maxbin-------------->No. of bins for rdf calculation (default 1000)"<<endl;
	cout<<"  -w, cellwidth----------->Size of each cell for cell-list (default 10.0)"<<endl;
	cout<<"  -f, neighborlist_factor->Factor to multiply the largest potential width (default 2.0)"<<endl;
	cout<<"  -q, bond_delta---------->Change in bond length allowed (default 0.0375)"<<endl;
	cout<<"  -a, andersen_freq------->One out of every andersen_freq collisions will be a ghost collision (default 10)"<<endl;
	cout<<"  -pf, print_freq--------->Print frequency for data, multiply by fpupdate_freq (default 0.001)"<<endl;
	cout<<"  -p, fpupdate_freq------->False position update frequency (default 1000)"<<endl;
	cout<<"  -b, burn_percent-------->Fraction of nsweep which are for burning (no thermostat, no virial, no rdf) (default 0.1)"<<endl;
	cout<<"  -pdb, filename_pdb------>PDB filename (Simul won't work without filename) (default test.pdb)"<<endl;
	cout<<"  -psf, filename_psf------>PSF filename (Simul won't work without filename) (default test.psf)"<<endl;
	cout<<"  -par, filename_par------>Parameter filename (Simul won't work without filename) (default test.par)"<<endl;
	cout<<"  -out, filename_out------>Partial output filename (Will be same as the name of pdb file if not given) (default test)"<<endl;
//	cout<<"  -n, nchains------------->No. of chains (default 8)"<<endl;
//	cout<<"  -s, sigma1-------------->Particle diameter (default 1.0)"<<endl;
//	cout<<"  -g, sigma2-------------->Square well diameter (default 2.0)"<<endl;
//	cout<<"  -m, mass---------------->Particle mass in amu (default 1.0)"<<endl;
//	cout<<"  -d, well_depth---------->Square well depth (positive means depth)  (default 1.0)"<<endl;
//	cout<<"  -b, bondlength---------->Average bond length (default 1.5)"<<endl;
//	cout<<"  -c, chainlength--------->Number of particles in chain (default 10)"<<endl;
//	cout<<"  -v, maxvel-------------->Maximum velocity allowed, to calculate nbrlist update frequency (default 10.0)"<<endl;
}

int ReadInput::ReadVariables(int argc, vector<string> arguments)//reading x, y and z coordinates as well as velocities
{
//Writing the default values of the parameters, which we use if no flag specified
	initialize.L = 20; initialize.nsweep = 1000; initialize.maxbin = 1000; initialize.andersen_freq = 10; 
	initialize.print_freq = 0.001; initialize.fpupdate_freq = 1000; initialize.burn_percent = 0.1;
	initialize.temperature = 30.0;	initialize.bond_delta = 0.0375; initialize.cellwidth = 10.0; 
	initialize.neighborlist_factor = 2.0; initialize.maxvel = 10.0;
	initialize.filename_pdb = "test.pdb"; initialize.filename_psf = "test.psf"; 
	initialize.filename_par = "test.par"; initialize.filename_out = "test";
	bool pdb = false, psf = false, out = false, par = false;
//Write the arguments to a string vector
//	vector<string> arguments;
//	for(int i = 0; i< argc; i++)
//		{arguments.push_back(string(argv[i]));}
	for(size_t i =0; i<arguments.size(); i++)
	{
		//Means this is not the only thing in the vector
		if(arguments[i].find("-") != arguments[i].npos)					
		{
			if(arguments[i] == "-h")
			{
				printHelpMessage();
				exit(1);
				return 0;
			}
			else if(arguments[i] == "-l")
			{
			//Checking if its the last argument or if the next is not a value for the variable
				if(i+1>=arguments.size() || arguments[i+1].find("-") != arguments[i+1].npos)
				{
					cerr<<"Value not written after the flag so exiting"<<endl;
					exit(1);
					return -1;
				}
				else
				{
					stringstream Lstream(arguments[i+1]);
					Lstream >> initialize.L;
					if(!Lstream && !Lstream.eof())
					{
						cout<<"No value written, write value after using a flag"<<endl;
						exit(1);
						return -1;
					}
				}
			}
			else if(arguments[i] == "-e")
			{
			//Checking if its the last argument or if the next is not a value for the variable
				if(i+1>=arguments.size() || arguments[i+1].find("-") != arguments[i+1].npos)
				{
					cerr<<"Value must be written after the flag, exiting due to no value returned"<<endl;
					exit(1);
					return -1;
				}
				else
				{
					stringstream nsweepstream(arguments[i+1]);
					nsweepstream >> initialize.nsweep;
					if(!nsweepstream && !nsweepstream.eof())
					{
						cout<<"No value written, write value after using a flag"<<endl;
						exit(1);
						return -1;
					}
				}
			}
			else if(arguments[i] == "-T")
			{
			//Checking if its the last argument or if the next is not a value for the variable
				if(i+1>=arguments.size() || arguments[i+1].find("-") != arguments[i+1].npos)
				{
					cerr<<"Value not written after the flag so exiting"<<endl;
					exit(1);
					return -1;
				}
				else
				{
					stringstream temperaturestream(arguments[i+1]);
					temperaturestream >> initialize.temperature;
					if(!temperaturestream && !temperaturestream.eof())
					{
						cout<<"No value written, write value after using a flag"<<endl;
						exit(1);
						return -1;
					}
				}
			}
			else if(arguments[i] == "-i")
			{
			//Checking if its the last argument or if the next is not a value for the variable
				if(i+1>=arguments.size() || arguments[i+1].find("-") != arguments[i+1].npos)
				{
					cerr<<"Value not written after the flag so exiting"<<endl;
					exit(1);
					return -1;
				}
				else
				{
					stringstream maxbinstream(arguments[i+1]);
					maxbinstream >> initialize.maxbin;
					if(!maxbinstream && !maxbinstream.eof())
					{
						cout<<"No value written, write value after using a flag"<<endl;
						exit(1);
						return -1;
					}
				}
			}
			else if(arguments[i] == "-w")
			{
			//Checking if its the last argument or if the next is not a value for the variable
				if(i+1>=arguments.size() || arguments[i+1].find("-") != arguments[i+1].npos)
				{
					cerr<<"Value not written after the flag so exiting"<<endl;
					exit(1);
					return -1;
				}
				else
				{
					stringstream cellwidthstream(arguments[i+1]);
					cellwidthstream >> initialize.cellwidth;
					if(!cellwidthstream && !cellwidthstream.eof())
					{
						cout<<"No value written, write value after using a flag"<<endl;
						exit(1);
						return -1;
					}
				}
			}
			else if(arguments[i] == "-f")
			{
			//Checking if its the last argument or if the next is not a value for the variable
				if(i+1>=arguments.size() || arguments[i+1].find("-") != arguments[i+1].npos)
				{
					cerr<<"Value not written after the flag so exiting"<<endl;
					exit(1);
					return -1;
				}
				else
				{
					stringstream neighborlist_factorstream(arguments[i+1]);
					neighborlist_factorstream >> initialize.neighborlist_factor;
					if(!neighborlist_factorstream && !neighborlist_factorstream.eof())
					{
						cout<<"No value written, write value after using a flag"<<endl;
						exit(1);
						return -1;
					}
				}
			}
			else if(arguments[i] == "-a")
			{
			//Checking if its the last argument or if the next is not a value for the variable
				if(i+1>=arguments.size() || arguments[i+1].find("-") != arguments[i+1].npos)
				{
					cerr<<"Value not written after the flag so exiting"<<endl;
					exit(1);
					return -1;
				}
				else
				{
					stringstream andersen_freqstream(arguments[i+1]);
					andersen_freqstream >> initialize.andersen_freq;
					if(!andersen_freqstream && !andersen_freqstream.eof())
					{
						cout<<"No value written, write value after using a flag"<<endl;
						exit(1);
						return -1;
					}
				}
			}
			else if(arguments[i] == "-q")
			{
			//Checking if its the last argument or if the next is not a value for the variable
				if(i+1>=arguments.size() || arguments[i+1].find("-") != arguments[i+1].npos)
				{
					cerr<<"Value not written after the flag so exiting"<<endl;
					exit(1);
					return -1;
				}
				else
				{
					stringstream bond_deltastream(arguments[i+1]);
					bond_deltastream >> initialize.bond_delta;
					if(!bond_deltastream && !bond_deltastream.eof())
					{
						cout<<"No value written, write value after using a flag"<<endl;
						exit(1);
						return -1;
					}
				}
			}
			else if(arguments[i] == "-pf")
			{
			//Checking if its the last argument or if the next is not a value for the variable
				if(i+1>=arguments.size() || arguments[i+1].find("-") != arguments[i+1].npos)
				{
					cerr<<"Value not written after the flag so exiting"<<endl;
					exit(1);
					return -1;
				}
				else
				{
					stringstream print_freqstream(arguments[i+1]);
					print_freqstream >> initialize.print_freq;
					if(!print_freqstream && !print_freqstream.eof())
					{
						cout<<"No value written, write value after using a flag"<<endl;
						exit(1);
						return -1;
					}
				}
			}
			else if(arguments[i] == "-p")
			{
			//Checking if its the last argument or if the next is not a value for the variable
				if(i+1>=arguments.size() || arguments[i+1].find("-") != arguments[i+1].npos)
				{
					cerr<<"Value not written after the flag so exiting"<<endl;
					exit(1);
					return -1;
				}
				else
				{
					stringstream fpupdate_freqstream(arguments[i+1]);
					fpupdate_freqstream >> initialize.fpupdate_freq;
					if(!fpupdate_freqstream && !fpupdate_freqstream.eof())
					{
						cout<<"No value written, write value after using a flag"<<endl;
						exit(1);
						return -1;
					}
				}
			}
			else if(arguments[i] == "-b")
			{
			//Checking if its the last argument or if the next is not a value for the variable
				if(i+1>=arguments.size() || arguments[i+1].find("-") != arguments[i+1].npos)
				{
					cerr<<"Value not written after the flag so exiting"<<endl;
					exit(1);
					return -1;
				}
				else
				{
					stringstream burn_percentstream(arguments[i+1]);
					burn_percentstream >> initialize.burn_percent;
					if(!burn_percentstream && !burn_percentstream.eof())
					{
						cout<<"No value written, write value after using a flag"<<endl;
						exit(1);
						return -1;
					}
				}
			}
			else if(arguments[i] == "-pdb")
			{
			//Checking if its the last argument or if the next is not a value for the variable
				if(i+1>=arguments.size() || arguments[i+1].find("-") != arguments[i+1].npos)
				{
					cerr<<"Value not written after the flag so exiting"<<endl;
					exit(1);
					return -1;
				}
				else
				{
					pdb = true;
					initialize.filename_pdb = arguments[i+1];
				}
			}
			else if(arguments[i] == "-psf")
			{
			//Checking if its the last argument or if the next is not a value for the variable
				if(i+1>=arguments.size() || arguments[i+1].find("-") != arguments[i+1].npos)
				{
					cerr<<"Value not written after the flag so exiting"<<endl;
					exit(1);
					return -1;
				}
				else
				{
					psf = true;
					initialize.filename_psf = arguments[i+1];
				}
			}
			else if(arguments[i] == "-par")
			{
			//Checking if its the last argument or if the next is not a value for the variable
				if(i+1>=arguments.size() || arguments[i+1].find("-") != arguments[i+1].npos)
				{
					cerr<<"Value not written after the flag so exiting"<<endl;
					exit(1);
					return -1;
				}
				else
				{
					par = true;
					initialize.filename_par = arguments[i+1];
				}
			}
			else if(arguments[i] == "-out")
			{
			//Checking if its the last argument or if the next is not a value for the variable
				if(i+1>=arguments.size() || arguments[i+1].find("-") != arguments[i+1].npos)
				{
					cerr<<"Value not written after the flag so exiting"<<endl;
					exit(1);
					return -1;
				}
				else
				{
					out = true;
					initialize.filename_out = arguments[i+1];
					cout<<"OUT filename = "<<initialize.filename_out<<endl;
				}
			}
			else
			{
				cout<<"Unrecognized flag, error error"<<endl;
				exit(1);
				return -1;
			}
		}
	}
        if(!pdb || !psf || !par)
        {
		cout<<"exiting because initialization files not there"<<endl;
		exit(1);
	}
}
