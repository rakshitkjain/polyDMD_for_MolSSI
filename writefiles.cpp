#include "writefiles.h"

void WriteFiles::DCDHeader(double Start_timestep, int N, std::string FileName)
{
	int zero=0, twentyfour=24, onesixtyfour=164, two=2;
	int number = 84;
	int num_timesteps=0;	//number of timesteps stored in dcd file, right now 0 
	int time_bw = 1;	//Time between saves, just storing for the sake of it (Maybe make DCD freq by time spent?)
	float delta = 1; 	//Trajectory timestep
	int nlines = 2;		//Number of 80 byte lines for the header
	char cord[4] = {'C','O','R','D'};

	std::ostringstream fnd;//FileName dummy
        fnd<<FileName<<".dcdheader";
        std::string FileName1 = fnd.str();

        ofstream outputFile;
        outputFile.open(FileName1, std::ios::out | std::ios::binary);

	outputFile.write(reinterpret_cast<const char*>(&number), sizeof(int));
//	for(int i = 0; i< 4; i++)
//	{
	outputFile.write(reinterpret_cast<const char*>(&cord), sizeof(char));
//	}
	outputFile.write(reinterpret_cast<const char*>(&zero), sizeof(int));
//Start_timestep: when you start recording the simulation (ISTART)
	outputFile.write(reinterpret_cast<const char*>(&Start_timestep), sizeof(float));	
	outputFile.write(reinterpret_cast<const char*>(&time_bw), sizeof(int));
	outputFile.write(reinterpret_cast<const char*>(&num_timesteps), sizeof(int));
//	outputFile.write(reinterpret_cast<const char*>(&dcd_sets), sizeof(int));
//	outputFile.write(reinterpret_cast<const char*>(&number), sizeof(int));

	for(int i = 0; i< 5; i++)
	{
		outputFile.write(reinterpret_cast<const char*>(&zero), sizeof(int));
	}
//making it a CHARMM format file
	outputFile.write(reinterpret_cast<const char*>(&delta), sizeof(float));
	outputFile.write(reinterpret_cast<const char*>(&zero), sizeof(int));
	for(int i = 0; i< 8; i++)
	{
		outputFile.write(reinterpret_cast<const char*>(&zero), sizeof(int));
	}
	outputFile.write(reinterpret_cast<const char*>(&twentyfour), sizeof(int));
	outputFile.write(reinterpret_cast<const char*>(&number), sizeof(int));
	outputFile.write(reinterpret_cast<const char*>(&onesixtyfour), sizeof(int));
	outputFile.write(reinterpret_cast<const char*>(&two), sizeof(int));
//      for(int i = 0; i< 9; i++)
//      {
//              outputFile.write(reinterpret_cast<const char*>(&zero), sizeof(int));
//      }
//      outputFile.write(reinterpret_cast<const char*>(&number), sizeof(int));

//	number = 164;
//      outputFile.write(reinterpret_cast<const char*>(&number), sizeof(int));
//      outputFile.write(reinterpret_cast<const char*>(&nlines), sizeof(int));
        for(int i = 0; i< 84; i++)
        {
                outputFile.write(reinterpret_cast<char*>(&zero), sizeof(char));
        }
	outputFile.write(reinterpret_cast<const char*>(&number), sizeof(int));

	number = 4;
        outputFile.write(reinterpret_cast<const char*>(&number), sizeof(int));
        outputFile.write(reinterpret_cast<const char*>(&N), sizeof(int));
        outputFile.write(reinterpret_cast<const char*>(&number), sizeof(int));
	outputFile.close();
}
//add routine to read number of sets if this previously exists, then modify it so that it can extract data from an existing dcdfile and write header and stuff on its own
void WriteFiles::DCDWriteStep(vector<Particle> &P, int N, std::string FileName)
{
	std::ostringstream fnd;//FileName dummy
        fnd<<FileName<<".dcd.temp";
        std::string FileName1 = fnd.str();
        ofstream outputFile;
//	int N4 = 4*N;
	int number = 48;
	double cos_alpha, cos_beta, cos_gamma;
	cos_alpha = cos(S.cellparameters[3]*PI/180.0);
	cos_beta = cos(S.cellparameters[4]*PI/180.0);
	cos_gamma = cos(S.cellparameters[5]*PI/180.0);
        outputFile.open(FileName1, std::ios::out | std::ios::binary | std::ios::app);

	outputFile.write(reinterpret_cast<const char*>(&number), sizeof(int));
	outputFile.write(reinterpret_cast<const char*>(&S.cellparameters[0]), sizeof(int));
	outputFile.write(reinterpret_cast<const char*>(&cos_gamma), sizeof(float));
	outputFile.write(reinterpret_cast<const char*>(&S.cellparameters[1]), sizeof(int));
	outputFile.write(reinterpret_cast<const char*>(&cos_beta), sizeof(float));
	outputFile.write(reinterpret_cast<const char*>(&cos_gamma), sizeof(float));
	outputFile.write(reinterpret_cast<const char*>(&S.cellparameters[2]), sizeof(int));
	
	outputFile.write(reinterpret_cast<const char*>(&number), sizeof(int));	
	for(int i = 0; i<4; i++)
	{
		outputFile.write(reinterpret_cast<const char*>(&N), sizeof(int));	
	}
	for(int i = 0; i<N; i++)
	{
		outputFile.write(reinterpret_cast<const char*>(&P[i].coordinate.x), sizeof(float));	
	}
	for(int i = 0; i<8; i++)
	{
		outputFile.write(reinterpret_cast<const char*>(&N), sizeof(int));	
	}
	for(int i = 0; i<N; i++)
	{
		outputFile.write(reinterpret_cast<const char*>(&P[i].coordinate.y), sizeof(float));	
	}
	for(int i = 0; i<8; i++)
	{
		outputFile.write(reinterpret_cast<const char*>(&N), sizeof(int));	
	}
	for(int i = 0; i<N; i++)
	{
		outputFile.write(reinterpret_cast<const char*>(&P[i].coordinate.z), sizeof(float));	
	}
	for(int i = 0; i<4; i++)
	{
		outputFile.write(reinterpret_cast<const char*>(&N), sizeof(int));	
	}
	outputFile.close();

	dcd_sets = dcd_sets + 1;
}
//Final DCD file combination
void WriteFiles::DCDFile(std::string FileName)
{
	std::ostringstream fnd;//FileName dummy
        fnd<<FileName<<".dcdheader";
        std::string FileName_header = fnd.str();

        fnd<<FileName<<".dcd.temp";
        std::string FileName_temp = fnd.str();

        fnd<<FileName<<".dcd";
        std::string FileName_dcd = fnd.str();

	std::ofstream of_dcd;	
	std::ifstream if_header(FileName_header, std::ios::binary);
	std::ifstream if_temp(FileName_temp, std::ios::binary);
	
	of_dcd.open(FileName_dcd, std::ios::binary | std::ios::app);

	of_dcd << if_header.rdbuf();
//	of_dcd.seekp(0, std::ios::end);
	of_dcd << if_temp.rdbuf();
	of_dcd.close();

	const char*header = FileName_header.c_str();
	const char*temp = FileName_temp.c_str();

	remove(header);
	remove(temp);
}
//PDB file to go along with DCD file, writing just once is fine.
void WriteFiles::PDBFile(vector<Particle> &P, std::string FileName)
{
	std::ostringstream fnd;//FileName dummy
        fnd<<FileName<<".pdb";
        std::string FileName1 = fnd.str();
        ofstream outputFile;
//Not appending right now, rewriting everytime
        outputFile.open(FileName1, std::ios::out);				

//REMARK original generated coordinate pdb file
//ATOM      1  N   CT2 A   1      34.512  35.455   4.217  1.00  0.00      CHI  N
//END
	outputFile<<"REMARK Writing PDB file for giung along with DCD"<<endl;
	for(int i = 0; i<S.N; i++)
	{
		outputFile<<setw(4)<<P[i].type<<"  "<<setw(5)<<i+1<<" "<<setw(4)<<P[i].name<<" "<<setw(3)<<P[i].resname<<" "<<P[i].chaintype<<setw(4)<<P[i].chainnumber<<"    "<<setw(8)<<setprecision(3)<<P[i].coordinate.x<<setw(8)<<setprecision(3)<<P[i].coordinate.y<<setw(8)<<setprecision(3)<<P[i].coordinate.z<<setw(6)<<setprecision(2)<<P[i].occupancy<<setw(6)<<setprecision(2)<<P[i].tempfactor<<"      "<<setw(4)<<P[i].moltype<<setw(2)<<P[i].symbol<<endl;
	}
	outputFile<<"END"<<endl;
	outputFile.close();
}
//For checking time variation of things, such as KE, etc 
void WriteFiles::TimeVarFileIni(std::string FileName)
{
	cout<<"Writing file for checking variation of parameters with time"<<endl;
	ofstream out;
	std::ostringstream fnd;							//FileName dummy
	fnd<<FileName<<".dat";
	std::string FileName1 = fnd.str();

	out.open(FileName1, ios::app);
		out<<"Time \t TE \t KE \t PE \t X momentum \t Y momentum \t Z momentum \t Temperature"<<endl;
	out.close();
}

void WriteFiles::TimeVarFile(double timestep, double ke, double pe, double te, VEL momentum, std::string FileName)
{
	ofstream out;
	std::ostringstream fnd;							//FileName dummy
	fnd<<FileName<<".dat";
	std::string FileName1 = fnd.str();

	double temperature = 2*ke/(3*S.N);
	out.open(FileName1, ios::app);
		out<<setw(8)<<timestep<<"\t"<<setw(8)<<te<<"\t"<<setw(8)<<ke<<"\t"<<setw(8)<<pe<<"\t"<<setw(8)<<momentum.vx<<"\t"<<setw(8)<<momentum.vy<<"\t"<<setw(8)<<momentum.vz<<"\t"<<setw(8)<<temperature<<endl;
	out.close();
}

void WriteFiles::RestartFile(double TIMESTEP, vector<Particle> &P, int N, std::string FileName)
{
	std::ostringstream fnd;							//FileName dummy
	fnd<<FileName<<".restart";
	std::string FileName1 = fnd.str();

      ofstream out;
      out.open(FileName1); 
	      out<<"ITEM: TIMESTEP"<<endl;
	      out<<TIMESTEP<<endl;
	      out<<"ITEM: NUMBER OF ATOMS"<<endl;
	      out<<N<<endl;
	      out<<"ITEM: BOX BOUNDS"<<endl;
	      out<<0.0<<"\t"<<S.L<<endl;
	      out<<0.0<<"\t"<<S.L<<endl;
	      out<<0.0<<"\t"<<S.L<<endl;
	      out<<"ITEM: ATOMS index type x y z"<<endl;
	      string name;
	      for(int i=0; i<N; i++)
	      {
		out<<setw(6)<<i+1<<"\t"<<setw(8)<<P[i].coordinate.x<<"\t"<<setw(8)<<P[i].coordinate.y<<"\t"<<setw(8)<<P[i].coordinate.z<<"\t"<<setw(8)<<P[i].velocity.vx;
		out<<"\t"<<setw(8)<<P[i].velocity.vy<<"\t"<<setw(8)<<P[i].velocity.vz<<endl;
	      }
      out.close();	
}
//Writing file in the LAMMPS trajectory format for particle visualisation
void WriteFiles::WriteDump(double TIMESTEP, vector<Particle> &P, int N, double L, std::string FileName)
{
	std::ostringstream fnd;//FileName dummy
	fnd<<FileName<<".dump";
	std::string FileName1 = fnd.str();
//	sprintf(FileName,"_nc_%d_cl_%d_T_%lf_L_%d_m_%lf_s1_%lf_s2_%lf_bl_%lf_del_%lf_d_%lf.dump", nchain,chnlen,temp,L,mass,sigma1,sigma2,bondlen,del,delta);      
	
	ofstream out;
      out.open(FileName1,ios::app);
 
      out<<"ITEM: TIMESTEP"<<endl;
      out<<TIMESTEP<<endl;
      out<<"ITEM: NUMBER OF ATOMS"<<endl;
      out<<N<<endl;
      out<<"ITEM: BOX BOUNDS"<<endl;
      out<<0.0<<"\t"<<L<<endl;
      out<<0.0<<"\t"<<L<<endl;
      out<<0.0<<"\t"<<L<<endl;
      out<<"ITEM: ATOMS index type x y z"<<endl;
      string name;
      for(int i=0; i<N; i++)
      {
            name="1";
            out<<setw(6)<<i+1<<"\t"<<name<<"\t"<<setw(8)<<P[i].coordinate.x<<"\t"<<setw(8)<<P[i].coordinate.y<<"\t"<<setw(8)<<P[i].coordinate.z<<endl;
      }
      out.close();
}
