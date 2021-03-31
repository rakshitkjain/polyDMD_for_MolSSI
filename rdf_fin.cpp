#include "header.h"

//this file is for rdf calculation in case the simulation needs to stop in between, so this will calculate the RDF from a coordinate file (DUMP) and use it to calculate rdfs

string DUMP;
int maxbin;
void ReadInput ( int argc, char *argv[] )
{
        options_description desc("Usage:\nDENSITY <options>");
	desc.add_options()
		("help,h", "print usage message")
		("DUMP,f", value<string>(&DUMP)->default_value("test.dat"), "dump file name")
		("maxbin,m", value<int>(&maxbin)->default_value(1000), "No. of bins for rdf calculation (default 1000)");
        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);
        if (vm.count("help"))
        {
          cout << desc << "\n";
          exit(1);
        }
        return;
}

void min_img(double &rx, double &ry, double &rz, double L)
{
	rx = rx -L*round(rx/L);
	ry = ry -L*round(ry/L);
	rz = rz -L*round(rz/L);
}

int main ( int argc, char *argv[] )
{
	cout<<"Reached";

	ReadInput(argc,argv);

	double L , delr;
	int N, n_sets;
	double dumpi, dumpd;
	char dumps[200];
	int count=0;//count gives the number of lines in the dump file
	vector<double> X;
	vector<double> Y;//storing particle coordinates
	vector<double> Z;
	double x,y,z, nideal, rlower, rupper;
    	int i, j, k, bin;
	double r2 = 1, r = 0, rx = 0, ry = 0, rz = 0;	

	vector<int>hist;//Used to tally the radial distribution function
	vector<double>gr;//Exact radial distribution function
	vector<double>rad;//Distance between two particles

	ifstream in;
	in.open(DUMP.c_str());
		while(!in.eof())
		{
			in.getline(dumps,256);
			count++;
		}
	in.close();

	cout<<"Number of lines in the dump file = "<<count<<endl;
	
	in.open(DUMP.c_str());    
		in.getline(dumps,256);
		in.getline(dumps,256);
		in.getline(dumps,256);
    
		in>>N;
    
		in.getline(dumps,256);
		in.getline(dumps,256);
		in>>dumpd>>L;
	in.close();
	n_sets = int(double(count)/double(N+9));//No. of sets stored in the dump file
	cout<<"N="<<N<<"\t L="<<L<<"\t Number of sets in the dump file = "<<n_sets<<endl;
	for(i = 0; i<=maxbin; i++)
	{
		hist.push_back(0);
		gr.push_back(0);
		rad.push_back(0);
	}

	delr = L/(2*double(maxbin));

/*	for(i = 0; i<N; i++)
	{
		X.push_back(0);
		Y.push_back(0);
		Z.push_back(0);
	}	*/
//Reading particle coordinates and storing in a single array//calculating hist value for the system
	int c=0;
	int sets = 0;
	in.open(DUMP.c_str());
		while(!in.eof())
		{				
			in.getline(dumps,256);
			in.getline(dumps,256);
			in.getline(dumps,256);
			in.getline(dumps,256);
			in.getline(dumps,256);
			in.getline(dumps,256);
			in.getline(dumps,256);
			in.getline(dumps,256);
			in.getline(dumps,256);
			c=c+9;
			while (c%(N+9) >= 9)
			{
				in>>dumpi>>dumpd>>x>>y>>z;
				c++;
				X.push_back(x);
				Y.push_back(y);
				Z.push_back(z);
			}
			in.getline(dumps,256);
/*			for(i=0; i<N; i++)
			{
				cout<<"i= "<<i<<"\t X = "<<X[sets*N + i]<<"\t Y= "<<Y[sets*N + i]<<"\t Z = "<<Z[sets*N + i]<<endl;
			}*/
			sets = sets + 1;
		}
	in.close();

	sets = 0;
	while(sets < n_sets)
	{
		for(i=0; i<N-1; i++)
		{
			for(j = i+1; j<N; j++)
			{
//				cout<<"i= "<<i<<"\t X = "<<X[sets*N + i]<<"\t Y= "<<Y[sets*N + i]<<"\t Z = "<<Z[sets*N + i]<<"j= "<<j<<"\t X = "<<X[sets*N + j]<<"\t Y= "<<Y[sets*N + j]<<"\t Z = "<<Z[sets*N + j]<<"\t rx = "<<rx<<"\t ry = "<<ry<<"\t rz = "<<rz<<"\t r2 = "<<r2<<"\t r ="<<r<<endl;
				rx = X[sets*N + i] - X[sets*N + j];
				ry = Y[sets*N + i] - Y[sets*N + j];
				rx = Z[sets*N + i] - Z[sets*N + j];
//				cout<<"i= "<<i<<"\t X = "<<X[sets*N + i]<<"\t Y= "<<Y[sets*N + i]<<"\t Z = "<<Z[sets*N + i]<<"j= "<<j<<"\t X = "<<X[sets*N + j]<<"\t Y= "<<Y[sets*N + j]<<"\t Z = "<<Z[sets*N + j]<<"\t rx = "<<rx<<"\t ry = "<<ry<<"\t rz = "<<rz<<"\t r2 = "<<r2<<"\t r ="<<r<<endl;
				min_img(rx,ry,rz, L);
				r2 = rx*rx + ry*ry + rz*rz;
				r = sqrt(r2);
//				cout<<"i= "<<i<<"\t X = "<<X[sets*N + i]<<"\t Y= "<<Y[sets*N + i]<<"\t Z = "<<Z[sets*N + i]<<"j= "<<j<<"\t X = "<<X[sets*N + j]<<"\t Y= "<<Y[sets*N + j]<<"\t Z = "<<Z[sets*N + j]<<"\t rx = "<<rx<<"\t ry = "<<ry<<"\t rz = "<<rz<<"\t r2 = "<<r2<<"\t r ="<<r<<endl;
/*				if(r<0.5)
				{
					cout<<"r<1 for set no. = "<<sets<<"\t and particles no. "<<i<<"\t and"<<j<<endl;
					exit(1);
				}*/
				bin = int(r/delr);
				hist[bin] = hist[bin] + 2;
//				r2 = 1, r = 0, rx = 0, ry = 0, rz = 0;
			}
		}
		sets = sets + 1;
	}
//Calculating the actual rdf
	for(i = 0; i<maxbin; i++)
	{
		rlower = double(i)*delr;
		rupper = double(i+1)*delr;
		nideal = 4*PI*double(N)*((rupper*rupper*rupper) - (rlower*rlower*rlower))/3;
		gr[i] = double(hist[i])/(double(n_sets)*double(N)*nideal);
		rad[i] = rlower + 0.5*delr;
	}
//Writing a file for the calculated RDF
	string name;
	name = "rdf_r"+DUMP+".txt";
	ofstream out;
	out.open(name);
	out<<"Distance between particles \t RDF"<<endl;
	for(i=0; i<maxbin; i++)
	{
		out<<setw(8)<<rad[i]<<"\t"<<setw(8)<<gr[i]<<endl;
	}
	out.close();
	return 0;
}
