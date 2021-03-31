 #include "rdf.h" 

void RDF::RDF_ini(double L, int maxbin)
{
	delta_r = L/(double(maxbin));
	cout<<"delta_r = "<<delta_r<<"\t maxbin = "<<maxbin<<endl;
	int i;
	for(i = 0; i<=maxbin; i++)
	{
		hist.push_back(0);
		gr.push_back(0);
		r.push_back(0);
	}
}

//Inputs needed for this: L, x,y,z coord, N
//The rdf function being called again and again to calculate the net change in the particle moves
void RDF::RDF_r(double L, int N, const vector<Particle> P, int maxbin)
{
	int i, j, bin;
	double r2, r;
	XYZ counter;
	for(i=0; i<N-1; i++)
	{
		for(j = i+1; j<N-1; j++)
		{
			counter.x = P[i].coordinate.x - P[j].coordinate.x;
			counter.y = P[i].coordinate.y - P[j].coordinate.y;
			counter.z = P[i].coordinate.z - P[j].coordinate.z;
			mi.S.min_img(counter, L);
			r = counter.norm();
			bin = std::floor(r/delta_r);
			if(bin > maxbin || bin < 0)					//Loop needed to find out if there is a wrong bin assignemnt for any particle pair
			{
				cout<<"i = "<<i<<"\t X = "<<P[i].coordinate.x<<"\t Y = "<<P[i].coordinate.y<<"\t Z = "<<P[i].coordinate.z<<"\t VX = "<<P[i].velocity.vx<<"\t VY = "<<P[i].velocity.vy<<"\t VZ = "<<P[i].velocity.vz<<endl;
				cout<<"j = "<<j<<"\t X = "<<P[j].coordinate.x<<"\t Y = "<<P[j].coordinate.y<<"\t Z = "<<P[j].coordinate.z<<"\t VX = "<<P[j].velocity.vx<<"\t VY = "<<P[j].velocity.vy<<"\t VZ = "<<P[j].velocity.vz<<endl;
				cout<<"Bin getting wrong value in rdf_r function, bin = "<<bin<<"\t maxbin = "<<maxbin<<"\t r = "<<r<<"\t delta_r = "<<delta_r<<"\t rx = "<<counter.x<<"\t ry = "<<counter.y<<"\t rz = "<<counter.z<<endl;
				cout<<"The culprit particles are: i = "<<i<<"\t j = "<<j<<endl;
				exit(1);
			}
			hist[bin] = hist[bin] + 2;
			rdf_steps = rdf_steps + 1;
		}
	}
}	
//inputs are L, maxbin, N, a, temp, mass
//rdf funtion at the end to calculate the final normalised rdf
void RDF::RDF_write(double L, int N, int maxbin, std::string FileName)
{
	int i, bin;
	double rlower, rupper, nideal;
	for(i = 0; i<maxbin; i++)
	{
		rlower = double(i)*delta_r;
		rupper = double(i+1)*delta_r;
		nideal = 4*PI*double(N)*((rupper*rupper*rupper) - (rlower*rlower*rlower))/3;
		//nideal = 1;
		gr[i] = double(hist[i])/(double(rdf_steps)*double(N)*nideal);			//Normalising wrt nideal
		r[i] = rlower + 0.5*delta_r;
	}
	std::ostringstream fnd;//FileName dummy
	fnd<<"_rdf"<<FileName<<".dat";
	std::string FileName1 = fnd.str();

	ofstream out;
	out.open(FileName1,ios::app);
		out<<"Distance between particles \t RDF"<<endl;
		for(i=0; i<maxbin; i++)
		{
			out<<setw(8)<<r[i]<<"\t"<<setw(8)<<gr[i]<<endl;	
		}
	out.close();
}

