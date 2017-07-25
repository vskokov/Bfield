
#include <iostream>
#include <random>
#include <numeric>
#include <vector>
#include <map>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/random.hpp>


const double Xsection=4.2;
double D=sqrt(Xsection/M_PI);

//double beta2;

unsigned seed = time(NULL)%1000;
boost::random::mt19937 generator(seed);

std::map<int,int> A;
std::map<int,double> a;
std::map<int,double> R;


void init_maps(void)
{
    // Gold
	A[79]=197;
 	a[79]=0.535;
	R[79]=6.38;

    // Cupper
	A[29]=63;
	a[29]=0.5977;
	R[29]=4.206;

    //Isobars
	A[40]=96;
	a[40]=0.48;
	R[40]=5.07;

	A[44]=96;
	a[44]=0.46;
	R[44]=5.14;
}



boost::random::uniform_real_distribution<double> unit_d(0.0,1.0);

class nucleon
{
	double x, y, z;
	int Q;

public:
	nucleon(double xc, double yc, double zc, int Qc)
	{
		x = xc;
		y = yc;
		z = zc;
		Q = Qc;
	}

	double get_x(void)
	{
		return x;
	}
	double get_y(void)
	{
		return y;
	}
	double get_z(void)
	{
		return z;
	}

	int get_Q(void)
	{
		return Q;
	}
};


class nucley
{
	int A, Z;
	double R, a, bx, beta2,  sqrts;
	std::vector<nucleon> nucleons;

	double BW(double r, double theta)
	{
		//double Rp = R*(1-beta2*0.25*1.261566261*(3.0*cos(theta)*cos(theta)-1)); //check it before using
        double Rp = R;
		return (1.0+exp(-R/a))/(1.0+exp((r-Rp)/a));
	}

	void random_rotation(double& vx, double& vy, double& vz)
	{
        // Use on your own risk
		double phi = 2.0*M_PI*unit_d(generator);
		double theta = acos(2.0*unit_d(generator)-1);
		double ux = sin(theta)*cos(phi);
		double uy = sin(theta)*sin(phi);
		double uz = cos(theta);
		double psi = 2.0*M_PI*unit_d(generator);
		double C = cos(psi);
		double S = sin(psi);
		double t = 1.0-cos(psi);

		double out_x = vx*(t*ux*ux+C)+
		               vy*(t*ux*uy-S*uz)+
		               vz*(t*ux*uz+S*uy);

		double out_y = vx*(t*ux*uy+S*uz)+
		               vy*(t*uy*uy+C)+
		               vz*(t*uy*uz-S*ux);

		double out_z = vx*(t*ux*uz-S*uy)+
		               vy*(t*uz*uy+S*ux)+
		               vz*(t*uz*uz+C);
		vx=out_x;
		vy=out_y;
		vz=out_z;
	}

	void sample(void)
	{
		int i=0;
		do
		{
			double x = 2.0*R*(unit_d(generator)-0.5)*2;
			double y = 2.0*R*(unit_d(generator)-0.5)*2;
			double z = 2.0*R*(unit_d(generator)-0.5)*2;
			double r = sqrt(x*x+y*y+z*z);
			double theta = atan2(y,x);

            
			double selector = unit_d(generator);
			if(selector<BW(r,theta))
			{
				//random_rotation(x,y,z); //check before using
				double gamma = sqrts/2.0;
				//accept
				if(i<Z)
				{
					nucleon TMP(x+bx,y,z/gamma,1);
					nucleons.push_back(TMP);
				}
				else
				{
					nucleon TMP(x+bx,y,z/gamma,0);
					nucleons.push_back(TMP);
				}
				i++;
				//std::cout << r<< "\n";
			}
		}
		while (i<A);
	}

public:
	nucley(int Ac, int Zc, double Rc, double ac, double bxc, double beta2, double sqrtsc)
	{
		A=Ac;
		Z=Zc;
		R=Rc;
		a=ac;
		bx=bxc;
		sqrts=sqrtsc;
		sample();
	}

	std::vector<double> position(int i)
	{
		std::vector<double> out;
		out.push_back(nucleons.at(i).get_x());
		out.push_back(nucleons.at(i).get_y());
		out.push_back(nucleons.at(i).get_z());
		return out;
	}

	std::vector<double> position_trans(int i)
	{
		std::vector<double> out;
		out.push_back(nucleons.at(i).get_x());
		out.push_back(nucleons.at(i).get_y());
		return out;
	}

	int charge(int i)
	{
		return(nucleons.at(i).get_Q());
	}

};


class collision
{
	nucley* A1;
	nucley* A2;
	int Ncoll;
	int Npart1;
	int Npart2;
	int A, Z;
	double sqrts,bx;
	std::vector<double> density;
	int den_N;
	double den_a;
	double xc;
	double yc;
	double step_a;

	void Ncoll_Npart_fun(void);
	void Nch_fun(void);

public:
	std::vector<double> center(void);
	int Nch;
	collision(int Ac, int Zc, double beta2, double sqrtsc)
	{
		A = Ac;
		Z = Zc;
		sqrts=sqrtsc;
		bx = 2*R[Z]*sqrt(unit_d(generator)); //impact parameter generation
		A1 = new nucley(A, Z, R[Z], a[Z], bx*0.5, beta2,  sqrts); //carefull with R and a
		A2 = new nucley(A, Z, R[Z], a[Z], -bx*0.5, beta2, sqrts);
		
		den_N=20; 
		step_a = 0.25;
		density.resize(den_N*den_N);
		
		for(int i=0;i<density.size();i++) density[i]=0.0;

		std::vector<double> c = center();
		xc = c[0];
		yc = c[1];
		
		Ncoll_Npart_fun();
		Nch_fun();
	}

	double int_to_x(int i)
	{
		return xc+step_a*(i-den_N/2); 
	}

	double int_to_y(int i)
	{
		return yc+step_a*(i-den_N/2); 
	}
	
	//bool sphaleron(std::vector<double> *pos, double *dilution);
	int Npart(void);

	void to_den(int i, int j, double v)
	{
		density[j*den_N+i] = v;
	}
	
	double from_den(int i, int j)
	{
		return density[j*den_N+i];
	}

	double eps2(void);
	double psi_RP(void);

	std::vector<double> MagneticField(double x, double y, double z);
	double get_b(void)
	{
		return bx;
	}
};

double distance_fun(std::vector<double> x, std::vector<double> y)
{
	//this is distance squared
	double sum=0;
	for(size_t i=0; i<x.size(); i++)
	{
		sum+=pow(x.at(i)-y.at(i),2);
	}
	return sum;
}

void collision::Ncoll_Npart_fun(void)
{
	Ncoll=0;
	Npart1=0;
	Npart2=0;
	std::vector<int> part1_hush((A));
	std::vector<int> part2_hush((A));

	for (int i=0; i<A; i++)
	{
		part1_hush.at(i) = 0 ;
		part2_hush.at(i) = 0 ;
	}

	for(int i1 = 0; i1 < A; i1++)
	{
		for(int i2 = 0; i2 < A; i2++)
		{
			std::vector<double> p1 = A1->position_trans(i1);
			std::vector<double> p2 = A2->position_trans(i2);
			double distance = distance_fun(p1,p2);
			if(distance<D)
			{
				//collision
				Ncoll++;
				part1_hush.at(i1) = 1;
				part2_hush.at(i2) = 1;

				double xb = 0.5*(p1[0]+p2[0]);
				double yb = 0.5*(p1[1]+p2[1]);

				int ix = den_N/2+int((xb-xc)/step_a+0.5); 
				int iy = den_N/2+int((yb-yc)/step_a+0.5); 
				if( (ix>-1)&&(ix<den_N) &&  (iy>-1)&&(iy<den_N) ) to_den(ix,iy, from_den(ix,iy)+2.0 );
			}
		}
	}


	Npart1 = std::accumulate(part1_hush.begin(), part1_hush.end(), 0);
	Npart2 = std::accumulate(part2_hush.begin(), part2_hush.end(), 0);
}


std::vector<double> collision::center(void)
{
	double x=0.0;
	double y=0.0;
	double N=0.0;
	for(int i1=0; i1 < A; i1++)
		for(int i2=0; i2 < A; i2++)
		{
			std::vector<double> p1 = A1->position_trans(i1);
			std::vector<double> p2 = A2->position_trans(i2);
			double distance = distance_fun(p1,p2);
			if(distance<D)
			{
				x+= p1.at(0)+p2.at(0);
				y+= p1.at(1)+p2.at(1);
				N+=2;
			}
		}
	x = x / N;
	y = y / N;
	
	std::vector<double> out;
	out.push_back(x);
	out.push_back(y);
	return out;
}



/*bool collision::sphaleron(std::vector<double> *pos, double *dilution)
{
	double max = 0;
	int im, jm;
	for(int i=0;i<den_N;i++)
	for(int j=0;j<den_N;j++)
	{
		double value = from_den(i,j);
		if(value>max) 
		{
			max = value;
			im = i;
			jm = j;
		}
	}
	
	double p = 500*unit_d(generator); 
	if(p<max) 
	{
		pos->at(0)=int_to_x(im)+0.5*step_a;
		pos->at(1)=int_to_y(jm)+0.5*step_a;
		*dilution = max/((double)(Npart1+Npart2));
		return true;
	}
                pos->at(0)=int_to_x(im)+0.5*step_a;
                pos->at(1)=int_to_y(jm)+0.5*step_a;
                *dilution = max/((double)(Npart1+Npart2));
	return false;
}
*/


double collision::eps2(void)
{
	double x=0.0;
	double y=0.0;
	double x2=0.0;
	double y2=0.0;
	double N=0.0;
	for(int i1=0; i1 < A; i1++)
		for(int i2=0; i2 < A; i2++)
		{
			std::vector<double> p1 = A1->position_trans(i1);
			std::vector<double> p2 = A2->position_trans(i2);
			double distance = distance_fun(p1,p2);
			if(distance<D)
			{
				x2 += p1.at(0)*p1.at(0) + p2.at(0)*p2.at(0);
				y2 += p1.at(1)*p1.at(1) + p2.at(1)*p2.at(1);
				x+= p1.at(0)+p2.at(0);
				y+= p1.at(1)+p2.at(1);
				N+=2;
			}
		}
	x = x / N;
	y = y / N;
	x2 = x2 / N;
	y2 = y2 / N;
	x2 = x2 - x*x;
	y2 = y2 - y*y;
	return (x2-y2)/(x2+y2);
}



double collision::psi_RP(void)
//Recation plane angle from ϵ_2
{
	double x=0.0;
	double y=0.0;
	double x2=0.0;
	double y2=0.0;
	double xy=0.0;
	double N=0.0;
	for(int i1=0; i1 < A; i1++)
		for(int i2=0; i2 < A; i2++)
		{
			std::vector<double> p1 = A1->position_trans(i1);
			std::vector<double> p2 = A2->position_trans(i2);
			double distance = distance_fun(p1,p2);
			if(distance<D)
			{
				x2 += p1.at(0)*p1.at(0) + p2.at(0)*p2.at(0);
				y2 += p1.at(1)*p1.at(1) + p2.at(1)*p2.at(1);
				xy+= p1.at(0)*p1.at(1) +  p2.at(0)*p2.at(1);
				x+= p1.at(0)+p2.at(0);
				y+= p1.at(1)+p2.at(1);
				N+=2;
			}
		}
	x = x / N;
	y = y / N;
	x2 = x2 / N;
	y2 = y2 / N;
	xy = xy / N;
	x2 = x2 - x*x;
	y2 = y2 - y*y;
	xy = xy - x*y; 
	
	double psi_RP = 0.5*(atan2(2*xy, x2-y2)+M_PI); 

	return psi_RP;
}



int collision::Npart(void)
{
        return (Npart1+Npart2)/2; 
}

void collision::Nch_fun(void)
{
	double x_hard = 0.13;
	double kval = 2.0;//from star fits to au+au
	double nbar = 0.9225*2.43*((1-x_hard) + 2.0*x_hard*Ncoll/(Npart1+Npart2));//matched to to Run 12 Au+Au
	double p_nbd = 1.0/(1.0+(nbar/kval));
	Nch = 0;
	boost::random::negative_binomial_distribution<int> NBD(kval,p_nbd);
	for(int i=0; i<(Npart1+Npart2)/2; i++)
	{
		Nch += int(NBD(generator)+0.5);
	}
}

std::vector<double> collision::MagneticField(double x, double y, double z)
{
	double Bx = 0.0;
	double By = 0.0;

	double v2 = 1.0 - pow(2.0/sqrts,2);
	double v = - sqrt(v2); //sign convention
	std::vector<double> O;
	O.push_back(x);
	O.push_back(y);
	O.push_back(z);


	double prefactor=1.0/137.0*(1.0-v2);
	for(int i=0; i<Z; i++)
	{
		std::vector<double> p1 = A1->position(i);
		std::vector<double> p2 = A2->position(i);

		double R2_1 = distance_fun(O,p1);
		double R2_2 = distance_fun(O,p2);

		p1.at(2) = 0.0;
		p2.at(2) = 0.0;
		O.at(2) = 0.0;

		double Rperp2_1 = distance_fun(O,p1);
		double Rperp2_2 = distance_fun(O,p2);

		double Rperpx_1  = O.at(0) - p1.at(0);
		double Rperpx_2  = O.at(0) - p2.at(0);

		double Rperpy_1  = O.at(1) - p1.at(1);
		double Rperpy_2  = O.at(1) - p2.at(1);

		if(R2_1>0.36)
		{
			Bx += 1.0/pow(sqrt(R2_1), 3.0)  * 1.0/pow(sqrt(1.0-Rperp2_1*v2/R2_1),3) * v * ( - Rperpy_1 );
			By += 1.0/pow(sqrt(R2_1), 3.0)  * 1.0/pow(sqrt(1.0-Rperp2_1*v2/R2_1),3) * v * (  Rperpx_1 );
		}

		if(R2_2>0.36)
		{
			Bx += 1.0/pow(sqrt(R2_2), 3.0)  * 1.0/pow(sqrt(1.0-Rperp2_2*v2/R2_2),3) * (-v) * ( - Rperpy_2 );
			By += 1.0/pow(sqrt(R2_2), 3.0)  * 1.0/pow(sqrt(1.0-Rperp2_2*v2/R2_2),3) * (-v) * (  Rperpx_2 );
		}
		//See equations (1) in https://arxiv.org/pdf/1111.1949.pdf
		//
	}

	Bx = Bx*prefactor;
	By = By*prefactor;

	std::vector<double> out;
	out.push_back(Bx);
	out.push_back(By);
	return out;
}


std::vector<double> Integrated_B(collision& BAM, std::vector<double> pos)
{
	int n=15;
	double sumx = 0.0;
	double sumy = 0.0;
	double L=0.25; 
	double step = L / n;
	for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
		{
			double x = pos.at(0)-L*0.5+i*step;
			double y = pos.at(1)-L*0.5+j*step;
			std::vector<double> B = BAM.MagneticField(x,y,0);
			sumx += B.at(0);
			sumy += B.at(1);
			B = BAM.MagneticField(x+step,y,0);
			sumx += B.at(0);
			sumy += B.at(1);
			B = BAM.MagneticField(x,y+step,0);
			sumx += B.at(0);
			sumy += B.at(1);
			B = BAM.MagneticField(x+step,y+step,0);
			sumx += B.at(0);
			sumy += B.at(1);
		}
	std::vector<double> Bint;
	Bint.push_back(0.25*sumx/n/n);
	Bint.push_back(0.25*sumy/n/n);
	return Bint;
}


int main(int)
{
	init_maps();
	int Z;
	
	std::vector<double> pos;
	pos.resize(2);
	double dilution;
    double translation_to_mpi2 = 2.044630589;
    int NoE=10000000;

	std::cin>>Z;
		
	for(int i=0; i<NoE; i++)
	{
		collision BAM(A[Z],Z, 0.0, 200.0);
		
		std::vector<double> center  = BAM.center();
		
        std::vector<double> Bint =  Integrated_B(BAM,pos);
		std::vector<double> Bcent =  BAM.MagneticField(0,0,0);
		std::cout <<    BAM.get_b() << " "; // impact parameter
		std::cout <<    translation_to_mpi2*Bcent.at(0) << " "; //B_x at the center
        std::cout <<    translation_to_mpi2*Bcent.at(1) << " "; //B_y at the center
		std::cout <<    translation_to_mpi2*Bint.at(0)  << " "; //B_x integrated 
		std::cout <<    translation_to_mpi2*Bint.at(1)  << " "; //B_y integrated 
		std::cout <<    -BAM.eps2() << " "; //epsilon_2
		std::cout <<    atan2(Bcent.at(1),Bcent.at(0))    << " "; //ψ_B defined by the field at the center 
		std::cout <<    atan2(Bint.at(1),Bint.at(0))    << " "; //ψ_B defined by the integrated field
		std::cout <<    BAM.psi_RP()    << " "; //ψ_2
		std::cout <<    BAM.Nch << " "; //NUmber of charged particles
        std::cout <<    center[0] << " "; // position of the center
        std::cout <<    center[1] << " "; // position of the center
        std::cout <<    BAM.Npart() << " "; //number of participants
        std::cout << "\n" << std::flush;
		
	}
}

