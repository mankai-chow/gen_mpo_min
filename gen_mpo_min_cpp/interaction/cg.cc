
#include "cg.h"

int normCG2Body_j(std::vector<double>& cg, int j1, int j2, int j)
{
    double c, norm = 1.0;
    double s1 = 0.5*(j1-1);
    double s2 = 0.5*(j2-1);
    double s =  0.5*(j -1);

    int m1 = (j+j1-j2-1)/2;  // s+s1-s2
    int m2 = j2-1;           // 2*s2
    cg[m1*j2+m2] = 1.0;

    while( m1<j1-1&&m2>0 )
    {   
        c = -std::sqrt( (s1-(m1-s1) )*(s1+m1-s1+1)/( (s2-(m2-s2)+1)*(s2+m2-s2) ) );  //std::sqrt((j1-m1)*(j1+m1+1)/((j2-m2+1)*(j2+m2)));
        cg[(m1+1)*j2+m2-1] = c*cg[m1*j2+m2];
        norm += (cg[(m1+1)*j2+m2-1]*cg[(m1+1)*j2+m2-1]);
        ++m1;
        --m2;
    }
    norm = sqrt(norm);

    m1 = (j+j1-j2-1)/2;  // s+s1-s2
    m2 = j2-1;           // 2*s2
    cg[m1*j2+m2] /= norm;
    while( m1<j1-1&&m2>0 )
    {
        ++m1;
        --m2;
        cg[m1*j2+m2] /= norm;
    }
    return 0;
}

void calCG2Body_j(std::vector<double>& cg, int j1, int j2, int j)
{
    for(int x=0;x<j1;++x)// m1 = x - s1;
    {
        for(int y=0;y<j2;++y)// m2 = y - s2;
        {
            calcg(cg, j1, j2, j, x, y);
        }
    }
    return;
}

double calcg(std::vector<double>& cg, int j1, int j2, int j, int x, int y)
{
    double s1 = (j1-1)*0.5;
    double s2 = (j2-1)*0.5;
    double s  = (j -1)*0.5;
    double m1 = x - s1;
    double m2 = y - s2;
    if( ((fabs(m1)-s1)>0.1)||((fabs(m2)-s2)>0.1) )
        return 0.0;
    if(cg[x*j2+y]>-1.1)
        return cg[x*j2+y];
    if(fabs(m1+m2)>0.1+s)
    {
        cg[x*j2+y] = 0.0;
        return 0.0;
    }
    double c  = std::sqrt( (s+m1+m2+1)*(s-m1-m2) );
    double c1 = std::sqrt( (s1+m1+1)*(s1-m1) )/c;
    double c2 = std::sqrt( (s2+m2+1)*(s2-m2) )/c;
    cg[x*j2+y] = calcg(cg, j1, j2, j, x+1, y)*c1 + calcg(cg, j1, j2, j, x, y+1)*c2;
    return cg[x*j2+y];
}

/*class cgj1j2j
{
private:
    int j1_; // stored 2*j1+1
    int j2_; // stored 2*j2+1
    int j_; //  stored 2*j +1
    bool physicalQ_ = false;
    std::vector<double> cg_;
public:
    cgj1j2j()==delete;
    cgj1j2j(int j1, int j2, int j);
    ~cgj1j2j();

    const double& cg(int m1, int m2);
}*/
cgj1j2j::cgj1j2j(int j1, int j2, int j)
{
    j1_ = j1;
    j2_ = j2;
    if( j-1<std::abs(j1-j2) || j>j1+j2 ) // |2*j1-2*j2| = 2*j
        printf("Not physical!!(j1=%d, j2=%d, j=%d)\n", j1, j2, j);
    j_ = j;
    physicalQ_ = true;
    cg_ = std::vector<double>( j1*j2, -5.0 );
    normCG2Body_j(cg_, j1_, j2_, j_);
    calCG2Body_j(cg_, j1_, j2_, j_);
    s1_ = 0.5*(j1_-1);
    s2_ = 0.5*(j2_-1);
    s_  = 0.5*(j_ -1);
}

std::ostream & operator<<(std::ostream & os, const cgj1j2j & obj)
{
    os << std::setiosflags(std::ios::showpos);
    os << "\n************************************************************************************\n";
    os << "ClebschGordan Coefficients:\t" << std::setprecision(2) ;
    os << "j1 = " << obj.s1_ << ", j2 = " << obj.s2_ << ", j = " << obj.s_;
    os << "\n------------------------------------------------------------------------------------\n       | ";
    os << setiosflags(std::ios::fixed) << std::setprecision(10) ;
    for(int y=0;y<obj.j2_;++y)
        os << y-(0.5*obj.j2_-0.5) << "\t";
    os << "|\n------------------------------------------------------------------------------------\n";
    for(int x=0;x<obj.j1_;++x)
    {
        os << std::defaultfloat << setiosflags(std::ios::fixed) << std::setprecision(4) ;
        os << x-(0.5*obj.j1_-0.5) << "| ";
        os << std::defaultfloat << setiosflags(std::ios::scientific) << std::setprecision(6);
        for(int y=0;y<obj.j2_;++y)
            os << obj.cg_[x*obj.j2_+y] << "\t";
        os << "|\n";
    }
    os << "------------------------------------------------------------------------------------\n\n";
    os << std::defaultfloat;
    return os;
}

std::ostream & operator<<(std::ostream & os, const cgj1j2 & obj)
{
    os << std::defaultfloat << setiosflags(std::ios::fixed) << std::setprecision(4);
    os << "(s1="<<obj.s1_<<", s2="<<obj.s2_<<", smin="<<obj.smin_<<", smax"<<obj.smax_<<")\n" << std::defaultfloat;
    for(auto& i : obj.cglist_)
        os << i;
    return os;
}

cgj1j2::cgj1j2(int j1, int j2)
{
    j1_ = j1;
    j2_ = j2;
    s1_ = 0.5*(j1_-1);
    s2_ = 0.5*(j2_-1);
    smin_ = std::abs(s1_-s2_);
    smax_ = std::abs(s1_+s2_);
    cglist_ = std::vector<cgj1j2j>( round(smax_-smin_)+1 );
    for(double s=smin_;s<smax_+0.1;++s)
    {
        int ind = round(s-smin_);
        int j = round(2*s+1);
        cglist_.at(ind) = std::move( cgj1j2j( j1_, j2_, j) );
    }
}

const cgj1j2j&
cgj1j2::cgs(double s)
{
    if( s-smin_<-0.1 || s-smax_>0.1 )
    {
        printf("Invalid value (s1=+%lf, s2=+%lf, smin=+%lf, smax=+%lf, s=+%lf), exit!\n", s1_, s2_, smin_, smax_, s);
        exit(10);
    }
    return cglist_.at(round(s-smin_));
}

/*int main()
{
    auto cg = cgj1j2j(9, 6, 4);
    std::cout << cg;
    cg = cgj1j2j(6, 9, 13);
    std::cout << cg;
    auto cglist = cgj1j2(6, 9);
    std::cout << cglist;
    std::cout << cglist.cgs(5.5);
    std::cout << cglist.cgj(10);
    std::cout << cglist.cgj(10).cgm(0.5,-1.0);
    std::cout << round(-1.0+4) << std::endl;
    return 0;
}*/


void cal2BodyAlist(const std::vector<double>& pseudo, std::vector<double>& A_list2, int N_o, int num_thread)
{
    double s = 0.5*N_o-0.5;
    for(double& i : A_list2)
        i = 0;
    auto cg = std::move(cgj1j2(N_o, N_o));

    //printCG2(cg);
    //# pragma omp parallel for num_threads(num_thread)
    for(double m1=-s;m1<s+0.1;++m1)
        for(double m2=-s;m2<s+0.1;++m2)
            for(double m3=-s;m3<s+0.1;++m3)
                for(double m4=-s;m4<s+0.1;++m4)
                    for(int j=0;j<std::min(N_o, (int)pseudo.size());++j)
                        A_list2[(size_t)round(m1+s)*N_o*N_o*N_o + (size_t)round(m2+s)*N_o*N_o + (size_t)round(m3+s)*N_o + (size_t)round(m4+s)] += cg.cgs(2*s-j).cgm(m1,m2)*cg.cgs(2*s-j).cgm(m4,m3)*pseudo[j];
    //printf("%lf %lf\n",cg->cg[8][(u-1-5)*u+6],cg->cg[8][(u-1-6)*u+5]);
    return;
}

// pairing operator \Delta*\Delta
void cal_deltadelta_lalb_Alist(double *A_list2, int l, int m, int la, int lb, int N_o, int num_thread)
{
	double s = 0.5*N_o-0.5;
    for(auto i=0;i<N_o*N_o*N_o*N_o;++i)
    { A_list2[i] = 0;}
	/*// my suggestion
	double comm = (std::abs(la-lb)&1)?-1:1;
	comm /= std::sqrt(2*l+1);*/
	// Prof. He's suggestion
	double comm = (std::abs(la-lb+m)&1)?-1:1;
	comm /= std::sqrt(2*l+1);

	auto cglalbl = std::move(cgj1j2j( round(2*la+1), round(2*lb+1), round(2*l+1) ));
	auto cgssla = std::move(cgj1j2j(N_o, N_o, round(2*la+1) ));
	auto cgsslb = std::move(cgj1j2j(N_o, N_o, round(2*lb+1) ));

	//# pragma omp parallel for num_threads(num_thread)
	for(double j1=-s;j1<s+0.1;++j1)
	{
		for(double j2=-s;j2<s+0.1;++j2)
		{
			for(double j3=-s;j3<s+0.1;++j3)
			{
				for(double j4=-s;j4<s+0.1;++j4)
				{
					if( std::abs(j1+j2)>la+0.1 || std::abs(j3+j4)>lb+0.1 )
						continue;
					if( std::abs(j1+j2-j3-j4-m)>1e-14 )
						continue;
					A_list2[(int)round(j1+s)*N_o*N_o*N_o+(int)round(j2+s)*N_o*N_o+(int)round(j3+s)*N_o+(int)round(j4+s)] = 
					comm * cglalbl.cgm(j1+j2,-j3-j4) * cgssla.cgm(j1, j2) * cgsslb.cgm(j4, j3);
					// Prof. He's suggestion
					A_list2[(int)round(j1+s)*N_o*N_o*N_o+(int)round(j2+s)*N_o*N_o+(int)round(j3+s)*N_o+(int)round(j4+s)] *= 
					((int)round(j3+j4)&1)?-1:1;
				}
			}
		}		
	}
}

void cal_deltadelta_lalb_Alist(std::vector<double>& A_list2, int l, int m, int la, int lb, int N_o, int num_thread)
{
    cal_deltadelta_lalb_Alist(A_list2.data(), l, m, la, lb, N_o, num_thread);
}

double f0(double v0, double v1, int N_o, int l)
{
    double s = 0.5*N_o-0.5;
    auto cg = cgj1j2j(2*N_o-1, 2*N_o-1, 2*l+1 );
    double resu = 2.0*v0-0.5*l*(l+1) * (4*s-1)/(s*(4*s+1))*v1;
    resu *= std::pow(-1, N_o-1)*(4*s+1)*cg.cgm(-2*s, 2*s);
    return resu;
}
double f1(           double v1, int N_o, int l)
{
    double s = 0.5*N_o-0.5;
    auto cg = cgj1j2j(2*N_o-3, 2*N_o-3, 2*l+1 );
    double resu = 2.0*v1;
    resu *= std::pow(-1, N_o-2)*(4*s-1);
    //printf("%lf\n", resu);
    //printf("%d, %lf, %lf\n", 2*N_o-2, 2*s-1, cg.cgm(1-2*s,2*s-1));
    return resu*cg.cgm(1-2*s,2*s-1);
}

void calDefectAlist(const std::vector<double>& theta, 
                    const std::vector<double>& h, 
                    std::vector<double>& A_list2, int N_o)
{
    const double s = 0.5*(N_o-1);
    auto cg = cgj1j2(N_o, N_o ); // <s,s|l>
    
    if(theta.size()!=h.size() || theta.size()==0)
    {
        std::cout << "Warning: theta.size() = " << theta.size() << ", but h.size() = " << h.size() << std::endl;
        return;
    }

    std::vector<double> pl(N_o, 0.0);
    
    for(int l=0;l<N_o;++l)
        for(size_t i=0;i<theta.size();++i)
            pl.at(l) += h.at(i)*std::legendre(l, std::cos(theta.at(i)));

    for(int x1=0;x1<N_o;++x1)
    {
        double m1 = x1 - s;
        double temp = 0.0;
        for(int l=0;l<N_o;++l)
        {
            temp += pl.at(l)*std::pow(-1, l+N_o-1+x1)*cg.cgs(l).cgm(s,-s)*cg.cgs(l).cgm(m1,-m1);
        }
        A_list2.at(x1) = (s+0.5)*temp;
    }
}