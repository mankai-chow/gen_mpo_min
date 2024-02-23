
#include "operator.h"
#include <chrono>
#include <iostream>
using namespace std;
using namespace std::chrono;
using namespace itensor;

MPO Hbulk_generator(SiteSet& sites, const std::vector<double>& pseudo, double t_hopping, int N_o, int num_thread)
{
    std::vector<double> Alist(N_o*N_o*N_o*N_o);
    cal2BodyAlist(pseudo, Alist, N_o, num_thread);

    auto ampo = AutoMPO(sites);
    for(auto m1=0;m1<N_o;++m1)
    {
        ampo +=  t_hopping, "Cdagup", m1+1, "Cdn", m1+1;
        ampo +=  t_hopping, "Cdagdn", m1+1, "Cup", m1+1;
        for(auto m2=0;m2<N_o;++m2)
        {
            for(auto m3=0;m3<N_o;++m3)
            {
                for(auto m4=0;m4<N_o;++m4)
                {
                    if(m1+m2!=m3+m4)
                        continue;
                    double cc = 0.0;
                    cc += 2*Alist[m1*N_o*N_o*N_o+m2*N_o*N_o+m3*N_o+m4];
                    ampo +=  cc,"Cdagup",m1+1,"Cdagdn",m2+1,"Cdn",m3+1,"Cup",m4+1;
                }
            }
        }
    }

    auto start = high_resolution_clock::now();
    MPO H = toMPO(ampo);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << duration.count() << endl;
    return H;
}

