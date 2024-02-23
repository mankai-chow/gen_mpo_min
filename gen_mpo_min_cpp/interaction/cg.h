
#ifndef __HEAD_CG_GENERAL
#define __HEAD_CG_GENERAL

#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>

int normCG2Body_j(std::vector<double>& cg, int j1, int j2, int j);
void calCG2Body_j(std::vector<double>& cg, int j1, int j2, int j);
double calcg(std::vector<double>& cg, int j1, int j2, int j, int x, int y);

class cgj1j2j
{
private:
    int j1_;
    int j2_;
    int j_;
    double s1_;
    double s2_;
    double s_;
    bool physicalQ_ = false;
    std::vector<double> cg_;
public:
    cgj1j2j(){};
    cgj1j2j(int j1, int j2, int j);
    ~cgj1j2j(){};

    const std::vector<double>& cg() const {return cg_;};
    const double& cgx(int x, int y) const {return cg_[x*j2_+y];};
    const double& cgm(double m1, double m2) const {return cg_[this->tox(m1)*j2_+this->toy(m2)];};

    friend std::ostream & operator<<(std::ostream & os, const cgj1j2j & obj);

    double tom1(int x) const {return (x-s1_);};
    double tom2(int y) const {return (y-s2_);};
    double tom( int z) const {return (z-s_ );};
    int tox(double m1) const {return round(m1+s1_);};
    int toy(double m2) const {return round(m2+s2_);};
    int toz(double m ) const {return round(m +s_ );};
};

class cgj1j2
{
private:
    int j1_;
    int j2_;
    double s1_;
    double s2_;
    double smin_;
    double smax_;
    std::vector<cgj1j2j> cglist_;
public:
    cgj1j2()=delete;
    cgj1j2(int j1, int j2);
    ~cgj1j2(){};

    const std::vector<cgj1j2j>& cglist(){return cglist_;};
    const cgj1j2j& cgj(int j){return this->cgs((double)(0.5*j-0.5));};
    const cgj1j2j& cgs(double s);

    friend std::ostream & operator<<(std::ostream & os, const cgj1j2 & obj);

    double tom1(int x){return (x-s1_);};
    double tom2(int y){return (y-s2_);};
    int tox(double m1){return round(m1+s1_);};
    int toy(double m2){return round(m2+s2_);};
};

void cal2BodyAlist(const std::vector<double>& pseudo, std::vector<double>& A_list2, int N_o, int num_thread);

void caldensityAlist(std::vector<double>& A_list2, int N_o, int l, int m, int num_thread);

void cal_deltadelta_lalb_Alist(std::vector<double>& A_list2, int l, int m, int la, int lb, int N_o, int num_thread);

void cal_deltadelta_lalb_Alist(double *A_list2, int l, int m, int la, int lb, int N_o, int num_thread);

double f0(double v0, double v1, int N_o, int l);
double f1(           double v1, int N_o, int l);

void calDefectAlist(const std::vector<double>& theta, 
                    const std::vector<double>& h, 
                    std::vector<double>& A_list2, int N_o);

#endif