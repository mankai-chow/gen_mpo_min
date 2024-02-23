#ifndef __HEAD_OPERATOR_H
#define __HEAD_OPERATOR_H

#include "itensor/all.h"
#include "cg.h"

using namespace itensor;
MPO Hbulk_generator(SiteSet& sites, const std::vector<double>& pseudo, double t_hopping, int N_o, int num_thread=8);

#endif