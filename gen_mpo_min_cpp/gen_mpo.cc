#include "itensor/all.h"

#include <cmath>
#include <functional>
#include <set>
#include <algorithm>
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>

#include "site/fermionBilayer.h"
#include "interaction/cg.h"
#include "interaction/operator.h"

using namespace itensor;

int main(int argc, char *argv[])
{
    int N_o = 24;
    std::vector<double> pseudo = {4.75, 1.};

    auto sites = FermionBilayer(N_o, {"ConserveNf", true, "ConserveLz", true});
    MPO H;
    H = Hbulk_generator(sites, pseudo, 3.16, N_o, 12); 
}
