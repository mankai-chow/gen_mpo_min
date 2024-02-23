//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#pragma once

#include "itensor/mps/siteset.h"
#include "itensor/util/str.h"

using namespace tinyformat;

namespace itensor {

class FermionBilayerSite;

using FermionBilayer = BasicSiteSet<FermionBilayerSite>;

class FermionBilayerSite
    {
    Index s;
    public:

    FermionBilayerSite(Index I) : s(I) { }

    FermionBilayerSite(Args const& args = Args::global())
        {
        auto ts = TagSet("Site,BFermi");//,Bilayer
        auto n = 1;
        if(args.defined("SiteNumber"))
          {
          n = args.getInt("SiteNumber");
          ts.addTags("n="+str(n));
          }
        else
          {
            Error("SiteNumber is unknow!");
          }
        auto conserveQNs = args.getBool("ConserveQNs",true);
        auto conserve_Nf = args.getBool("ConserveNf", conserveQNs);
        auto conserve_Lz = args.getBool("ConserveLz", conserveQNs);
        //auto conserve_LzSz = args.getBool("ConserveLzSz", false);
        if(not conserveQNs)
            {
            s = Index(4,ts);
            }
        else if(conserve_Lz && conserve_Nf)
            {
                s = Index(QN({"Nf",0,-1}, {"Lz",0}),1,
                          QN({"Nf",1,-1}, {"Lz",n}),1,
                          QN({"Nf",1,-1}, {"Lz",n}),1,
                          QN({"Nf",2,-1}, {"Lz",2*n}),1,Out,ts);
            }
        else if( (!conserve_Lz) && conserve_Nf)
            {
                s = Index(QN({"Nf",0,-1}),1,
                          QN({"Nf",1,-1}),1,
                          QN({"Nf",1,-1}),1,
                          QN({"Nf",2,-1}),1,Out,ts);
            }
        /*
        if(conserve_LzSz)
            {
                s = Index(QN({"Nf",0,-1}, {"Sz",0}, {"Lz",0}),1,
                          QN({"Nf",1,-1}, {"Sz",+1}, {"Lz",n}),1,
                          QN({"Nf",1,-1}, {"Sz",-1}, {"Lz",n}),1,
                          QN({"Nf",2,-1}, {"Sz",0}, {"Lz",2*n}),1,Out,ts);
            }*/
        }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
        {
        if(state == "0" || state == "Emp") 
            {
            return s(1);
            }
        else 
        if(state == "+" || state == "Up") 
            {
            return s(2);
            }
        else 
        if(state == "-" || state == "Dn") 
            {
            return s(3);
            }
        else 
        if(state == "S" || state == "UpDn") 
            {
            return s(4);
            }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IndexVal{};
        }

	ITensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        IndexVal Em(s(1)),
                 EmP(sP(1)),
                 Up(s(2)),
                 UpP(sP(2)),
                 Dn(s(3)),
                 DnP(sP(3)),
                 UD(s(4)),
                 UDP(sP(4));

        ITensor Op(dag(s),sP);

        if(opname == "Nup")
            {
            Op.set(Up,UpP,1);
            Op.set(UD,UDP,1);
            }
        else
        if(opname == "Ndn")
            {
            Op.set(Dn,DnP,1);
            Op.set(UD,UDP,1);
            }
        else
        if(opname == "Nupdn")
            {
            Op.set(UD,UDP,1);
            }
        else
        if(opname == "Ntot" || opname == "N" )
            {
            Op.set(Up,UpP,1);
            Op.set(Dn,DnP,1);
            Op.set(UD,UDP,2);
            }
        else
        if(opname == "Cup")
            {
            Op.set(Up,EmP,1); 
            Op.set(UD,DnP,1); 
            }
        else
        if(opname == "Cdagup")
            {
            Op.set(Em,UpP,1); 
            Op.set(Dn,UDP,1);
            }
        else
        if(opname == "Cdn")
            {
            Op.set(Dn,EmP,1); 
            Op.set(UD,UpP,-1); 
            }
        else
        if(opname == "Cdagdn")
            {
            Op.set(Em,DnP,1); 
            Op.set(Up,UDP,-1);
            }
        else
        if(opname == "Aup")
            {
            Op.set(Up,EmP,1); 
            Op.set(UD,DnP,1); 
            }
        else
        if(opname == "Adagup")
            {
            Op.set(Em,UpP,1); 
            Op.set(Dn,UDP,1);
            }
        else
        if(opname == "Adn")
            {
            Op.set(Dn,EmP,1); 
            Op.set(UD,UpP,1); 
            }
        else
        if(opname == "Adagdn")
            {
            Op.set(Em,DnP,1); 
            Op.set(Up,UDP,1);
            }
        else
        if(opname == "FermiPhase" || opname == "F")
            {
            Op.set(Em,EmP,+1); 
            Op.set(Up,UpP,-1);
            Op.set(Dn,DnP,-1);
            Op.set(UD,UDP,+1);
            }
        else
        if(opname == "Fup")
            {
            Op.set(Em,EmP,+1); 
            Op.set(Up,UpP,-1);
            Op.set(Dn,DnP,+1);
            Op.set(UD,UDP,-1);
            }
        else
        if(opname == "Fdn")
            {
            Op.set(Em,EmP,+1); 
            Op.set(Up,UpP,+1);
            Op.set(Dn,DnP,-1);
            Op.set(UD,UDP,-1);
            }
        else
        if(opname == "Sz")
            {
            Op.set(Up,UpP,+0.5); 
            Op.set(Dn,DnP,-0.5);
            }
        else
        if(opname == "S+")
            {
            Op.set(Dn,UpP,1); 
            }
        else
        if(opname == "S-")
            {
            Op.set(Up,DnP,1); 
            }
        else
        if(opname == "S2")
            {
            //S dot S on-site
            Op.set(Up,UpP,0.75); 
            Op.set(Dn,DnP,0.75);
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }

    //
    // Deprecated, for backwards compatibility
    //

    FermionBilayerSite(int n, Args const& args = Args::global())
        {
        *this = FermionBilayerSite({args,"SiteNumber=",n});
        }

    };

//
// Deprecated, for backwards compatability
//

//using SpinlessSite = FermionFQHESite;

//using Spinless = BasicSiteSet<SpinlessSite>;


// order {1, up, dn, updn}

ITensor splitgate_electron(const Index& pIdx, const double& alpha, const double& beta)
{
    ITensor pi({sim(addTags(pIdx, "RegionA")), sim(addTags(pIdx, "RegionB")), dag(pIdx)});
    pi.set({1,1,1}, 1.0);       // |0>
    pi.set({2,1,2}, alpha);     // |up,A>
    pi.set({1,2,2}, beta);      // |up,B>
    pi.set({3,1,3}, alpha);     // |dn,A>
    pi.set({1,3,3}, beta);      // |dn,B>
    pi.set({4,1,4}, alpha*alpha);       // |updn,A>
    pi.set({1,4,4}, beta*beta);         // |updn,B>
    pi.set({2,3,4}, alpha*beta);
    pi.set({3,2,4}, -alpha*beta);
    return pi;
}

ITensor swapgate_electron(const Index& i1, const Index& i2)
{
    ITensor swapGate = ITensor({dag(i1), dag(i2), prime(i1), prime(i2)});

//  0
    int i = 1;
    // 0  0
    swapGate.set({i,1,i,1}, +1.0);
    // 0  up
    swapGate.set({i,2,i,2}, +1.0);
    // 0  dn
    swapGate.set({i,3,i,3}, +1.0);
    // 0  updn
    swapGate.set({i,4,i,4}, +1.0);

//  up
    i = 2;
    // up  0
    swapGate.set({i,1,i,1}, +1.0);
    // up  up
    swapGate.set({i,2,i,2}, -1.0);
    // up  dn
    swapGate.set({i,3,i,3}, -1.0);
    // up  updn
    swapGate.set({i,4,i,4}, +1.0);

//  dn
    i = 3;
    // dn  0
    swapGate.set({i,1,i,1}, +1.0);
    // dn  up
    swapGate.set({i,2,i,2}, -1.0);
    // dn  dn
    swapGate.set({i,3,i,3}, -1.0);
    // dn  updn
    swapGate.set({i,4,i,4}, +1.0);

//  updn
    i = 4;
    // updn  0
    swapGate.set({i,1,i,1}, +1.0);
    // updn  up
    swapGate.set({i,2,i,2}, +1.0);
    // updn  dn
    swapGate.set({i,3,i,3}, +1.0);
    // updn  updn
    swapGate.set({i,4,i,4}, +1.0);
    return swapGate;
}

} //namespace itensor
