///////////////////////////////////////////////////////////////////////////
//                           ParticleTools.h                             //
//=======================================================================//
//                                                                       //
//        Miscellaneous tools for particles: output 4vec info to terminal//
//        get the correct gen muons from a gen event, some other tools   //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef PARTICLETOOLS
#define PARTICLETOOLS

#include "ParticleTools.h"
#include "TLorentzVector.h"
#include "VarSet.h"

class ParticleTools
{
    public: 
        ParticleTools(){};
        ~ParticleTools(){};

        static TString output4vecInfo(TLorentzVector& v);
        static TLorentzVector getMotherPtEtaPhiM(float pt0, float eta0, float phi0, float m0, float pt1, float eta1, float phi1, float m1);
        static bool isValid4Vec(TLorentzVector& v);
        static float dR(float eta1, float phi1, float eta2, float phi2);
};
#endif
