// ParticleTools.h

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

        static TLorentzVector getMotherPtEtaPhiM(float pt0, float eta0, float phi0, float m0, float pt1, float eta1, float phi1, float m1);
        static _TrackInfo getGenMuDY(bool m, bool postFSR, VarSet& vars);
        static TLorentzVector convertTrackMuTo4Vec(_TrackInfo& track);
        static bool isValid4Vec(TLorentzVector& v);
        static float dR(float eta1, float phi1, float eta2, float phi2);
        static bool isMuPairGenMatchedDY(float dRcut, float ptCut, VarSet& vars);
         
};
#endif
