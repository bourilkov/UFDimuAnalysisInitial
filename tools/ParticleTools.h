// ParticleTools.h

#ifndef PARTICLETOOLS
#define PARTICLETOOLS

#include "ParticleTools.h"
#include "TLorentzVector.h"
#include "VarSet.h"

double static constexpr MASS_MUON = 0.105658367;

class ParticleTools
{
    public: 
        ParticleTools(){};
        ~ParticleTools(){};

        static TLorentzVector getMotherPtEtaPhiM(float pt0, float eta0, float phi0, float m0, float pt1, float eta1, float phi1, float m1);
        static _TrackInfo getGenMuDY(bool m, bool postFSR, VarSet& vars);
        static TLorentzVector convertTrackMuTo4Vec(_TrackInfo& track);
        static bool isValid4Vec(TLorentzVector& v);
};
#endif
