// BranchSet.h
// Branches for the ttree, so that we can load the data

#ifndef ADD_BRANCHSET
#define ADD_BRANCHSET

#include "TBranch.h"

class BranchSet
{
    public:
        BranchSet(){};
        ~BranchSet(){};

        TBranch* nVertices = 0;
        TBranch* eventInfo = 0;
        TBranch* mht = 0;

        TBranch* recoDimuCands = 0;
        TBranch* recoMuons = 0;
        TBranch* recoElectrons = 0;
        TBranch* jets = 0;

        TBranch* eff_wgt = 0;
        TBranch* pu_wgt = 0;
        TBranch* nPU = 0;
        TBranch* gen_wgt = 0;

        int getEntry(int i)
        {
            if(nVertices != 0) nVertices->GetEntry(i);
            if(eventInfo != 0) eventInfo->GetEntry(i);
            if(mht != 0) mht->GetEntry(i);

            if(recoDimuCands != 0) recoDimuCands->GetEntry(i);
            if(recoMuons != 0) recoMuons->GetEntry(i);
            if(recoElectrons != 0) recoElectrons->GetEntry(i);
            if(jets != 0) jets->GetEntry(i);

            if(eff_wgt != 0) eff_wgt->GetEntry(i);
            if(pu_wgt != 0) pu_wgt->GetEntry(i);
            if(nPU != 0) nPU->GetEntry(i);
            if(gen_wgt != 0) gen_wgt->GetEntry(i);
        }
};

#endif
