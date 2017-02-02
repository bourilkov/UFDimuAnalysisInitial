//CollectionCleaner.hxx

#ifndef ADD_COLLECTIONCLEANER
#define ADD_COLLECTIONCLEANER

#include "VarSet.h"
#include "TLorentzVector.h"
#include "ParticleTools.h"
#include <vector>
#include <iostream>

class CollectionCleaner
{
    public:
        CollectionCleaner(){};
        ~CollectionCleaner(){};
        void cleanByDR(std::vector<TLorentzVector>& cleanThis, std::vector<TLorentzVector>& fromThis, float dRmin)
        {
        // remove items from cleanThis if they are too close in dR to any item in fromThis
            for(unsigned int i=0; i<fromThis.size(); i++)
            {   
                for(unsigned int j=0; j<cleanThis.size(); j++)
                {
                    if(cleanThis[j].DeltaR(fromThis[i]) < dRmin) 
                    {
                        // remove the item and make decrement j
                        // since the next item will fall back into
                        // the same index (j-- then j++ leaves j the same)
                        cleanThis.erase(cleanThis.begin()+j); 
                        j--;
                    }
                }
            }   

        };
};

#endif
