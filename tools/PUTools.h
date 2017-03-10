///////////////////////////////////////////////////////////////////////////
//                           PUTools.h                                   //
//=======================================================================//
//                                                                       //
//     Make histograms for the StandAloneLumiReweighter, which gives     //
//     us PU weights. But we usually just do this at analyzer level now, //
//     rather than use the stand alone tool to calculate the weight      //
//     at plotting time.                                                 //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef PUTOOLS
#define PUTOOLS

#include <map>
#include "TString.h"
#include "Sample.h"

class PUTools
{
    public:
        PUTools(){};
        ~PUTools(){};
        static void savePUHisto(Sample* s, TString savedir);
        static void makePUHistos(std::map<std::string, Sample*> samples);
};
#endif
