//PUTools.h

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
