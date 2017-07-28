///////////////////////////////////////////////////////////////////////////
//                       CategorySelection.h                             //
//=======================================================================//
//                                                                       //
// Categorizer objects categorize each event. Define the cuts for the    //
// categories and the evaluate function to determine the category        //
// structure. We keep track of the different categories in the           //
// in the categorizer via categoryMap<TString, Category>. Each category  //
// tracks its historams with histoMap<TString, TH1D*> and some TLists.   //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ADD_CATEGORYSELECTION
#define ADD_CATEGORYSELECTION

#include "VarSet.h"
#include "TH1D.h"
#include "TList.h"
#include "TXMLEngine.h"
#include <map>
#include <utility>
#include <iostream>

//////////////////////////////////////////////////////////////////////////
//// ______________________Category_____________ _______________________//
//////////////////////////////////////////////////////////////////////////

class Category
{
   public:
       Category(){};
       ~Category(){};

       Category(TString key)
       {
           this->key = key;
           this->name = key;
           this->hide = false;
       }
    
       Category(TString key, bool hide)
       {
           this->key = key;
           this->name = key;
           this->hide = hide;
       }
    
       Category(TString key, bool hide, bool isTerminal)
       {
           this->key = key;
           this->name = key;
           this->hide = hide;
           this->isTerminal = isTerminal;
       }
    
       // map sample name to vector of events in the category
       std::map<TString, std::vector< std::pair<int, long long int> > > eventsMap;

       // Categorizer.evaluate will determine whether an event falls into this category
       // if an event falls into this category the boolean value will be set to true
       bool inCategory = false; 

       // is a final category used for limit setting
       bool isTerminal = false;

       // make a histogram of the category or not
       bool hide = false;

       // Usually we want to make a bunch of histograms for each category, so we keep them here
       std::map<TString, TH1D*> histoMap;   // access histos by name
       TList* histoList = new TList();      // histos ordered by xsec so that we can create the ratio and stack plot
       TList* signalList = new TList();     // signal histos
       TList* bkgList = new TList();        // bkg histos
       TList* dataList = new TList();       // data histo

       // Book-keeping
       TString key;
       TString name = "";
};

//////////////////////////////////////////////////////////////////////////
//// ______________________Categorizer__________ _______________________//
//////////////////////////////////////////////////////////////////////////

class Categorizer
{
// This class evaluates the event and determines which category or categories it belongs to

    public:
        // the categories the event may fall into
        std::map<TString, Category> categoryMap;

        // set up the categories and map them to a tstring
        virtual void initCategoryMap() = 0;

        // Determine which category the event belongs to
        virtual void evaluate(VarSet& vars) = 0;

        // reset the boolean values for the categories
        void reset()
        {
            for(auto &entry : categoryMap)
                entry.second.inCategory = false;
        };

        // output the category selection results
        void outputResults()
        {
            for(auto &entry : categoryMap)
                if(!entry.second.hide && entry.second.name != entry.second.key) 
                    std::cout << "    (" << entry.second.name << ") " << entry.first << ": " << entry.second.inCategory << std::endl;
                else if(!entry.second.hide)
                    std::cout << "    " << entry.first << ": " << entry.second.inCategory << std::endl;

            std::cout << std::endl;
        };
};

//////////////////////////////////////////////////////////////////////////
//// ______________________XMLCategorizer_______________________________//
//////////////////////////////////////////////////////////////////////////

// A decision node or terminal node for an XMLCategorizer that reads in an XML
// Decision Tree as the categorization.
class CategoryNode
{
    public: 
        CategoryNode(){};
        CategoryNode(CategoryNode* cmother, CategoryNode* cleft, CategoryNode* cright, 
                     TString ckey, double csplitVar, TString csplitVarName, double csplitVal, double csignificanceSquared)
        {
            mother = cmother;
            left = cleft;
            right = cright;
            key = ckey;
            name = key;
            splitVar = csplitVar;
            splitVarName = csplitVarName;
            splitVal = csplitVal;
            significanceSquared = csignificanceSquared;
        };
        ~CategoryNode(){};

        void theMiracleOfChildBirth();
        CategoryNode* filterEventToDaughter(VarSet& vars);

        void output()
        {
            std::cout << Form("/// %s \n  # splitVarName : %s \n  # splitVal     : %7.3f \n  # significance2: %5.3f \n\n", 
                              name.Data(), splitVarName.Data(), splitVal, significanceSquared);
        };

        CategoryNode* mother = 0;
        CategoryNode* left = 0;
        CategoryNode* right = 0;

        TString key;
        TString name;
        int splitVar;
        TString splitVarName;
        double splitVal;
        double significanceSquared;
};

// XMLCategorizer reads in an XML Decision Tree as the categorization.
class XMLCategorizer : public Categorizer
{

    public:
        XMLCategorizer();
        XMLCategorizer(TString xmlfile);
        ~XMLCategorizer(){};

        CategoryNode* rootNode = 0;
        void initCategoryMap();
        void evaluate(VarSet& vars);
        void evaluateRecursive(VarSet& vars, CategoryNode* cnode);
        void loadFromXML(TString filename);
        void loadFromXMLRecursive(TXMLEngine* xml, XMLNodePointer_t xnode, CategoryNode* cnode);
        CategoryNode* filterEvent(VarSet& vars);
        CategoryNode* filterEventRecursive(VarSet& vars);
};

//////////////////////////////////////////////////////////////////////////
//// ________XML Categorizer + Other Cuts_______________________________//
//////////////////////////////////////////////////////////////////////////

class CategorySelectionHybrid : public XMLCategorizer
{
// Should actually make a class for this in the .cxx file. 
// I define the functions here so that categorize.cxx doesn't break.
    public:
        CategorySelectionHybrid(){ initCategoryMap(); }; 
        CategorySelectionHybrid(TString xmlfile){ initCategoryMap(); loadFromXML(xmlfile); }; 

        // Determine which category the event belongs to
        void evaluate(VarSet& vars){ evaluateRecursive(vars, rootNode); };
        void initCategoryMap(){};
};

//////////////////////////////////////////////////////////////////////////
//// ________BDT Based Categories Run2__________________________________//
//////////////////////////////////////////////////////////////////////////

class CategorySelectionBDT : public Categorizer
{
    public:
        CategorySelectionBDT(); 

        // Determine which category the event belongs to
        void evaluate(VarSet& vars);
        void initCategoryMap();
};

//////////////////////////////////////////////////////////////////////////
//// ______________________Run1Categorizer______________________________//
//////////////////////////////////////////////////////////////////////////

class CategorySelectionRun1 : public Categorizer
{
// The run1 H->MuMu category selection
    public:
        CategorySelectionRun1(); 
        CategorySelectionRun1(float cLeadPtMin, float cSubleadPtMin, float cMETMax, float cDijetMassMinVBFT, 
                              float cDijetDeltaEtaMinVBFT, float cDijetMassMinGGFT,
                              float cDimuPtMinGGFT, float cDimuPtMin01T); 

        // Preselection
        float cLeadPtMin;
        float cSubleadPtMin;
        float cMETMax;

        // VBF Tight
        float cDijetMassMinVBFT;
        float cDijetDeltaEtaMinVBFT;
     
        // GGF Tight
        float cDijetMassMinGGFT;
        float cDimuPtMinGGFT;

        // 01Tight
        float cDimuPtMin01T; 

        // Determine which category the event belongs to
        // result stored in isVBFTight, isGGFTight, etc 
        void evaluate(VarSet& vars);
        void initCategoryMap();
};

//////////////////////////////////////////////////////////////////////////
//// ______________________SynchCategorizer______________________________//
//////////////////////////////////////////////////////////////////////////

class CategorySelectionSynch : public Categorizer
{
// Category selection for synchronization purposes
    public:
        CategorySelectionSynch(); 
        CategorySelectionSynch(float cLeadPtMin, float cSubleadPtMin, float cMETMax, float cDijetMassMinVBFT, 
                              float cDijetDeltaEtaMinVBFT, float cDijetMassMinGGFT,
                              float cDimuPtMinGGFT, float cDimuPtMin01T); 

        // Preselection
        float cLeadPtMin;
        float cSubleadPtMin;
        float cMETMax;

        // VBF Tight
        float cDijetMassMinVBFT;
        float cDijetDeltaEtaMinVBFT;
     
        // GGF Tight
        float cDijetMassMinGGFT;
        float cDimuPtMinGGFT;

        // 01Tight
        float cDimuPtMin01T; 

        // Determine which category the event belongs to
        // result stored in isVBFTight, isGGFTight, etc 
        void evaluate(VarSet& vars);
        void initCategoryMap();
};

//////////////////////////////////////////////////////////////////////////
//// ______________________Run2Categorizer______________________________//
//////////////////////////////////////////////////////////////////////////

class LotsOfCategoriesRun2 : public Categorizer
{
// Adrian's new categories for Run2, oh boy

    public:
        LotsOfCategoriesRun2(); 
        //LotsOfCategoriesRun2(); 

        // Preselection
        int  c_pre_numExtraLeptonsMax;
        bool c_pre_useTau;
        int  c_pre_numBJetsMin;

        // b-jet categories
        
        // left side of b-jet test, more b-jets
        // tth (2 extra leptons), tth/bbh (0 extra leptons)
        // catch_all_bucket (1 extra lepton or doesn't fall into other categories some other way)
        int c_1b_numExtraLeptons_tth; 
        int c_1b_numExtraLeptons_tth_bbh; 

            // ttH categories
            int c_tth_numExtraElectrons;

            // tth-bbh categories
            float c_tth_bbh_mbb;
            float c_tth_bbh_mt_bMET;
            float c_tth_bbh_MET;

        // right side of b-jet test, fewer b-jets
        // 0 extra leptons, 1 or 2 extra leptons
        int c_0b_numExtraLeptonsMin; 

            // nonVlH (0 extra leptons)
            int c_0b_nonVlH_njetsMin;
            
                // >= c_0b_nonVlH_njets
                float c_0b_nonVlH_2j_mjj_min_vbfTight;
                float c_0b_nonVlH_2j_dEtajj_min_vbfTight;

                float c_0b_nonVlH_2j_mjj_min_vbfLoose;
                float c_0b_nonVlH_2j_dEtajj_min_vbfLoose;

                float c_0b_nonVlH_2j_mjj_min_VhH;
                float c_0b_nonVlH_2j_mjj_max_VhH;
                float c_0b_nonVlH_2j_dEtajjMuMu_max_VhH;
            
                // < c_0b_nonVlH_njets
                float c_0b_nonVlH_01j_MET_min_ZvvH;              
                float c_0b_nonVlH_01j_dimuPt_min_gfTight;

                // muon geometry categorization for gf categories
                float c_geo_bmax;
                float c_geo_omax;
                float c_geo_emax;

            // VlH (1 or 2 extra leptons)
            float c_0b_VlH_MET_min;

                // VlH according to the V decays
                int c_0b_VlH_We_num_e;    
                int c_0b_VlH_We_num_mu;    

                int c_0b_VlH_Wmu_num_e;    
                int c_0b_VlH_Wmu_num_mu;    

                int c_0b_VlH_Ztautau_num_e;    
                int c_0b_VlH_Ztautau_num_mu;    

                int c_0b_VlH_Zmumu_num_e;
                int c_0b_VlH_Zmumu_num_mu;

                int c_0b_VlH_Zee_num_e;  
                int c_0b_VlH_Zee_num_mu;  

        // list of geometric categories
        std::vector<TString> geometricNames;

        // Determine which category the event belongs to
        void evaluate(VarSet& vars);
        void evaluateMuonGeometry(VarSet& vars);
        void initCategoryMap();
};

//////////////////////////////////////////////////////////////////////////
//// ______________________FEWZCategorizer______________________________//
//////////////////////////////////////////////////////////////////////////

class CategorySelectionFEWZ : public Categorizer
{
// Categories to compare DY, Data, and DY-FEWZ

    public:
        CategorySelectionFEWZ(); 
        CategorySelectionFEWZ(bool useRecoMu, bool useRecoJets); 
        CategorySelectionFEWZ(bool useReco, bool useRecoJets, float massSplit, float etaCentralSplit, float jetPtMin, float jetEtaMax); 

        bool useRecoMu;
        bool useRecoJets;

        // Selections
        float cMassSplit;
        float cEtaCentralSplit;
        float cJetPtMin;
        float cJetEtaMax;

        // Determine which category the event belongs to
        void evaluate(VarSet& vars);
        void initCategoryMap();
};


#endif
